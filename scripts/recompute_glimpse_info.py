#!/usr/bin/env python3
"""
Recompute cohort-wide GLIMPSE2 INFO/AF and INFO/INFO on the FULL sample set.

Why this is needed
------------------
GLIMPSE2_phase computes INFO/AF and INFO/INFO from the target samples present
in *that particular run* (see GLIMPSE phase/src/io/genotype_writer.cpp). When a
large cohort is phased in sample batches and then recombined with `bcftools
merge`, those two INFO fields only reflect a single batch (typically the first
one), not the full cohort. `GLIMPSE2_ligate` merely re-aligns phase and copies
INFO verbatim, so the error persists into the final file.

The per-sample GT/DS/GP fields are correct regardless of batching (imputation is
independent per target sample). Only the cross-sample aggregate statistics need
to be recomputed, once, on the full cohort. This tool does exactly that.

Formula (identical to GLIMPSE, accumulated over ALL samples)
------------------------------------------------------------
For each biallelic site, using each diploid sample's posterior probabilities
gp0/gp1/gp2 (probability of carrying 0/1/2 ALT alleles):

    ds_i    = gp1_i + 2*gp2_i                          # per-sample ALT dosage E[g]
    ds_sum  = sum_i ds_i
    ds2_sum = sum_i ds_i^2      = sum_i (E[g_i])^2
    e2_sum  = sum_i (gp1_i + 4*gp2_i) = sum_i E[g_i^2]

    n_hap   = 2 * N_samples                             # all-diploid cohort (NIPT)
    AF      = ds_sum / n_hap
    INFO    = 1 - (e2_sum - ds2_sum) / (n_hap * AF * (1-AF))   for 0<AF<1, else 1
    INFO    = max(INFO, 0), rounded to `--decimals` (3 by default)

Note: (e2_sum - ds2_sum) == sum_i Var(g_i), and n_hap*AF*(1-AF) is the expected
binomial variance under HWE, so INFO is the IMPUTE2-style imputation R^2 score.
INFO/RAF (reference-panel frequency) is independent of target samples and is
left untouched.

Two sub-commands
----------------
recompute  (single-file method)
    Recompute AF/INFO from the FULL-cohort per-sample GP of one merged file, using
    the formula above. Correct for common variants, but note that GLIMPSE stores
    GP floored to 3 decimals (phase/src/io/genotype_writer.cpp). For very rare
    variants that quantisation biases the recovered INFO downward and can collapse
    a true 0.2 to 0. Use only when nothing but the final merged file is available.

aggregate  (per-batch method, accurate for rare variants)
    Combine the per-batch AF/INFO that GLIMPSE already computed at full internal
    precision, instead of re-deriving them from the lossy 3-decimal GP. Because
    the per-sample posterior-variance sum is additive across samples, the exact
    cohort statistics are recovered from each batch's written AF_b/INFO_b:

        n_hap   = sum_b (2 * N_b)                         # total haplotypes
        ds_sum  = sum_b AF_b * (2*N_b)                     # AF_b is float32 (near-exact)
        Var_b   = (1 - INFO_b) * (2*N_b) * AF_b*(1-AF_b)   # per-batch variance sum
        var_sum = sum_b Var_b
        AF      = ds_sum / n_hap
        INFO    = 1 - var_sum / (n_hap * AF*(1-AF))        for 0<AF<1, else 1

    This is exact up to INFO_b's 3-decimal rounding (~1e-3), so it stays accurate
    even in the rare-variant tail where `recompute` fails. It needs the per-batch
    files (identical site sets) plus one full-cohort file whose records and GP are
    kept; only that file's AF/INFO are overwritten.

Design
------
Both sub-commands use pysam with per-region multiprocessing (contig-level, or
position-range splitting derived from a lightweight POS-only scan for single-
chromosome / single-region inputs), following the same parallel pattern as
phase_graft_and_dosage.py. Each worker writes its region to a temporary file and
the temporaries are concatenated in genomic order.

Author: Shujia Huang

Examples
--------
    # Single-file recompute (final merged file only; rare-variant INFO approximate)
    python recompute_glimpse_info.py recompute -i final.raw.vcf.gz -o final.vcf.gz -t 16

    # Per-batch aggregate (accurate; run per chunk before the batch files are deleted)
    python recompute_glimpse_info.py aggregate \
        --records-from imputed_chunk_0.bcf \
        --batch-inputs imputed_chunk_0_b000.bcf imputed_chunk_0_b001.bcf \
        -o imputed_chunk_0.fixed.bcf -t 8
"""

import os
import sys
import atexit
import shutil
import logging
import argparse
import tempfile
import multiprocessing

import pysam


# Setup professional logging environment
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s [%(filename)s:%(lineno)d] - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger("GLIMPSE_INFO_Recompute")


def get_pysam_write_mode(filename: str) -> str:
    """
    Determine the optimal pysam VariantFile write mode based on file extension.
    Supports standard VCF (.vcf), compressed BGZF VCF (.vcf.gz), and binary BCF (.bcf).
    """
    fn_lower = filename.lower()
    if fn_lower.endswith('.bcf'):
        return 'wb'
    elif fn_lower.endswith('.vcf.gz'):
        return 'wz'
    elif fn_lower.endswith('.vcf'):
        return 'w'
    else:
        raise ValueError(
            f"Unsupported output extension: '{filename}'. "
            f"Use '.vcf', '.vcf.gz', or '.bcf'."
        )


def _ensure_info_headers(header, af_tag: str, info_tag: str) -> None:
    """
    Ensure the AF/INFO fields exist in the header (they normally do in GLIMPSE
    output). Adding is guarded so the shared header object stays untouched when
    the tags are already present, allowing records to be written without
    cross-header translation.
    """
    if af_tag not in header.info:
        header.info.add(af_tag, "A", "Float",
                        "ALT allele frequency computed from DS/GP across all target samples")
    if info_tag not in header.info:
        header.info.add(info_tag, "A", "Float",
                        "IMPUTE info quality score for diploid samples")


def _recompute_record(rec, n_hap: int, af_tag: str, info_tag: str, decimals: int) -> bool:
    """
    Recompute AF and INFO for a single biallelic record from its full-cohort GP,
    overwriting rec.info[af_tag] and rec.info[info_tag] in place.

    Returns True if the record was updated, False if it was left unchanged
    (non-biallelic sites are skipped, matching GLIMPSE which only emits INFO for
    n_allele == 2).
    """
    alts = rec.alts
    if not alts or len(alts) != 1:
        return False

    ds_sum = 0.0
    ds2_sum = 0.0
    e2_sum = 0.0

    # Single pass over all samples. Missing GP (e.g. a Ref/Ref default) is
    # treated as zero dosage, contributing nothing to any of the sums, exactly
    # as GLIMPSE handles its Ref/Ref default (gp0 = 1).
    for sample in rec.samples.values():
        gp = sample.get("GP")
        if not gp or len(gp) < 3:
            continue
        g1 = gp[1]
        g2 = gp[2]
        if g1 is None or g2 is None:
            continue
        ds = g1 + 2.0 * g2
        ds_sum += ds
        ds2_sum += ds * ds
        e2_sum += g1 + 4.0 * g2

    if n_hap <= 0:
        return False

    af = ds_sum / n_hap
    if 0.0 < af < 1.0:
        info = 1.0 - (e2_sum - ds2_sum) / (n_hap * af * (1.0 - af))
    else:
        info = 1.0
    if info < 0.0:
        info = 0.0

    # INFO fields are Number=A; write single-element tuples for the ALT allele.
    rec.info[af_tag] = (round(af, 6),)
    rec.info[info_tag] = (round(info, decimals),)
    return True


def _ensure_index(vcf_path: str) -> bool:
    """
    Ensure a tabix/CSI index exists for the input VCF/BCF file.
    Attempts auto-generation if missing; required for region-based parallel access.
    """
    for suffix in ('.tbi', '.csi'):
        if os.path.exists(vcf_path + suffix):
            return True

    logger.info(f"No index found for '{vcf_path}', attempting auto-generation ...")
    try:
        path_lower = vcf_path.lower()
        if path_lower.endswith('.bcf'):
            pysam.index(vcf_path)
        elif path_lower.endswith('.vcf.gz'):
            pysam.tabix_index(vcf_path, preset='vcf', force=True)
        else:
            logger.error(
                f"Cannot auto-index '{vcf_path}': only .vcf.gz and .bcf are supported. "
                f"Please bgzip-compress and index the file manually."
            )
            return False
        logger.info(f"Successfully created index for '{vcf_path}'")
        return True
    except Exception as e:
        logger.error(f"Failed to create index for '{vcf_path}': {e}")
        return False


def _chunk_contigs(contigs, num_workers):
    """
    Distribute contigs into approximately equal *contiguous* blocks so that
    merging temporary outputs in worker order preserves genomic sort order.
    """
    n = len(contigs)
    if n <= num_workers:
        return [[c] for c in contigs]

    base = n // num_workers
    remainder = n % num_workers
    chunks = []
    start = 0
    for i in range(num_workers):
        size = base + (1 if i < remainder else 0)
        if size > 0:
            chunks.append(contigs[start:start + size])
            start += size
    return [c for c in chunks if c]


def _find_populated_contigs(vcf_path, candidate_contigs):
    """
    Return the subset of candidate_contigs that actually contain variant records.
    Uses index-backed random access so this is O(num_active_contigs), not
    O(file_size). Contigs listed in the header but absent from the body (common
    in single-chromosome extracts) are excluded.
    """
    populated = []
    with pysam.VariantFile(vcf_path) as vf:
        for contig in candidate_contigs:
            try:
                it = vf.fetch(contig)
            except (ValueError, OSError):
                logger.debug(f"Contig '{contig}' not accessible via index; skipping.")
                continue
            try:
                next(it)
                populated.append(contig)
            except StopIteration:
                pass  # contig in index but region is empty
    return populated


def _scan_contig_spans(vcf_path, contigs):
    """
    Single lightweight pass over the file to obtain, for each populated contig,
    the actual [min_pos, max_pos, record_count]. Only CHROM/POS are touched; the
    per-sample FORMAT fields are never decoded, so this pass is far cheaper than
    the recompute pass itself.

    This makes position-range splitting independent of whether the header declares
    contig lengths (GLIMPSE phase output without --contigs-fai carries none), and
    it bounds every range to where records actually are, so a single-region input
    is split across workers without generating empty ranges.
    """
    target = set(contigs)
    spans = {}  # contig -> [min_pos, max_pos, count]
    with pysam.VariantFile(vcf_path) as vf:
        for rec in vf:
            c = rec.contig
            if c not in target:
                continue
            p = rec.pos  # 1-based
            s = spans.get(c)
            if s is None:
                spans[c] = [p, p, 1]
            else:
                if p < s[0]:
                    s[0] = p
                if p > s[1]:
                    s[1] = p
                s[2] += 1
    return spans


def _split_span(contig, lo, hi, n):
    """
    Split the 1-based inclusive interval [lo, hi] into n contiguous, non-overlapping
    sub-ranges with no gaps. Returns a list of (contig, start, end) tuples. n is
    clamped so a range is never finer than 1 bp.
    """
    span = hi - lo + 1
    n = max(1, min(n, span))
    if n == 1:
        return [(contig, lo, hi)]
    base = span // n
    rem = span % n
    ranges = []
    start = lo
    for j in range(n):
        size = base + (1 if j < rem else 0)
        end = start + size - 1
        ranges.append((contig, start, end))
        start = end + 1
    return ranges


def _allocate_workers(counts, num_workers):
    """
    Distribute num_workers slots across contigs proportionally to their record
    counts, giving every populated contig at least one slot (largest-remainder
    method). Called only when len(counts) < num_workers, so the allocation can
    always sum up to num_workers.
    """
    k = len(counts)
    alloc = [1] * k
    remaining = num_workers - k
    total = sum(counts)
    if remaining <= 0 or total <= 0:
        return alloc
    quotas = [remaining * cnt / total for cnt in counts]
    base = [int(q) for q in quotas]
    for i in range(k):
        alloc[i] += base[i]
    leftover = remaining - sum(base)
    order = sorted(range(k), key=lambda i: quotas[i] - base[i], reverse=True)
    for i in range(leftover):
        alloc[order[i]] += 1
    return alloc


def _build_region_chunks(vcf_path, num_workers):
    """
    Build parallel work chunks. Each chunk is a list of items that are EITHER:
      - a contig name (str)             -> process the entire contig
      - a (contig, start, end) tuple    -> process a 1-based inclusive range

    Splitting is driven purely by which contigs actually contain records and by
    their real position spans (never by header contig lengths):
      * >= num_workers populated contigs -> split at contig granularity.
      * fewer contigs than workers (incl. single-chromosome / single-region) ->
        scan the real (min, max, count) per contig and split each contig's span
        by position, so every worker receives a non-empty slice.
    """
    with pysam.VariantFile(vcf_path) as vf:
        header_contigs = list(vf.header.contigs.keys())

    # _find_populated_contigs preserves header order, which is the genomic sort
    # order, so temp files merged in chunk order stay correctly sorted.
    contigs = _find_populated_contigs(vcf_path, header_contigs)

    if not contigs:
        logger.warning("No contigs with variant records found; nothing to parallelise.")
        return []

    # Enough populated contigs: split at contig granularity (no scan needed).
    if len(contigs) >= num_workers:
        return _chunk_contigs(contigs, num_workers)

    # Fewer contigs than workers: derive real spans and split by position so that
    # single-chromosome and single-region inputs still use every worker.
    spans = _scan_contig_spans(vcf_path, contigs)
    contigs = [c for c in contigs if c in spans and spans[c][2] > 0]
    if not contigs:
        logger.warning("No records found while scanning spans; nothing to parallelise.")
        return []

    counts = [spans[c][2] for c in contigs]
    alloc = _allocate_workers(counts, num_workers)

    chunks = []
    for contig, n_slots in zip(contigs, alloc):
        lo, hi, cnt = spans[contig]
        n_slots = min(n_slots, cnt)  # never create more ranges than records
        for rng in _split_span(contig, lo, hi, n_slots):
            chunks.append([rng])

    while len(chunks) > num_workers:
        chunks[-2].extend(chunks[-1])
        chunks.pop()

    return [c for c in chunks if c]


def _process_region_worker(args):
    """
    Multiprocessing worker: recompute AF/INFO for an assigned set of regions and
    write the results to a temporary output file. Runs in a child process.
    """
    input_path, temp_output_path, regions, write_mode, af_tag, info_tag, decimals = args

    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s [Worker-%(process)d] - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    worker_logger = logging.getLogger("GLIMPSE_INFO_Recompute.Worker")

    processed = 0
    updated = 0
    try:
        with pysam.VariantFile(input_path) as vcf_in:
            header = vcf_in.header  # shared object: records write without translation
            _ensure_info_headers(header, af_tag, info_tag)
            n_hap = 2 * len(header.samples)

            with pysam.VariantFile(temp_output_path, write_mode, header=header) as vcf_out:
                for item in regions:
                    if isinstance(item, tuple):
                        contig, start, end = item
                        label = f"{contig}:{start}-{end}"
                    else:
                        contig, start, end = item, None, None
                        label = contig

                    try:
                        if start is not None and end is not None:
                            # pysam.fetch() uses 0-based half-open coordinates;
                            # our ranges are 1-based inclusive, so subtract 1 from
                            # start to avoid dropping the first variant of a chunk.
                            it = vcf_in.fetch(contig, start - 1, end)
                        else:
                            it = vcf_in.fetch(contig)
                    except (ValueError, OSError):
                        worker_logger.debug(f"Region '{label}' not found in index; skipping.")
                        continue

                    for rec in it:
                        # Boundary guard: tabix bins may return records just
                        # outside the assigned range; skip them to prevent
                        # duplicate output across adjacent chunks.
                        if start is not None and (rec.pos < start or rec.pos > end):
                            continue
                        if _recompute_record(rec, n_hap, af_tag, info_tag, decimals):
                            updated += 1
                        vcf_out.write(rec)
                        processed += 1

        worker_logger.info(f"Completed region '{regions[0] if len(regions) == 1 else '%d regions' % len(regions)}': "
                           f"{processed} variants written, {updated} AF/INFO recomputed.")
    except Exception:
        worker_logger.exception("Worker encountered a fatal error.")
        raise

    return processed


def _merge_temp_files(temp_files, output_path, write_mode):
    """
    Concatenate temporary VCF/BCF chunks into the final output, in order.
    Writes the header once and streams records from each temp file.
    """
    logger.info(f"Merging {len(temp_files)} temporary chunk(s) into '{output_path}' ...")
    if not temp_files:
        logger.error("No temporary files to merge.")
        return

    with pysam.VariantFile(temp_files[0]) as first:
        header = first.header

    merged = 0
    with pysam.VariantFile(output_path, write_mode, header=header) as out_vcf:
        for temp_file in temp_files:
            with pysam.VariantFile(temp_file) as tf:
                for rec in tf:
                    out_vcf.write(rec)
                    merged += 1
    logger.info(f"Merge complete: {merged} total variants written to output.")


def _recompute_parallel(input_in, output_out, write_mode, num_workers,
                        af_tag, info_tag, decimals, tmp_dir):
    """
    Orchestrate parallel recomputation: ensure index, build region chunks,
    process in parallel, then merge temporary outputs in genomic order.
    """
    logger.info(f"Launching parallel recomputation with {num_workers} workers ...")

    if not _ensure_index(input_in):
        logger.error("Cannot proceed in parallel: input index is missing.")
        logging.shutdown()
        sys.exit(1)

    chunks = _build_region_chunks(input_in, num_workers)
    actual_workers = len(chunks)
    if actual_workers == 0:
        logger.warning("No contigs with variant records found; writing empty output.")
        with pysam.VariantFile(input_in) as vf:
            header = vf.header
        with pysam.VariantFile(output_out, write_mode, header=header):
            pass
        return

    n_contig_chunks = sum(1 for c in chunks if all(isinstance(x, str) for x in c))
    n_region_chunks = actual_workers - n_contig_chunks
    logger.info(
        f"Work distributed across {actual_workers} worker(s): "
        f"{n_contig_chunks} contig-level, {n_region_chunks} position-range chunk(s)."
    )

    # Keep temporaries beside the output to avoid tmpfs size limits on HPC nodes.
    base_dir = tmp_dir if tmp_dir else os.path.dirname(os.path.abspath(output_out))
    temp_dir = tempfile.mkdtemp(prefix="glimpse_info_", dir=base_dir)
    atexit.register(shutil.rmtree, temp_dir, ignore_errors=True)

    if output_out.endswith('.bcf'):
        temp_ext = '.bcf'
    elif output_out.endswith('.vcf.gz'):
        temp_ext = '.vcf.gz'
    else:
        temp_ext = '.vcf'

    worker_args = []
    for i, chunk in enumerate(chunks):
        temp_path = os.path.join(temp_dir, f"chunk_{i:04d}{temp_ext}")
        worker_args.append((input_in, temp_path, chunk, write_mode, af_tag, info_tag, decimals))

    try:
        with multiprocessing.Pool(processes=actual_workers) as pool:
            results = pool.map(_process_region_worker, worker_args)
        total = sum(r for r in results if isinstance(r, int))
        logger.info(f"All workers finished. Total variants processed: {total}")
    except Exception:
        logger.exception("Parallel processing aborted due to a worker failure.")
        logging.shutdown()
        sys.exit(1)

    temp_files = [os.path.join(temp_dir, f"chunk_{i:04d}{temp_ext}") for i in range(len(chunks))]
    temp_files = [f for f in temp_files if os.path.exists(f)]
    _merge_temp_files(temp_files, output_out, write_mode)

    shutil.rmtree(temp_dir, ignore_errors=True)
    logger.info("Parallel recomputation finished successfully.")


def _recompute_sequential(input_in, output_out, write_mode, af_tag, info_tag, decimals):
    """Single-process streaming recomputation (used when threads == 1)."""
    with pysam.VariantFile(input_in) as vcf_in:
        header = vcf_in.header
        _ensure_info_headers(header, af_tag, info_tag)
        n_hap = 2 * len(header.samples)

        processed = 0
        updated = 0
        with pysam.VariantFile(output_out, write_mode, header=header) as vcf_out:
            for rec in vcf_in:
                if _recompute_record(rec, n_hap, af_tag, info_tag, decimals):
                    updated += 1
                vcf_out.write(rec)
                processed += 1
                if processed % 50000 == 0:
                    logger.info(f"Processed {processed} variants ...")
    logger.info(f"Sequential run finished: {processed} variants written, {updated} AF/INFO recomputed.")


def recompute_info(input_in: str, output_out: str, threads: int,
                   af_tag: str, info_tag: str, decimals: int, tmp_dir: str) -> None:
    """
    Recompute cohort-wide INFO/AF and INFO/INFO. Dispatches to the parallel path
    when threads > 1, otherwise streams sequentially.
    """
    write_mode = get_pysam_write_mode(output_out)
    logger.info(f"Input : {input_in}")
    logger.info(f"Output: {output_out} (mode: '{write_mode}')")

    with pysam.VariantFile(input_in) as vf:
        n_samples = len(vf.header.samples)
    if n_samples == 0:
        logger.error("Input contains no samples; cannot recompute cohort statistics.")
        logging.shutdown()
        sys.exit(1)
    logger.info(f"Cohort: {n_samples} samples, n_hap = {2 * n_samples} (assumed all-diploid).")

    if threads > 1:
        _recompute_parallel(input_in, output_out, write_mode, threads,
                            af_tag, info_tag, decimals, tmp_dir)
    else:
        _recompute_sequential(input_in, output_out, write_mode, af_tag, info_tag, decimals)


# =============================================================================
# Sub-command 2: aggregate
# -----------------------------------------------------------------------------
# Instead of re-deriving AF/INFO from the merged file's lossy 3-decimal GP
# (the `recompute` path above), combine the per-batch AF_b/INFO_b that GLIMPSE
# already wrote at full internal precision. See the module docstring for the
# derivation. This is the accurate path for rare variants and is meant to run
# per chunk, right after `bcftools merge` and BEFORE the per-batch files are
# deleted, so their AF/INFO are still available.
# =============================================================================

def _get_info_scalar(rec, tag):
    """
    Read a single INFO value. GLIMPSE's AF/INFO are Number=A, which pysam returns
    as a 1-tuple; this normalises both tuple and plain-scalar returns to a float,
    or None when the tag (or the record) is absent.
    """
    if rec is None:
        return None
    val = rec.info.get(tag)
    if val is None:
        return None
    if isinstance(val, (tuple, list)):
        return val[0] if len(val) > 0 else None
    return val


def _aggregate_record(rec, batch_recs, n_haps, n_hap_total, af_tag, info_tag, decimals):
    """
    Overwrite rec.info[af_tag]/rec.info[info_tag] with the exact cohort AF/INFO
    reconstructed from the per-batch AF_b/INFO_b. Returns True on success, or
    False if a batch value is unusable (then the record keeps GLIMPSE's value).

    Key identity: the per-sample posterior-variance sum is additive across
    samples, so summing each batch's variance term Var_b gives the exact cohort
    numerator. AF_b is stored by GLIMPSE as float32 (no rounding), so ds_sum is
    near-exact; the only loss is INFO_b's 3-decimal rounding inside Var_b.
    """
    ds_sum = 0.0
    var_sum = 0.0
    for brec, nhap in zip(batch_recs, n_haps):
        af_b = _get_info_scalar(brec, af_tag)
        if af_b is None:
            # Without this batch's AF the cohort dosage mass is incomplete -> abort.
            return False
        # ds_sum_b = AF_b * (2*N_b): recover the batch dosage sum.
        ds_sum += af_b * nhap
        # Var_b is non-zero only for a polymorphic batch; GLIMPSE emits INFO_b = 1
        # (=> Var_b = 0) whenever AF_b is exactly 0 or 1, so we skip those.
        if 0.0 < af_b < 1.0:
            info_b = _get_info_scalar(brec, info_tag)
            if info_b is None:
                info_b = 1.0
            # Invert GLIMPSE's definition INFO_b = 1 - Var_b/(nhap_b*AF_b*(1-AF_b)).
            var_sum += (1.0 - info_b) * nhap * af_b * (1.0 - af_b)

    if n_hap_total <= 0:
        return False

    af = ds_sum / n_hap_total
    if 0.0 < af < 1.0:
        info = 1.0 - var_sum / (n_hap_total * af * (1.0 - af))
    else:
        info = 1.0
    if info < 0.0:
        info = 0.0

    # INFO fields are Number=A; write single-element tuples for the ALT allele.
    rec.info[af_tag] = (round(af, 6),)
    rec.info[info_tag] = (round(info, decimals),)
    return True


def _consume_batch_records_at(pos, ref, alts, batch_iters, batch_cur):
    """
    Advance each per-batch iterator to `pos` and, if its head record matches the
    site (same POS/REF/ALT), consume and return it. `batch_cur` is the per-iterator
    look-ahead buffer and is mutated in place.

    Returns (batch_records, all_matched): batch_records[k] is batch k's matching
    record (or None), and all_matched is True only when every batch supplied one.
    Since all batches are imputed against the same reference chunk their site sets
    are identical, so all_matched is True in the normal case; the per-position
    advance also keeps things aligned across boundary-skipped or multiallelic rows.
    """
    brecs = []
    all_matched = True
    for k in range(len(batch_iters)):
        cur = batch_cur[k]
        # Skip any batch records preceding pos (e.g. tabix bin spill at chunk edges,
        # or a main record that was boundary-skipped just before this one).
        while cur is not None and cur.pos < pos:
            cur = next(batch_iters[k], None)
        if (cur is not None and cur.pos == pos
                and cur.ref == ref and cur.alts == alts):
            brecs.append(cur)
            batch_cur[k] = next(batch_iters[k], None)  # consume the matched record
        else:
            brecs.append(None)
            batch_cur[k] = cur                         # keep head for the next site
            all_matched = False
    return brecs, all_matched


def _aggregate_region_worker(args):
    """
    Multiprocessing worker for `aggregate`.

    For its assigned region(s) it walks the full-cohort records-from file and, in
    lockstep, the matching records of every per-batch file, then writes each
    records-from record with AF/INFO replaced by the cohort aggregate. Runs in a
    child process and writes to its own temporary file.
    """
    (records_from_path, batch_paths, temp_output_path, regions, write_mode,
     n_haps, af_tag, info_tag, decimals) = args

    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s [Worker-%(process)d] - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    worker_logger = logging.getLogger("GLIMPSE_INFO_Recompute.Aggregate.Worker")

    n_hap_total = sum(n_haps)
    processed = 0
    updated = 0
    try:
        main_vcf = pysam.VariantFile(records_from_path)
        batch_vcfs = [pysam.VariantFile(p) for p in batch_paths]
        try:
            header = main_vcf.header  # shared object -> records write without translation
            _ensure_info_headers(header, af_tag, info_tag)

            with pysam.VariantFile(temp_output_path, write_mode, header=header) as vcf_out:
                for item in regions:
                    if isinstance(item, tuple):
                        contig, start, end = item
                    else:
                        contig, start, end = item, None, None

                    # One iterator per file over the SAME region. pysam.fetch() uses
                    # 0-based half-open coordinates; our ranges are 1-based inclusive,
                    # hence start-1.
                    try:
                        if start is not None:
                            main_it = main_vcf.fetch(contig, start - 1, end)
                            batch_iters = [bv.fetch(contig, start - 1, end) for bv in batch_vcfs]
                        else:
                            main_it = main_vcf.fetch(contig)
                            batch_iters = [bv.fetch(contig) for bv in batch_vcfs]
                    except (ValueError, OSError):
                        worker_logger.debug(f"Region '{contig}' not accessible via index; skipping.")
                        continue

                    # Look-ahead buffer holding each batch iterator's current head.
                    batch_cur = [next(it, None) for it in batch_iters]

                    for rec in main_it:
                        pos = rec.pos
                        # Boundary guard: drop records outside the assigned range so
                        # adjacent chunks never emit the same site twice. The next
                        # in-range record re-syncs the batch iterators via advance.
                        if start is not None and (pos < start or pos > end):
                            continue

                        alts = rec.alts
                        brecs, all_matched = _consume_batch_records_at(
                            pos, rec.ref, alts, batch_iters, batch_cur)

                        # Only biallelic sites present in every batch are aggregated;
                        # anything else keeps GLIMPSE's original value untouched.
                        if alts is not None and len(alts) == 1 and all_matched:
                            if _aggregate_record(rec, brecs, n_haps, n_hap_total,
                                                 af_tag, info_tag, decimals):
                                updated += 1
                        vcf_out.write(rec)
                        processed += 1
        finally:
            main_vcf.close()
            for bv in batch_vcfs:
                bv.close()

        worker_logger.info(f"Aggregate region done: {processed} variants written, "
                           f"{updated} AF/INFO aggregated.")
    except Exception:
        worker_logger.exception("Aggregate worker encountered a fatal error.")
        raise

    return processed


def _aggregate_sequential(records_from, batch_inputs, output_out, write_mode,
                          n_haps, af_tag, info_tag, decimals):
    """Single-process streaming aggregation (used when threads == 1)."""
    n_hap_total = sum(n_haps)
    main_vcf = pysam.VariantFile(records_from)
    batch_vcfs = [pysam.VariantFile(p) for p in batch_inputs]
    processed = 0
    updated = 0
    try:
        header = main_vcf.header
        _ensure_info_headers(header, af_tag, info_tag)
        batch_iters = [iter(bv) for bv in batch_vcfs]
        batch_cur = [next(it, None) for it in batch_iters]
        with pysam.VariantFile(output_out, write_mode, header=header) as vcf_out:
            for rec in main_vcf:
                brecs, all_matched = _consume_batch_records_at(
                    rec.pos, rec.ref, rec.alts, batch_iters, batch_cur)
                if rec.alts is not None and len(rec.alts) == 1 and all_matched:
                    if _aggregate_record(rec, brecs, n_haps, n_hap_total,
                                         af_tag, info_tag, decimals):
                        updated += 1
                vcf_out.write(rec)
                processed += 1
                if processed % 50000 == 0:
                    logger.info(f"Aggregated {processed} variants ...")
    finally:
        main_vcf.close()
        for bv in batch_vcfs:
            bv.close()
    logger.info(f"Sequential aggregate finished: {processed} variants written, "
                f"{updated} AF/INFO aggregated.")


def _aggregate_parallel(records_from, batch_inputs, output_out, write_mode, num_workers,
                        n_haps, af_tag, info_tag, decimals, tmp_dir):
    """
    Parallel `aggregate`: ensure indexes, split the records-from file into region
    chunks (reusing the same splitter as `recompute`), aggregate each region in a
    worker, then merge the temporaries in genomic order.
    """
    logger.info(f"Launching parallel aggregation with up to {num_workers} workers ...")

    # Region fetch needs an index on the records-from file AND on every batch file.
    if not _ensure_index(records_from):
        logger.error("Cannot proceed: records-from index is missing.")
        logging.shutdown()
        sys.exit(1)
    for bp in batch_inputs:
        if not _ensure_index(bp):
            logger.error(f"Cannot proceed: index missing for batch file '{bp}'.")
            logging.shutdown()
            sys.exit(1)

    chunks = _build_region_chunks(records_from, num_workers)
    actual_workers = len(chunks)
    if actual_workers == 0:
        logger.warning("No contigs with variant records found; writing empty output.")
        with pysam.VariantFile(records_from) as vf:
            header = vf.header
        with pysam.VariantFile(output_out, write_mode, header=header):
            pass
        return

    base_dir = tmp_dir if tmp_dir else os.path.dirname(os.path.abspath(output_out))
    temp_dir = tempfile.mkdtemp(prefix="glimpse_agg_", dir=base_dir)
    atexit.register(shutil.rmtree, temp_dir, ignore_errors=True)

    if output_out.endswith('.bcf'):
        temp_ext = '.bcf'
    elif output_out.endswith('.vcf.gz'):
        temp_ext = '.vcf.gz'
    else:
        temp_ext = '.vcf'

    worker_args = []
    for i, chunk in enumerate(chunks):
        temp_path = os.path.join(temp_dir, f"chunk_{i:04d}{temp_ext}")
        worker_args.append((records_from, batch_inputs, temp_path, chunk, write_mode,
                            n_haps, af_tag, info_tag, decimals))

    try:
        with multiprocessing.Pool(processes=actual_workers) as pool:
            results = pool.map(_aggregate_region_worker, worker_args)
        total = sum(r for r in results if isinstance(r, int))
        logger.info(f"All aggregate workers finished. Total variants processed: {total}")
    except Exception:
        logger.exception("Parallel aggregation aborted due to a worker failure.")
        logging.shutdown()
        sys.exit(1)

    temp_files = [os.path.join(temp_dir, f"chunk_{i:04d}{temp_ext}") for i in range(len(chunks))]
    temp_files = [f for f in temp_files if os.path.exists(f)]
    _merge_temp_files(temp_files, output_out, write_mode)

    shutil.rmtree(temp_dir, ignore_errors=True)
    logger.info("Parallel aggregation finished successfully.")


def aggregate_info(records_from: str, batch_inputs, output_out: str, threads: int,
                   af_tag: str, info_tag: str, decimals: int, tmp_dir: str) -> None:
    """
    Reconstruct cohort AF/INFO by combining the per-batch AF_b/INFO_b, then write
    them onto the records/GP of the full-cohort `records_from` file. Dispatches to
    the parallel path when threads > 1, otherwise streams sequentially.
    """
    write_mode = get_pysam_write_mode(output_out)
    if not batch_inputs:
        logger.error("aggregate requires at least one batch input file.")
        logging.shutdown()
        sys.exit(1)
    for p in [records_from] + list(batch_inputs):
        if not os.path.exists(p):
            logger.error(f"Input file not found: {p}")
            logging.shutdown()
            sys.exit(1)

    # Per-batch haplotype counts (2 * N_b), read once from each file's header.
    n_haps = []
    for bp in batch_inputs:
        with pysam.VariantFile(bp) as bf:
            n_haps.append(2 * len(bf.header.samples))
    n_hap_total = sum(n_haps)

    with pysam.VariantFile(records_from) as vf:
        n_records_samples = len(vf.header.samples)

    logger.info(f"Records-from : {records_from} ({n_records_samples} samples)")
    logger.info(f"Batches      : {len(batch_inputs)} files, per-batch n_hap = {n_haps}")
    logger.info(f"Cohort n_hap : {n_hap_total} (= 2 x {n_hap_total // 2} samples)")
    logger.info(f"Output       : {output_out} (mode: '{write_mode}')")

    # Sanity check: the merged records-from should hold exactly the union of batch
    # samples. A mismatch usually means the batch file list is wrong.
    if n_hap_total != 2 * n_records_samples:
        logger.warning(
            f"Sum of batch haplotypes ({n_hap_total}) != 2 x records-from samples "
            f"({2 * n_records_samples}). Aggregation uses the batch total as the "
            f"denominator; verify the batch list matches the merged cohort."
        )

    if threads > 1:
        _aggregate_parallel(records_from, batch_inputs, output_out, write_mode, threads,
                            n_haps, af_tag, info_tag, decimals, tmp_dir)
    else:
        _aggregate_sequential(records_from, batch_inputs, output_out, write_mode,
                              n_haps, af_tag, info_tag, decimals)


def _add_common_args(sp):
    """Attach the options shared by both sub-commands to a sub-parser."""
    sp.add_argument("-t", "--threads", type=int, default=1,
                    help="Number of parallel worker processes. When > 1 the work is split "
                         "by contig (or by position for single-chromosome / single-region "
                         "input) and merged back in genomic order.")
    sp.add_argument("--af-tag", type=str, default="AF",
                    help="INFO tag name for the ALT allele frequency.")
    sp.add_argument("--info-tag", type=str, default="INFO",
                    help="INFO tag name for the IMPUTE info score.")
    sp.add_argument("--decimals", type=int, default=3,
                    help="Decimal places for the rounded INFO score (GLIMPSE uses 3).")
    sp.add_argument("--tmp-dir", type=str, default=None,
                    help="Directory for intermediate chunk files (default: output directory).")


def main():
    parser = argparse.ArgumentParser(
        description="Recompute cohort-wide GLIMPSE2 INFO/AF and INFO/INFO for "
                    "sample-batched imputation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    sub = parser.add_subparsers(dest="command", required=True, help="Sub-command to run.")

    # --- recompute: single-file GP method (approximate for rare variants) ---
    p_rc = sub.add_parser(
        "recompute", formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Recompute AF/INFO from the full-cohort per-sample GP of ONE merged file. "
             "Simple, but rare-variant INFO is approximate (GLIMPSE stores GP floored "
             "to 3 decimals).")
    p_rc.add_argument("-i", "--input", required=True, type=str,
                      help="Input imputed file (.vcf/.vcf.gz/.bcf) with FORMAT/GP.")
    p_rc.add_argument("-o", "--output", required=True, type=str,
                      help="Output file; format inferred from extension.")
    _add_common_args(p_rc)

    # --- aggregate: per-batch method (accurate even for rare variants) ---
    p_ag = sub.add_parser(
        "aggregate", formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Combine per-batch AF/INFO (computed by GLIMPSE at full precision) into "
             "exact cohort AF/INFO. Accurate for rare variants. Run per chunk BEFORE "
             "the per-batch files are deleted.")
    p_ag.add_argument("--records-from", required=True, type=str,
                      help="Full-cohort file (e.g. the per-chunk bcftools-merge output) "
                           "whose records/GP/samples are kept; only its AF/INFO are "
                           "overwritten.")
    p_ag.add_argument("--batch-inputs", nargs="+", default=None, type=str,
                      help="Per-batch imputed files (identical site sets), each carrying "
                           "its own AF/INFO and samples.")
    p_ag.add_argument("--batch-list", default=None, type=str,
                      help="Alternative to --batch-inputs: a text file with one per-batch "
                           "file path per line ('#' comments allowed).")
    p_ag.add_argument("-o", "--output", required=True, type=str,
                      help="Output file; format inferred from extension.")
    _add_common_args(p_ag)

    args = parser.parse_args()

    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if args.decimals < 0:
        parser.error("--decimals must be >= 0")

    if args.command == "recompute":
        recompute_info(args.input, args.output, args.threads,
                       args.af_tag, args.info_tag, args.decimals, args.tmp_dir)
    elif args.command == "aggregate":
        batch_inputs = list(args.batch_inputs) if args.batch_inputs else []
        if args.batch_list:
            with open(args.batch_list) as fh:
                batch_inputs.extend(line.strip() for line in fh
                                    if line.strip() and not line.startswith("#"))
        if not batch_inputs:
            parser.error("aggregate requires --batch-inputs or --batch-list")
        aggregate_info(args.records_from, batch_inputs, args.output, args.threads,
                       args.af_tag, args.info_tag, args.decimals, args.tmp_dir)


if __name__ == "__main__":
    main()
