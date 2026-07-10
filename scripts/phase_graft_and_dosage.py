"""
https://github.com/ShujiaHuang/genotools/blob/main/scripts/phase_graft_and_dosage.py

Production-ready VCF/BCF Harmonization and Dosage Calculation Tool.
Synchronizes physical haplotypes (Phased GT) with observational variant data,
automatically detecting and supporting .vcf, .vcf.gz (BGZF), and .bcf formats.

Author: Shujia Huang
Date  : 2026-07-09


Examples:

python phase_graft_and_dosage.py -u HD_unphased.vcf.gz -p HD_phased.bcf -o HD_harmonized.bcf -t 4
python phase_graft_and_dosage.py -u HD_unphased.vcf.gz -p HD_phased.vcf.gz -o HD_harmonized.vcf.gz -t 4

"""

import sys
import os
import atexit
import shutil
import logging
import argparse
import tempfile
import multiprocessing
from typing import Tuple, Optional

import pysam


# Setup professional logging environment
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s [%(filename)s:%(lineno)d] - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger("VCF_Harmonizer")


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
            f"Unsupported file extension for output: '{filename}'. "
            f"Explicitly use '.vcf', '.vcf.gz', or '.bcf'."
        )


# Precomputed lookup table for 10^(-PL/10).
# PL values in GATK output are typically small integers (0–255); replacing math.pow
# with a dict lookup eliminates ~27 billion C-calls for large-scale datasets.
_PL_LUT: dict = {v: 10.0 ** (-v / 10.0) for v in range(256)}


def calculate_gp_ds(pl: Tuple[int, ...]) -> Tuple[Optional[Tuple[float, ...]], Optional[float]]:
    """
    Convert Phred-scaled likelihoods (PL) to Genotype Probabilities (GP) and Dosage (DS).
    Uses a precomputed lookup table for common PL values to avoid per-call math.pow.
    """
    if not pl:
        return None, None

    try:
        likelihoods = [_PL_LUT.get(val) if 0 <= val < 256 else 10.0 ** (-val / 10.0) for val in pl]
        sum_l = sum(likelihoods)
        
        if sum_l == 0:
            return None, None
            
        gps = tuple(l / sum_l for l in likelihoods)
        
        # Biallelic diploid case (Standard GATK joint-calling output layout)
        if len(pl) == 3:
            ds = gps[1] + 2.0 * gps[2]
            return tuple(round(p, 4) for p in gps), round(ds, 4)
            
        # Biallelic haploid case (e.g., Male sex chromosomes)
        elif len(pl) == 2:
            return tuple(round(p, 4) for p in gps), round(gps[1], 4)
            
        # True multi-allelic site: calculate GP but omit DS to prevent model bias
        elif len(pl) > 3:
            return tuple(round(p, 4) for p in gps), None

        # Unexpected PL length (e.g., malformed single-value PL); cannot compute valid GP/DS
        else:
            return None, None
            
    except Exception as e:
        logger.debug(f"Numerical error parsing PL {pl}: {e}")
        return None, None


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
                f"Cannot auto-index '{vcf_path}': only .vcf.gz and .bcf formats are supported. "
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
    Returns a list of lists, each sublist being one worker's assigned contigs.
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

    Uses index-backed random access (``fetch(contig)``) so this is
    O(num_active_contigs), not O(file_size).  Contigs listed in the VCF header
    but absent from the body (common in single-chromosome extracts) are excluded.
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


def _build_region_chunks(vcf_path, num_workers):
    """
    Build work chunks for parallel processing.

    Each chunk is a list of items that are EITHER:
      - a contig name (str)  →  process the entire contig
      - a (contig, start, end) tuple  →  process a position range within a contig

    Contig-level vs position-range splitting is decided dynamically based on
    which contigs *actually contain records*, not on the header contig list
    (which often lists the full reference genome even for single-chromosome extracts).
    """
    with pysam.VariantFile(vcf_path) as vf:
        header_contigs = list(vf.header.contigs.keys())
        contig_lengths = {}
        for c in header_contigs:
            try:
                length = vf.header.contigs[c].length
                contig_lengths[c] = length if length is not None else 0
            except (AttributeError, ValueError):
                contig_lengths[c] = 0

    # ── Determine which contigs actually contain variant records ──
    populated = set(_find_populated_contigs(vcf_path, header_contigs))
    contigs = [c for c in header_contigs if c in populated]

    if not contigs:
        logger.warning("No contigs with variant records found; nothing to parallelise.")
        return []

    # ── Contig-level chunking when there are enough populated contigs ──
    if len(contigs) >= num_workers:
        return _chunk_contigs(contigs, num_workers)

    # ── Fewer contigs than workers: split by position range ──
    # Some contigs may have unknown/zero length in the header but still contain
    # variants.  Those are kept as whole-contig chunks (1 worker each).
    # Active contigs (positive length) are split by position proportionally.
    inactive = [c for c in contigs if contig_lengths.get(c, 0) <= 0]
    active = [(c, contig_lengths[c]) for c in contigs if contig_lengths[c] > 0]

    if not active:
        return _chunk_contigs(contigs, num_workers)

    # Reserve one worker slot for each inactive contig
    slots_for_active = max(1, num_workers - len(inactive))
    total_length = sum(l for _, l in active)

    # Pre-calculate how many splits each active contig gets
    split_map: dict = {}
    remaining_slots = slots_for_active
    for idx, (contig, length) in enumerate(active):
        remaining_active = len(active) - idx
        if remaining_active == 1:
            n = max(1, remaining_slots)
        else:
            fraction = length / max(1, total_length)
            n = max(1, min(remaining_slots - (remaining_active - 1),
                          int(fraction * slots_for_active + 0.5)))
        split_map[contig] = n
        remaining_slots -= n

    # Build chunks in VCF-header contig order (preserves genomic sort order)
    chunks = []
    for contig in contigs:
        if contig in split_map:
            length = contig_lengths[contig]
            n_splits = min(split_map[contig], length)  # never split finer than 1 bp
            if n_splits <= 1:
                chunks.append([contig])
            else:
                split_len = max(1, length // n_splits)
                for i in range(n_splits):
                    start = i * split_len + 1
                    end = (i + 1) * split_len if i < n_splits - 1 else length
                    chunks.append([(contig, start, end)])
        else:
            # Inactive contig (zero/None length) — whole-contig chunk
            chunks.append([contig])

    # Guard: merge excess chunks if we overshot
    while len(chunks) > num_workers:
        chunks[-2].extend(chunks[-1])
        chunks.pop()

    return [c for c in chunks if c]


def _process_region_worker(args):
    """
    Multiprocessing worker: processes a subset of contigs and writes to a
    temporary output file.  Executed in an independent child process.
    """
    unphased_path, phased_path, temp_output_path, contigs, write_mode = args

    # Configure per-worker logging
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s [Worker-%(process)d] - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    worker_logger = logging.getLogger("VCF_Harmonizer.Worker")

    worker_logger.info(
        f"Assigned {len(contigs)} region(s): "
        f"{str(contigs[:5])[:-1]}{', ...' if len(contigs) > 5 else ']'}"
    )

    processed = 0
    try:
        with pysam.VariantFile(unphased_path) as unphased_vcf, \
             pysam.VariantFile(phased_path) as phased_vcf:

            header = unphased_vcf.header
            if "GT" not in header.formats:
                header.formats.add("GT", "1", "String", "Genotype")
            if "GP" not in header.formats:
                header.formats.add("GP", "G", "Float", "Genotype Probabilities derived strictly from PL")
            if "DS" not in header.formats:
                header.formats.add("DS", "1", "Float", "Alternate Allele Dosage calculated from GP")

            unphased_samples = list(header.samples)

            with pysam.VariantFile(temp_output_path, write_mode, header=header) as out_vcf:
                for item in contigs:
                    # Each item is either a contig name (str) or a position-range tuple
                    if isinstance(item, tuple):
                        contig, start, end = item
                        label = f"{contig}:{start}-{end}"
                    else:
                        contig, start, end = item, None, None
                        label = contig

                    try:
                        if start is not None and end is not None:
                            # pysam.fetch() uses 0-based half-open coordinates.
                            # _build_region_chunks emits 1-based [start, end];
                            # subtract 1 from start to avoid silently dropping
                            # the first variant of every position-range chunk.
                            u_iter = unphased_vcf.fetch(contig, start - 1, end)
                            p_iter = phased_vcf.fetch(contig, start - 1, end)
                        else:
                            u_iter = unphased_vcf.fetch(contig)
                            p_iter = phased_vcf.fetch(contig)
                    except ValueError:
                        worker_logger.debug(f"Region '{label}' not found in index; skipping.")
                        continue

                    u_rec = next(u_iter, None)
                    p_rec = next(p_iter, None)

                    while u_rec is not None and p_rec is not None:
                        if u_rec.contig != p_rec.contig or u_rec.pos != p_rec.pos:
                            worker_logger.warning(
                                f"Sync mismatch at {u_rec.contig}:{u_rec.pos} vs "
                                f"{p_rec.contig}:{p_rec.pos}; aborting region '{label}'."
                            )
                            break

                        for sample in unphased_samples:
                            u_fmt = u_rec.samples[sample]
                            p_fmt = p_rec.samples[sample]

                            gt = p_fmt.get("GT")
                            if gt is not None:
                                u_fmt["GT"] = gt

                            u_fmt.phased = getattr(p_fmt, "phased", False)

                            pl = u_fmt.get("PL")
                            if pl is not None:
                                gp, ds = calculate_gp_ds(pl)
                                if gp is not None:
                                    u_fmt["GP"] = gp
                                if ds is not None:
                                    u_fmt["DS"] = ds

                        out_vcf.write(u_rec)
                        processed += 1

                        u_rec = next(u_iter, None)
                        p_rec = next(p_iter, None)

                    if u_rec is not None or p_rec is not None:
                        worker_logger.warning(
                            f"Region '{label}': record count mismatch between input files."
                        )

        # Build concise region labels for the completion log
        region_labels = [f"{c[0]}:{c[1]}-{c[2]}" if isinstance(c, tuple) else c for c in contigs[:5]]
        suffix = f", ... ({len(contigs)} total)" if len(contigs) > 5 else ""
        worker_logger.info(
            f"Completed. Processed {processed} variants across "
            f"{', '.join(region_labels)}{suffix}."
        )
    except Exception:
        worker_logger.exception("Worker encountered a fatal error.")
        raise

    return processed


def _merge_temp_files(temp_files, output_path, write_mode):
    """
    Merge multiple temporary VCF/BCF files into a single final output file.
    Writes the header once and streams records from each temp file in order.
    """
    logger.info(f"Merging {len(temp_files)} temporary chunk(s) into '{output_path}' ...")

    if not temp_files:
        logger.error("No temporary files to merge.")
        return

    # Read header from the first temp file
    with pysam.VariantFile(temp_files[0]) as first:
        header = first.header

    merged_count = 0
    with pysam.VariantFile(output_path, write_mode, header=header) as out_vcf:
        for temp_file in temp_files:
            with pysam.VariantFile(temp_file) as tf:
                for rec in tf:
                    out_vcf.write(rec)
                    merged_count += 1

    logger.info(f"Merge complete: {merged_count} total variants written to output.")


def _harmonize_parallel(unphased_in, phased_in, output_out, write_mode, num_workers):
    """
    Orchestrate parallel harmonization:
      1. Ensure indices exist
      2. Build region chunks (auto-splits large contigs by position when
         there are fewer contigs than workers, e.g. single-chromosome VCF)
      3. Process in parallel via multiprocessing
      4. Merge temporary output files
    """
    logger.info(f"Launching parallel harmonization with {num_workers} workers ...")

    # ── 1. Ensure region-query indices exist ──
    if not _ensure_index(unphased_in) or not _ensure_index(phased_in):
        logger.error("Cannot proceed with parallel processing: required indices are missing.")
        logging.shutdown()
        sys.exit(1)

    # ── 2. Gather contigs and distribute work (auto-splits by position when needed) ──
    chunks = _build_region_chunks(unphased_in, num_workers)
    actual_workers = len(chunks)

    if actual_workers == 0:
        logger.warning("No contigs with variant records found; writing empty output.")
        with pysam.VariantFile(unphased_in) as vf:
            header = vf.header
        with pysam.VariantFile(output_out, write_mode, header=header) as _:
            pass
        return

    # Count statistics for logging
    n_contig_chunks = sum(1 for c in chunks if all(isinstance(x, str) for x in c))
    n_region_chunks = sum(1 for c in chunks if any(isinstance(x, tuple) for x in c))
    logger.info(
        f"Work distributed across {actual_workers} worker(s): "
        f"{n_contig_chunks} contig-level chunk(s), {n_region_chunks} position-range chunk(s)."
    )

    # ── 3. Prepare temp directory for intermediate outputs ──
    #     Place under the output file's directory instead of system /tmp to
    #     avoid tmpfs size limits on HPC nodes.
    output_dir = os.path.dirname(os.path.abspath(output_out))
    temp_dir = tempfile.mkdtemp(prefix="vcf_harmonize_", dir=output_dir)
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
        worker_args.append((unphased_in, phased_in, temp_path, chunk, write_mode))

    # ── 4. Execute workers ──
    try:
        with multiprocessing.Pool(processes=actual_workers) as pool:
            results = pool.map(_process_region_worker, worker_args)
            
        total = sum(r for r in results if isinstance(r, int))
        if total == 0:
            logger.error("All workers returned zero variants; aborting.")
            logging.shutdown()
            sys.exit(1)
        logger.info(f"All workers finished successfully. Total variants processed: {total}")
    except Exception:
        logger.exception("Parallel processing aborted due to a worker failure.")
        logging.shutdown()
        sys.exit(1)

    # ── 5. Merge intermediate files into final output ──
    temp_files = [
        os.path.join(temp_dir, f"chunk_{i:04d}{temp_ext}")
        for i in range(len(chunks))
    ]
    temp_files = [f for f in temp_files if os.path.exists(f)]
    _merge_temp_files(temp_files, output_out, write_mode)

    # ── 6. Cleanup ──
    shutil.rmtree(temp_dir, ignore_errors=True)
    logger.info("Parallel harmonization finished successfully.")


def harmonize_cohorts(unphased_in: str, phased_in: str, output_out: str, threads: int = 1) -> None:
    """
    High-throughput synchronous streaming merge between the unphased baseline
    and phased reference data. Injects exact phase states, computes accurate GP/DS values,
    and dynamically stream-compresses to the target output format.

    When *threads > 1*, automatically splits input by chromosome for parallel
    multi-process execution and merges results upon completion.
    """
    write_mode = get_pysam_write_mode(output_out)
    logger.info(f"Input Unphased Source: {unphased_in}")
    logger.info(f"Input Phased Reference: {phased_in}")
    logger.info(f"Output Destination: {output_out} (Mode: '{write_mode}')")

    # ── Pre-flight: validate sample consistency ──
    with pysam.VariantFile(unphased_in) as unphased_check, \
         pysam.VariantFile(phased_in) as phased_check:
        unphased_samples = list(unphased_check.header.samples)
        phased_samples = list(phased_check.header.samples)
        if set(unphased_samples) != set(phased_samples):
            logger.error(
                f"Sample mismatch error! Unphased samples: {unphased_samples[:5]}..., "
                f"Phased samples: {phased_samples[:5]}..."
            )
            logging.shutdown()
            sys.exit(1)
        logger.info(f"Cohort validation complete ({len(unphased_samples)} samples matched).")

    # ── Parallel path ──
    if threads > 1:
        _harmonize_parallel(unphased_in, phased_in, output_out, write_mode, threads)
        return

    # ── Sequential path (original logic for threads=1) ──
    with pysam.VariantFile(unphased_in, "r") as unphased_vcf, \
         pysam.VariantFile(phased_in, "r") as phased_vcf:

        header = unphased_vcf.header
        if "GT" not in header.formats:
            header.formats.add("GT", "1", "String", "Genotype")
        if "GP" not in header.formats:
            header.formats.add("GP", "G", "Float", "Genotype Probabilities derived strictly from PL")
        if "DS" not in header.formats:
            header.formats.add("DS", "1", "Float", "Alternate Allele Dosage calculated from GP")

        logger.info("Executing sequential stream ...")

        processed_count = 0
        with pysam.VariantFile(output_out, write_mode, header=header) as out_vcf:
            u_rec = next(unphased_vcf, None)
            p_rec = next(phased_vcf, None)

            while u_rec is not None and p_rec is not None:
                if u_rec.contig != p_rec.contig or u_rec.pos != p_rec.pos:
                    logger.critical(
                        f"CRITICAL DESYNCHRONIZATION DETECTED!\n"
                        f"Unphased iterator at -> {u_rec.contig}:{u_rec.pos}\n"
                        f"Phased iterator at   -> {p_rec.contig}:{p_rec.pos}\n"
                        f"Please ensure both input files contain identical variants and sorting orders."
                    )
                    logging.shutdown()
                    sys.exit(1)

                for sample in unphased_samples:
                    u_fmt = u_rec.samples[sample]
                    p_fmt = p_rec.samples[sample]

                    gt = p_fmt.get("GT")
                    if gt is not None:
                        u_fmt["GT"] = gt

                    phased_flag = getattr(p_fmt, "phased", False)
                    u_fmt.phased = phased_flag

                    pl = u_fmt.get("PL")
                    if pl is not None:
                        gp, ds = calculate_gp_ds(pl)
                        if gp is not None:
                            u_fmt["GP"] = gp
                        if ds is not None:
                            u_fmt["DS"] = ds

                out_vcf.write(u_rec)
                processed_count += 1

                if processed_count % 50000 == 0:
                    logger.info(f"Successfully processed {processed_count} genomic variants.")

                u_rec = next(unphased_vcf, None)
                p_rec = next(phased_vcf, None)

            if u_rec is not None or p_rec is not None:
                logger.error("Input files contain a different number of records; processing stopped early.")
                logging.shutdown()
                sys.exit(1)

    logger.info(f"Execution finished. Total processed variants: {processed_count}.")


def main():
    parser = argparse.ArgumentParser(
        description="Enterprise-grade VCF/BCF format harmonizer for advanced genomic analysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-u", "--unphased", required=True, type=str,
                        help="Input unphased variant file (.vcf, .vcf.gz, or .bcf) containing raw PL data.")
    parser.add_argument("-p", "--phased", required=True, type=str,
                        help="Input phased reference file (.vcf, .vcf.gz, or .bcf) with high-confidence phase states.")
    parser.add_argument("-o", "--output", required=True, type=str,
                        help="Output harmonized variant file. Format automatically inferred from extension (.vcf, .vcf.gz, .bcf).")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of parallel worker processes. When > 1, input is automatically "
                             "split by chromosome for multi-process execution and results are merged.")
    
    args = parser.parse_args()
    harmonize_cohorts(args.unphased, args.phased, args.output, args.threads)


if __name__ == "__main__":
    main()



