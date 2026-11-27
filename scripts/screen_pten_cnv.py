#!/usr/bin/env python3
"""
Fast multi-process CNV (microdeletion/microduplication) screener for the chr10
PTEN (10q23.31) locus across large-scale low-pass NIPT BAM/CRAM cohorts.

Uses ONLY pysam for BAM/CRAM and FASTA operations, combined with:
  1) Library-size normalization (autosomal mapped reads from the index; CRAI
     indexes do not record per-contig counters, so CRAM input falls back to the
     number of usable reads inside the target region)
  2) Per-sample GC-content stratification correction (when --ref-fasta is given)
  3) Cohort-wide bin median & robust MAD Z-score normalization

Author: Shujia Huang and Gemini & Qoder
Date: 2026-11-27
"""
import argparse
import csv
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pysam

# GRCh38 coordinates for PTEN and adjacent 10q23.31 context
PTEN_CORE_START = 87863625
PTEN_CORE_END = 87971930
DEFAULT_SCAN_START = 80000000
DEFAULT_SCAN_END = 96000000
AUTOSOMES_CLEAN = {str(i) for i in range(1, 23)}

# Minimum usable reads per bin when the library size is estimated from the scan
# window (CRAM input, whose .crai index carries no per-contig read counters).
MIN_REGION_READS_PER_BIN = 5


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Screen chr10 PTEN CNVs from NIPT BAM/CRAM files using pysam.")
    p.add_argument("--bam-list", required=True, type=Path, help="File with 'sample_id<TAB>bam_path' or one bam_path per line.")
    p.add_argument("-o", "--out-prefix", required=True, type=Path, help="Output prefix for TSV reports.")
    p.add_argument("--ref-fasta", type=Path, default=None, help="GRCh38 reference FASTA (required for CRAM and GC correction).")
    p.add_argument("--chrom", type=str, default="chr10", help="Target chromosome name (default: chr10).")
    p.add_argument("--scan-start", type=int, default=DEFAULT_SCAN_START, help="1-based flanking scan start (default: 80000000).")
    p.add_argument("--scan-end", type=int, default=DEFAULT_SCAN_END, help="1-based flanking scan end (default: 96000000).")
    p.add_argument("--bin-size", type=int, default=50000, help="Sliding window size in bp (default: 50000).")
    p.add_argument("--step-size", type=int, default=25000, help="Sliding window step size in bp (default: 25000).")
    p.add_argument("--min-mapq", type=int, default=30, help="Minimum read mapping quality (default: 30).")
    p.add_argument("--max-n-frac", type=float, default=0.20, help="Max allowed 'N' base fraction per bin (default: 0.20).")
    p.add_argument("--z-thresh", type=float, default=4.0, help="Absolute Z-score threshold for candidate CNV calling (default: 4.0).")
    p.add_argument("-t", "--threads", type=int, default=16, help="Number of parallel worker processes (default: 16).")
    p.add_argument("--min-library-reads", type=int, default=10000, help="Minimum library size of a sample: autosomal mapped reads from the index, or region reads for CRAM (default: 10000).")
    p.add_argument("--min-valid-samples", type=int, default=10, help="Minimum number of QC-passing samples required for cohort normalization (default: 10).")
    p.add_argument("--save-matrix", action="store_true", help="Save normalized copy-ratio matrix (.tsv.gz) for downstream segmentation.")
    return p.parse_args()


def clean_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.startswith("chr") else chrom


def resolve_chrom_name(chrom: str, available_contigs: Tuple[str, ...]) -> str:
    if chrom in available_contigs:
        return chrom
    
    alt = clean_chrom(chrom) if chrom.startswith("chr") else f"chr{chrom}"
    if alt in available_contigs:
        return alt
    
    raise ValueError(f"Chromosome '{chrom}' not found in header contigs.")


def load_sample_bams(bam_list_path: Path) -> List[Tuple[str, str]]:
    samples: List[Tuple[str, str]] = []
    seen: set = set()
    with bam_list_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) >= 2:
                sid, bpath = parts[0].strip(), parts[1].strip()
            else:
                bpath = parts[0].strip()
                sid = Path(bpath).name.split(".")[0]
            
            if sid in seen:
                sid = f"{sid}_{len(samples)}"
            
            seen.add(sid)
            samples.append((sid, bpath))
            
    if not samples:
        sys.exit(f"[ERROR] No valid BAM/CRAM entries found in {bam_list_path}")
        
    return samples


def build_bins(
    chrom: str,
    scan_start: int,
    scan_end: int,
    bin_size: int,
    step_size: int,
    ref_fasta: Optional[Path],
    max_n_frac: float,
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], np.ndarray]:
    """
    Generate valid sliding bins across [scan_start, scan_end], compute GC content
    via pysam.FastaFile if provided, and identify indices overlapping PTEN core.
    """
    starts: List[int] = []
    ends: List[int] = []
    gc_vals: List[float] = []

    fasta = pysam.FastaFile(str(ref_fasta)) if ref_fasta else None
    fa_chrom = chrom
    if fasta:
        try:
            fa_chrom = resolve_chrom_name(chrom, tuple(fasta.references))
        except ValueError:
            fasta.close()
            sys.exit(f"[ERROR] Chromosome '{chrom}' not found in --ref-fasta '{ref_fasta}'. "
                     "Use a reference that contains the target contig, or drop --ref-fasta.")

    for pos in range(scan_start, scan_end - bin_size + 1, step_size):
        b_start = pos
        b_end = pos + bin_size - 1
        if fasta:
            seq = fasta.fetch(fa_chrom, b_start - 1, b_end).upper()
            n_count = seq.count("N")
            if len(seq) == 0 or (n_count / len(seq)) > max_n_frac:
                continue
            valid_bases = len(seq) - n_count
            gc = (seq.count("G") + seq.count("C")) / valid_bases if valid_bases > 0 else 0.5
            gc_vals.append(gc)
        starts.append(b_start)
        ends.append(b_end)

    if fasta:
        fasta.close()

    bin_starts = np.array(starts, dtype=np.int64)
    bin_ends = np.array(ends, dtype=np.int64)
    gc_array = np.array(gc_vals, dtype=np.float64) if ref_fasta else None

    # Identify bins overlapping the PTEN coding/genomic locus
    pten_mask = (bin_ends >= PTEN_CORE_START) & (bin_starts <= PTEN_CORE_END)
    pten_indices = np.where(pten_mask)[0]
    if len(pten_indices) == 0:
        sys.exit("[ERROR] No valid bins overlap the PTEN core locus. Check --scan-start/--scan-end.")

    return bin_starts, bin_ends, gc_array, pten_indices


def index_library_depth(bam: pysam.AlignmentFile) -> Optional[int]:
    """
    Sum of mapped reads on autosomes as recorded by the BAM/CRAM index.

    Returns None when the index carries no per-contig counters: CRAI indexes do not
    store them, so pysam reports 0 mapped reads for every contig of a CRAM file
    (see pysam.AlignmentFile.get_index_statistics). Callers must fall back to a
    region-based library size in that case.
    """
    try:
        stats = bam.get_index_statistics()
    except (ValueError, NotImplementedError):
        return None

    total = sum(stat.mapped for stat in stats if clean_chrom(stat.contig) in AUTOSOMES_CLEAN)
    return total if total > 0 else None


def process_single_bam(
    task: Tuple[int, str, str],
    chrom: str,
    scan_start: int,
    scan_end: int,
    bin_starts: np.ndarray,
    bin_ends: np.ndarray,
    gc_stratum_ids: Optional[np.ndarray],
    ref_fasta_str: Optional[str],
    min_mapq: int,
    min_library_reads: int,
) -> Tuple[int, str, Optional[np.ndarray], int, str, str]:
    """
    Worker function: count high-MAPQ reads in target region bins via pysam, apply
    library-size and per-sample GC normalization.

    Returns a tuple of
    (task_index, sample_id, normalized_counts, library_reads, library_source, status),
    where library_source is 'index' (autosomal mapped reads) or 'region' (usable
    reads in the scan window, used when the index has no counters, e.g. CRAM).
    """
    idx, sid, bam_path = task
    library_reads = 0
    library_source = "index"

    try:
        bam = pysam.AlignmentFile(bam_path, reference_filename=ref_fasta_str)
    except Exception as exc:
        return (idx, sid, None, 0, library_source, f"OpenError: {exc}")

    try:
        if not bam.has_index():
            return (idx, sid, None, 0, library_source,
                    "NoIndex: missing .bai/.csi/.crai index (run samtools index)")

        # 1. Library size from the index (fast, exact for BAM); None for CRAM
        library_reads = index_library_depth(bam)
        if library_reads is not None and library_reads < min_library_reads:
            return (idx, sid, None, library_reads, library_source,
                    f"LowLibrary: {library_reads} autosomal mapped reads < {min_library_reads}")

        # 2. Fetch reads across the target chr10 region
        bam_chrom = resolve_chrom_name(chrom, bam.references)
        read_positions: List[int] = []
        for aln in bam.fetch(bam_chrom, max(0, scan_start - 1), scan_end):
            if (
                aln.is_unmapped
                or aln.is_secondary
                or aln.is_supplementary
                or aln.is_duplicate
                or aln.is_qcfail
                or aln.mapping_quality < min_mapq
            ): continue
            read_positions.append(aln.reference_start + 1)  # 1-based start

        if not read_positions:
            # No usable read in the whole scan window: a zero bin profile would look
            # like a genome-wide deletion, so reject the sample instead of calling it.
            return (idx, sid, None, library_reads, library_source,
                    "NoSignal: no usable read in the scan window")

        if library_reads is None:
            # CRAI indexes record no per-contig counters, so fall back to the usable
            # reads of the scan window. Both steps below (cohort bin median and
            # per-sample median re-centering) remove any per-sample scale factor,
            # so CNV detection keeps the same relative sensitivity.
            library_source = "region"
            library_reads = len(read_positions)
            # A window-based library size cannot use --min-library-reads directly:
            # the read count scales with the window size and the sequencing depth.
            # Require a minimum coverage per bin instead, never stricter than the flag.
            min_region_reads = min(min_library_reads, max(MIN_REGION_READS_PER_BIN * len(bin_starts), 200))
            if library_reads < min_region_reads:
                return (idx, sid, None, library_reads, library_source,
                        f"LowLibrary: {library_reads} usable region reads < {min_region_reads}")

        # 3. Assign reads to sliding bins via binary search
        r_arr = np.sort(np.array(read_positions, dtype=np.int64))
        # For each bin [b_start, b_end], count reads with start in [b_start, b_end]
        left_idx = np.searchsorted(r_arr, bin_starts, side="left")
        right_idx = np.searchsorted(r_arr, bin_ends, side="right")
        counts = (right_idx - left_idx).astype(np.float64)

        # Normalize by library size (reads per million library reads)
        norm_counts = (counts / library_reads) * 1e6

        # 4. Per-sample GC-stratum median correction (if GC strata provided)
        if gc_stratum_ids is not None:
            corrected = norm_counts.copy()
            positive = norm_counts[norm_counts > 0]
            global_med = float(np.median(positive)) if positive.size else 1.0
            for stratum in np.unique(gc_stratum_ids):
                smask = gc_stratum_ids == stratum
                s_vals = norm_counts[smask]
                # Zero-count bins are not informative about the stratum level, so
                # the stratum median is taken over positive bins only.
                s_pos = s_vals[s_vals > 0]
                s_med = float(np.median(s_pos)) if s_pos.size else 0.0
                if s_med > 0:
                    corrected[smask] = (s_vals / s_med) * global_med
                    
            norm_counts = corrected

        return (idx, sid, norm_counts, library_reads, library_source, "OK")

    except Exception as exc:
        msg = f"ReadError: {type(exc).__name__}: {exc}"
        if bam.is_cram and not ref_fasta_str:
            msg += " (CRAM needs --ref-fasta)"
        
        return (idx, sid, None, library_reads, library_source, msg)
    
    finally:
        bam.close()



def estimate_cnv_span(
    z_row: np.ndarray,
    ratio_row: np.ndarray,
    bin_starts: np.ndarray,
    bin_ends: np.ndarray,
    pten_indices: np.ndarray,
    z_thresh: float,
) -> Tuple[str, int, int, int]:
    """Expand contiguous abnormal bins outward from PTEN core to estimate CNV boundaries."""
    pten_z = z_row[pten_indices]
    finite = np.isfinite(pten_z)
    if finite.any():
        # Start from the strongest PTEN bin, ignoring bins without usable signal
        center_idx = int(pten_indices[finite][int(np.argmax(np.abs(pten_z[finite])))])
    else:
        center_idx = int(pten_indices[0])
    direction = -1 if z_row[center_idx] < 0 else 1
    call_type = "DEL" if direction == -1 else "DUP"

    # Expand left and right while bins remain consistently shifted (|Z| >= 2.5)
    soft_thresh = min(2.5, z_thresh * 0.6)
    left = center_idx
    while left > 0 and (z_row[left - 1] * direction) >= soft_thresh:
        left -= 1
    right = center_idx
    while right < len(z_row) - 1 and (z_row[right + 1] * direction) >= soft_thresh:
        right += 1

    est_start = int(bin_starts[left])
    est_end = int(bin_ends[right])
    est_size = est_end - est_start + 1
    return call_type, est_start, est_end, est_size


def fmt_num(value: float, digits: int) -> str:
    """Format a float for TSV output, writing NA for non-finite values."""
    return f"{value:.{digits}f}" if np.isfinite(value) else "NA"


ALL_COLUMNS = [
    "sample_id", "status", "library_reads", "library_norm_source", "flanking_chr10_mad",
    "pten_mean_copy_ratio", "pten_mean_zscore", "pten_peak_zscore",
    "pten_abnormal_bins", "total_pten_bins", "cnv_call",
    "est_cnv_start", "est_cnv_end", "est_cnv_size_bp",
]


def qc_failed_row(
    sid: str, status: str, library_reads: int, library_source: str, n_pten_bins: int, cnv_call: str = "FAILED_QC"
) -> Dict[str, object]:
    """Per-sample report row for a sample that was not evaluated against the cohort."""
    return {
        "sample_id": sid, "status": status,
        "library_reads": library_reads, "library_norm_source": library_source,
        "flanking_chr10_mad": "NA", "pten_mean_copy_ratio": "NA",
        "pten_mean_zscore": "NA", "pten_peak_zscore": "NA",
        "pten_abnormal_bins": "NA", "total_pten_bins": n_pten_bins,
        "cnv_call": cnv_call, "est_cnv_start": "NA",
        "est_cnv_end": "NA", "est_cnv_size_bp": "NA",
    }


def write_report(path: Path, rows: List[Dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=ALL_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def print_status_summary(statuses: List[str]) -> None:
    """Print how many samples failed for each distinct reason (most frequent first)."""
    failures = Counter(s for s in statuses if s != "OK")
    if not failures:
        print("    all samples passed QC: the cohort itself is too small, add samples or lower "
              "--min-valid-samples", file=sys.stderr, flush=True)
        return

    # Statuses carry per-sample detail (read counts, error text), so group them by code
    groups: Dict[str, List[object]] = {}
    for status, count in failures.items():
        code = status.split(":")[0]
        total, example = groups.get(code, [0, status])
        groups[code] = [total + count, example]
    for code, (count, example) in sorted(groups.items(), key=lambda item: -item[1][0]):
        print(f"    {count:>8} sample(s): {example}", file=sys.stderr, flush=True)


def main() -> None:
    args = parse_args()
    if args.ref_fasta and not args.ref_fasta.exists():
        sys.exit(f"[ERROR] Reference FASTA not found: {args.ref_fasta}")

    sample_bams = load_sample_bams(args.bam_list)
    total_input = len(sample_bams)
    print(f"[INFO] Loaded {total_input} BAM/CRAM entries from {args.bam_list}", file=sys.stderr)

    bin_starts, bin_ends, gc_array, pten_indices = build_bins(
        args.chrom, args.scan_start, args.scan_end, args.bin_size, args.step_size, args.ref_fasta, args.max_n_frac
    )
    n_bins = len(bin_starts)
    print(f"[INFO] Constructed {n_bins} valid bins ({len(pten_indices)} overlapping PTEN core).", file=sys.stderr)

    # Discretize GC content into 2% strata (0..50) for robust per-sample median correction
    gc_stratum_ids = np.floor(gc_array * 50).astype(np.int32) if gc_array is not None else None
    ref_fasta_str = str(args.ref_fasta) if args.ref_fasta else None

    raw_matrix = np.full((total_input, n_bins), np.nan, dtype=np.float64)
    lib_reads = np.zeros(total_input, dtype=np.int64)
    lib_source = ["NA"] * total_input
    status_list = ["Pending"] * total_input

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    all_out = Path(f"{args.out_prefix}.all_samples.tsv")
    cand_out = Path(f"{args.out_prefix}.candidates.tsv")

    tasks = [(i, sid, bpath) for i, (sid, bpath) in enumerate(sample_bams)]
    completed = 0
    with ProcessPoolExecutor(max_workers=args.threads) as pool:
        futures = {
            pool.submit(
                process_single_bam,
                t,
                args.chrom,
                args.scan_start,
                args.scan_end,
                bin_starts,
                bin_ends,
                gc_stratum_ids,
                ref_fasta_str,
                args.min_mapq,
                args.min_library_reads,
            ): t[0]
            for t in tasks
        }
        for fut in as_completed(futures):
            task_idx = futures[fut]
            try:
                idx, sid, norm_counts, n_reads, depth_source, status = fut.result()
            except Exception as exc:
                # A crashed worker (segfault, OOM, broken pool) must not abort the cohort
                status_list[task_idx] = f"WorkerError: {type(exc).__name__}: {exc}"
                completed += 1
                continue
            status_list[idx] = status
            lib_reads[idx] = n_reads
            lib_source[idx] = depth_source
            if norm_counts is not None:
                raw_matrix[idx, :] = norm_counts
            completed += 1
            if completed % 5000 == 0 or completed == total_input:
                print(f"[PROGRESS] Processed {completed}/{total_input} samples...", file=sys.stderr)

    valid_mask = np.array([s == "OK" for s in status_list])
    valid_indices = np.where(valid_mask)[0]
    n_valid = len(valid_indices)
    if n_valid < args.min_valid_samples:
        # Write the QC report before exiting, otherwise the reason for every failed
        # sample is lost and the run is impossible to diagnose.
        write_report(all_out, [
            qc_failed_row(
                sid, status_list[i], int(lib_reads[i]), lib_source[i], len(pten_indices),
                "FAILED_QC" if status_list[i] != "OK" else "NOT_EVALUATED",
            )
            for i, (sid, _) in enumerate(sample_bams)
        ])
        print(f"[ERROR] Too few valid samples ({n_valid}/{total_input} < {args.min_valid_samples}) "
              f"for cohort normalization.", file=sys.stderr)
        print(f"[INFO] Per-sample status written to {all_out}", file=sys.stderr)
        print_status_summary(status_list)
        sys.exit(1)

    mat = raw_matrix[valid_indices, :]

    # Cohort-wide bin median normalization (eliminates residual mappability & GC bias)
    bin_medians = np.median(mat, axis=0)
    bin_medians[bin_medians <= 0] = np.nan
    copy_ratio = mat / bin_medians

    # Re-center each sample by its own median copy ratio across the flanking chr10 region
    sample_medians = np.nanmedian(copy_ratio, axis=1, keepdims=True)
    sample_medians = np.where(np.isfinite(sample_medians) & (sample_medians > 0), sample_medians, 1.0)
    copy_ratio = copy_ratio / sample_medians

    # Robust bin-wise MAD standard deviation & Z-score matrix.
    # Bins with a near-zero MAD (all-zero or perfectly flat bins) must not produce
    # an infinite Z, nor be silently dropped: use a floor derived from the cohort.
    bin_mad = np.nanmedian(np.abs(copy_ratio - 1.0), axis=0)
    robust_sd = bin_mad * 1.4826
    finite_sd = robust_sd[np.isfinite(robust_sd) & (robust_sd > 0)]
    sd_floor = max(float(np.median(finite_sd)) * 0.1, 1e-6) if finite_sd.size else 1e-6
    robust_sd = np.where(np.isfinite(robust_sd) & (robust_sd > sd_floor), robust_sd, sd_floor)
    z_matrix = (copy_ratio - 1.0) / robust_sd

    # Evaluate PTEN core metrics per valid sample
    candidate_rows: List[Dict[str, object]] = []
    report_rows: List[Dict[str, object]] = []

    v_ptr = 0
    for orig_idx, (sid, _) in enumerate(sample_bams):
        if not valid_mask[orig_idx]:
            report_rows.append(qc_failed_row(
                sid, status_list[orig_idx], int(lib_reads[orig_idx]), lib_source[orig_idx], len(pten_indices)
            ))
            continue

        cr_row = copy_ratio[v_ptr, :]
        z_row = z_matrix[v_ptr, :]
        v_ptr += 1

        flank_mad = float(np.nanmedian(np.abs(cr_row - 1.0)))
        pten_cr = float(np.nanmean(cr_row[pten_indices]))
        pten_z_vals = z_row[pten_indices]
        finite_z = np.isfinite(pten_z_vals)
        if finite_z.any():
            abn_bins = int(np.sum(np.abs(pten_z_vals[finite_z]) >= args.z_thresh))
            pten_mean_z = float(np.mean(pten_z_vals[finite_z]))
            pten_peak_z = float(pten_z_vals[finite_z][int(np.argmax(np.abs(pten_z_vals[finite_z])))])
        else:
            abn_bins = 0
            pten_mean_z = pten_peak_z = float("nan")

        # Call candidate CNV if >=2 PTEN bins exceed Z-threshold or mean |Z| >= z_thresh
        is_del = (pten_mean_z <= -args.z_thresh or (abn_bins >= 2 and pten_peak_z <= -args.z_thresh)) and pten_cr <= 0.80
        is_dup = (pten_mean_z >= args.z_thresh or (abn_bins >= 2 and pten_peak_z >= args.z_thresh)) and pten_cr >= 1.20

        if is_del or is_dup:
            ctype, c_start, c_end, c_size = estimate_cnv_span(
                z_row, cr_row, bin_starts, bin_ends, pten_indices, args.z_thresh
            )
        elif not np.isfinite(pten_cr):
            # No usable signal in the PTEN bins: report it instead of calling NORMAL
            ctype, c_start, c_end, c_size = "NO_SIGNAL", 0, 0, 0
        else:
            ctype, c_start, c_end, c_size = "NORMAL", 0, 0, 0

        row_dict: Dict[str, object] = {
            "sample_id": sid,
            "status": "OK",
            "library_reads": int(lib_reads[orig_idx]),
            "library_norm_source": lib_source[orig_idx],
            "flanking_chr10_mad": fmt_num(flank_mad, 4),
            "pten_mean_copy_ratio": fmt_num(pten_cr, 4),
            "pten_mean_zscore": fmt_num(pten_mean_z, 3),
            "pten_peak_zscore": fmt_num(pten_peak_z, 3),
            "pten_abnormal_bins": abn_bins,
            "total_pten_bins": len(pten_indices),
            "cnv_call": ctype,
            "est_cnv_start": c_start if ctype not in ("NORMAL", "NO_SIGNAL") else "NA",
            "est_cnv_end": c_end if ctype not in ("NORMAL", "NO_SIGNAL") else "NA",
            "est_cnv_size_bp": c_size if ctype not in ("NORMAL", "NO_SIGNAL") else "NA",
        }
        report_rows.append(row_dict)
        if ctype not in ("NORMAL", "NO_SIGNAL"):
            candidate_rows.append(row_dict)

    write_report(all_out, report_rows)
    write_report(cand_out, candidate_rows)

    if args.save_matrix:
        import gzip
        mat_out = Path(f"{args.out_prefix}.copy_ratio_matrix.tsv.gz")
        bin_headers = [f"{args.chrom}:{s}-{e}" for s, e in zip(bin_starts, bin_ends)]
        with gzip.open(mat_out, "wt", encoding="utf-8") as mfh:
            mfh.write("sample_id\t" + "\t".join(bin_headers) + "\n")
            v_ptr = 0
            for orig_idx, (sid, _) in enumerate(sample_bams):
                if valid_mask[orig_idx]:
                    vals_str = "\t".join(f"{v:.4f}" for v in copy_ratio[v_ptr, :])
                    mfh.write(f"{sid}\t{vals_str}\n")
                    v_ptr += 1

    n_no_signal = sum(1 for row in report_rows if row.get("cnv_call") == "NO_SIGNAL")
    print(f"[DONE] Valid samples evaluated: {n_valid}/{total_input}")
    if n_no_signal:
        print(f"[WARN] {n_no_signal} sample(s) had no usable PTEN-bin signal (called NO_SIGNAL).", file=sys.stderr)
    if n_valid < total_input:
        print(f"[INFO] {total_input - n_valid} sample(s) skipped, by reason:", file=sys.stderr)
        print_status_summary(status_list)
        
    print(f"[SUMMARY] Candidate PTEN CNVs detected: {len(candidate_rows)} (written to {cand_out})")


if __name__ == "__main__":
    main()



