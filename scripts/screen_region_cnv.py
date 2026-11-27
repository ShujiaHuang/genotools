#!/usr/bin/env python3
"""
Fast multi-process CNV (microdeletion/microduplication) screener for a target
region of one chromosome, defaulting to PTEN (10q23.31) on GRCh38, across
large-scale low-pass NIPT BAM/CRAM cohorts.

The target locus is a parameter rather than an assumption: --target takes a
samtools/bcftools style region (default chr10:87863625-87971930, PTEN), while the
surrounding scan window (--scan) only supplies the bins used for normalization and
is fitted around the target automatically when left unset.

Uses ONLY pysam for BAM/CRAM and FASTA operations, combined with:
  1) Library-size normalization (mapped reads from the index; CRAI indexes do not
     record per-contig counters, so CRAM input falls back to the number of usable
     reads inside the scan window, reported as library_norm_source=region)
  2) Per-sample GC-content stratification correction (when --ref-fasta is given)
  3) Cohort-wide bin median & robust MAD Z-score normalization
  4) Optional figures (--plot-dir, needs matplotlib): one cohort overview figure plus
     per-sample profiles for the called candidates; --plot-samples narrows that to the
     samples it names ('all' for every sample, 'candidates,<id>' to keep the candidates)

Author: Shujia Huang and Gemini & Qoder
Date: 2026-11-27
"""
import argparse
import csv
import re
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pysam

# Built-in target: PTEN (10q23.31) on GRCh38, used only when --target is omitted. Any
# other locus is one region string away, e.g. --target chr17:41190000-41310000.
DEFAULT_TARGET_REGION = "chr10:87863625-87971930"
# Historical normalization window around the default target. Keeping it as the default
# means a run without --scan reproduces the original PTEN numbers.
DEFAULT_SCAN_START = 80000000
DEFAULT_SCAN_END = 96000000
# Flanking width used to fit the scan window around any other target.
DEFAULT_FLANK_BP = 8_000_000

# samtools/bcftools style interval: 'chr' or 'chr:start-end' (1-based, inclusive).
REGION_RE = re.compile(r"^([^:\s]+)(?::([\d,]+)(?:-([\d,]+))?)?$")

# Half-width of the zoomed panels around the target centre. The scan window is far wider
# than the target (~16 Mb vs ~110 kb for PTEN), so a bin-level event is invisible on a
# full-window axis; the figures therefore show both views.
PLOT_FLANK_BP = 1_000_000

# Copy-ratio gates: a Z-score tail alone is not enough to call a CNV, the mean copy
# ratio of the target bins has to move at least this far as well.
DEL_COPY_RATIO_MAX = 0.80
DUP_COPY_RATIO_MIN = 1.20

# Minimum usable reads per bin when the library size is estimated from the scan
# window (CRAM input, whose .crai index carries no per-contig read counters).
MIN_REGION_READS_PER_BIN = 5


def positive_int(value: str) -> int:
    """argparse type for counts that must be at least 1 (bin/step size, threads)."""
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected an integer, got '{value}'")
    if number < 1:
        raise argparse.ArgumentTypeError(f"must be >= 1, got {number}")
    return number


def non_negative_int(value: str) -> int:
    """argparse type for thresholds where 0 is meaningful (MAPQ, library size)."""
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected an integer, got '{value}'")
    if number < 0:
        raise argparse.ArgumentTypeError(f"must be >= 0, got {number}")
    return number


def positive_float(value: str) -> float:
    """argparse type for thresholds that a zero or negative value would degenerate."""
    try:
        number = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a number, got '{value}'")
    if not number > 0:
        raise argparse.ArgumentTypeError(f"must be > 0, got {number}")
    return number


def fraction(value: str) -> float:
    """argparse type for fractions in (0, 1]."""
    try:
        number = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a number, got '{value}'")
    if not 0 < number <= 1:
        raise argparse.ArgumentTypeError(f"must be in (0, 1], got {number}")
    return number


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Screen CNVs in a target region (default: chr10 PTEN) from NIPT BAM/CRAM files using pysam."
    )
    p.add_argument("--bam-list", required=True, type=Path, help="File with 'sample_id<TAB>bam_path' or one bam_path per line.")
    p.add_argument("-o", "--out-prefix", required=True, type=Path, help="Output prefix for TSV reports.")
    p.add_argument("--ref-fasta", type=Path, default=None, help="GRCh38 reference FASTA (required for CRAM and GC correction).")
    p.add_argument("--target", metavar="REGION", default=None, help="Region to screen, samtools style 'chr:start-end' (default: chr10:87863625-87971930, PTEN).")
    p.add_argument("--scan", metavar="REGION", default=None, help="Region that supplies the bins used for normalization; must contain --target (default: chr10:80000000-96000000 when it covers the target, otherwise target +/- 8 Mb clamped to the contig).")
    p.add_argument("--bin-size", type=positive_int, default=50000, help="Sliding window size in bp (default: 50000).")
    p.add_argument("--step-size", type=positive_int, default=25000, help="Sliding window step size in bp (default: 25000).")
    p.add_argument("--min-mapq", type=non_negative_int, default=30, help="Minimum read mapping quality (default: 30).")
    p.add_argument("--max-n-frac", type=fraction, default=0.20, help="Max allowed 'N' base fraction per bin, in (0, 1] (default: 0.20).")
    p.add_argument("--z-thresh", type=positive_float, default=4.0, help="Absolute Z-score threshold for candidate CNV calling (default: 4.0).")
    p.add_argument("-t", "--threads", type=positive_int, default=16, help="Number of parallel worker processes (default: 16).")
    p.add_argument("--min-library-reads", type=non_negative_int, default=10000, help="Minimum library size of a sample: mapped reads from the index, or usable reads in the scan window for CRAM (default: 10000).")
    p.add_argument("--min-valid-samples", type=positive_int, default=10, help="Minimum number of QC-passing samples required for cohort normalization (default: 10).")
    p.add_argument("--save-matrix", action="store_true", help="Save normalized copy-ratio matrix (.tsv.gz) for downstream segmentation.")
    p.add_argument("--plot-dir", type=Path, default=None, help="Write a cohort overview figure into this directory, plus a PNG for each called candidate (or only the samples named by --plot-samples); needs matplotlib, off by default.")
    p.add_argument("--plot-samples", metavar="ID[,ID...]|all|candidates", default=None,
                   help="Per-sample figures to write when --plot-dir is given: a comma-separated id list draws ONLY those samples, 'all' draws every QC-passing sample, and the token 'candidates' keeps the called samples in the list (default without this flag: candidates only).")
    return p.parse_args()


def parse_region(
    text: str,
    what: str,
    contig_lengths: Optional[Dict[str, int]] = None,
) -> Tuple[str, int, int]:
    """Parse a samtools/bcftools style interval into (chrom, start, end), 1-based inclusive.

    Accepts 'chr:start-end' and a bare 'chr' (the whole contig, which needs its length);
    digit groups may be comma separated, as in 'chr10:87,863,625-87,971,930'.
    """
    match = REGION_RE.match(text.strip())
    if not match:
        sys.exit(f"[ERROR] Cannot parse {what} '{text}': use 'chr:start-end', "
                 f"e.g. chr10:87863625-87971930.")

    chrom, start_text, end_text = match.groups()
    if start_text is None:
        length = (contig_lengths or {}).get(chrom)
        if length is None:
            sys.exit(f"[ERROR] {what} '{text}' has no coordinates and the length of '{chrom}' "
                     f"is unknown. Give an explicit interval, e.g. {chrom}:1-1000000.")
        return chrom, 1, length
    if end_text is None:
        sys.exit(f"[ERROR] {what} '{text}' is missing the end coordinate: use 'chr:start-end'.")

    start, end = int(start_text.replace(",", "")), int(end_text.replace(",", ""))
    if start < 1 or end < start:
        sys.exit(f"[ERROR] Invalid {what} '{text}': need 1 <= start <= end.")
    return chrom, start, end


def resolve_target_region(
    target: Optional[str],
    contig_lengths: Optional[Dict[str, int]] = None,
) -> Tuple[str, int, int]:
    """Chromosome and 1-based interval to screen; the built-in default is PTEN on chr10."""
    region_text = target or DEFAULT_TARGET_REGION
    chrom, start, end = parse_region(region_text, "--target", contig_lengths)

    if contig_lengths and chrom not in contig_lengths:
        known = ", ".join(list(contig_lengths)[:5])
        sys.exit(f"[ERROR] Chromosome '{chrom}' from --target '{region_text}' is not in the "
                 f"alignment header. Use the contig naming of the files, e.g. one of: {known} ...")

    length = (contig_lengths or {}).get(chrom)
    if length is not None and end > length:
        sys.exit(f"[ERROR] --target '{region_text}' ends beyond the length of '{chrom}' "
                 f"({length} bp).")
    return chrom, start, end


def resolve_scan_window(
    chrom: str,
    target_start: int,
    target_end: int,
    scan: Optional[str] = None,
    contig_lengths: Optional[Dict[str, int]] = None,
    flank_bp: int = DEFAULT_FLANK_BP,
) -> Tuple[int, int, str]:
    """Normalization window that must contain the target region.

    Without --scan the window is the historical 80-96 Mb PTEN window when it covers the
    target, and otherwise a +/- flank_bp flank around it (clamped to the contig), so any
    locus works without extra flags. Returns (start, end, source) where source is
    'explicit', 'default' or 'fitted'.
    """
    if scan is None:
        length = (contig_lengths or {}).get(chrom)
        if DEFAULT_SCAN_START <= target_start and target_end <= DEFAULT_SCAN_END:
            # The historical window still has to be clamped: on a contig shorter than 96 Mb
            # (chr16 is ~90.3 Mb) it would run past the chromosome end.
            scan_end = DEFAULT_SCAN_END if length is None else min(DEFAULT_SCAN_END, length)
            return DEFAULT_SCAN_START, scan_end, "default"

        scan_start = max(1, target_start - flank_bp)
        scan_end = target_end + flank_bp
        length = (contig_lengths or {}).get(chrom)
        if length is not None:
            scan_start = min(scan_start, length)
            scan_end = min(scan_end, length)
        source = "fitted"
    else:
        scan_chrom, scan_start, scan_end = parse_region(scan, "--scan", contig_lengths)
        if scan_chrom != chrom:
            sys.exit(f"[ERROR] --scan '{scan}' is on '{scan_chrom}' but --target is on '{chrom}'; "
                     f"both regions must be on the same chromosome.")
        length = (contig_lengths or {}).get(chrom)
        if length is not None and scan_end > length:
            sys.exit(f"[ERROR] --scan '{scan}' ends beyond the length of '{chrom}' ({length} bp).")
        source = "explicit"

    if target_start < scan_start or target_end > scan_end:
        sys.exit(f"[ERROR] Target region {chrom}:{target_start}-{target_end} is not inside the "
                 f"scan window {chrom}:{scan_start}-{scan_end}. Widen --scan: the window only "
                 f"supplies the bins used for normalization.")
        
    return scan_start, scan_end, source


def unique_sample_id(sid: str, used: set) -> str:
    """Return a sample id no other entry of this list already uses.

    A duplicate used to be renamed to '<id>_<sample count>', but that candidate can be
    taken as well (a list holding 'A_2', 'A', 'A' produced 'A_2' twice), and two rows
    sharing one id silently merge in every downstream table. So keep counting until the
    id is really free, and say so, because reports and figures are keyed by this name.
    """
    if sid not in used:
        used.add(sid)
        return sid

    suffix = 1
    while f"{sid}_{suffix}" in used:
        suffix += 1
    renamed = f"{sid}_{suffix}"
    print(f"[WARN] Duplicate sample id '{sid}': renamed to '{renamed}' for this run. Use "
          f"unique ids in --bam-list to keep the reports unambiguous.", file=sys.stderr)
    used.add(renamed)
    return renamed


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
            
            samples.append((unique_sample_id(sid, seen), bpath))
            
    if not samples:
        sys.exit(f"[ERROR] No valid BAM/CRAM entries found in {bam_list_path}")
        
    return samples


def build_bins(
    chrom: str,
    target_start: int,
    target_end: int,
    scan_start: int,
    scan_end: int,
    bin_size: int,
    step_size: int,
    ref_fasta: Optional[Path],
    max_n_frac: float,
    contig_lengths: Optional[Dict[str, int]] = None,
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], np.ndarray]:
    """
    Generate valid sliding bins across [scan_start, scan_end], compute GC content
    via pysam.FastaFile if provided, and return the indices of the bins overlapping
    the target region [target_start, target_end].

    contig_lengths (from the alignment header) is optional; when given, the FASTA
    length of the target contig is checked against it so a GRCh37 BAM paired with a
    GRCh38 reference fails loudly instead of silently scanning the wrong bases.
    """
    starts: List[int] = []
    ends: List[int] = []
    gc_vals: List[float] = []

    fasta = pysam.FastaFile(str(ref_fasta)) if ref_fasta else None
    if fasta and chrom not in fasta.references:
        fasta.close()
        sys.exit(f"[ERROR] Chromosome '{chrom}' not found in --ref-fasta '{ref_fasta}'. "
                 "Use a reference that contains the target contig, or drop --ref-fasta.")

    if fasta and contig_lengths:
        header_length = contig_lengths.get(chrom)
        fasta_length = fasta.get_reference_length(chrom)
        if header_length is not None and fasta_length != header_length:
            fasta.close()
            sys.exit(f"[ERROR] '{chrom}' is {header_length} bp in the alignment header but "
                     f"{fasta_length} bp in --ref-fasta '{ref_fasta}': the reads and the "
                     f"reference come from different builds (e.g. GRCh37 vs GRCh38), so the "
                     f"nucleotides and coordinates would not match.")

    for pos in range(scan_start, scan_end - bin_size + 1, step_size):
        b_start = pos
        b_end = pos + bin_size - 1
        if fasta:
            seq = fasta.fetch(chrom, b_start - 1, b_end).upper()
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

    # Identify bins overlapping the target locus
    target_mask = (bin_ends >= target_start) & (bin_starts <= target_end)
    target_indices = np.where(target_mask)[0]
    if len(target_indices) == 0:
        sys.exit(f"[ERROR] No valid bin overlaps the target region {chrom}:{target_start}-"
                 f"{target_end}. Check --target and --max-n-frac "
                 f"(bins with too many N bases are dropped).")

    return bin_starts, bin_ends, gc_array, target_indices


def index_library_depth(bam: pysam.AlignmentFile) -> Optional[int]:
    """
    Mapped reads recorded by the BAM/CRAM index.

    Only a per-sample constant factor matters here: the downstream cohort bin median
    and per-sample median re-centering strip any scale, so counting every contig is
    equivalent to counting autosomes only and survives non-standard contig naming.

    Returns None when the index carries no per-contig counters: CRAI indexes do not
    store them, so pysam reports 0 mapped reads for every contig of a CRAM file
    (see pysam.AlignmentFile.get_index_statistics). Callers must fall back to a
    region-based library size in that case.
    """
    try:
        stats = bam.get_index_statistics()
        """
        stats 的返回值格式：
        chr1   248956422   1599044   0
        │        │          │      └─ unmapped(落在该 contig 上的未比对 reads, 通常恒为 0)
        │        │          └──────── mapped: chr1 上 1,599,044 条比对 reads
        │        └─────────────────── chr1 长度 248,956,422 → GRCh38
        └──────────────────────────── contig 名
        """
    except (ValueError, NotImplementedError):
        return None

    total = sum(stat.mapped for stat in stats)
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
    where library_source is 'index' (mapped reads from the index) or 'region' (usable
    reads in the scan window, used when the index has no counters, e.g. CRAM), and 'NA'
    when nothing could be measured (unreadable file, missing index).
    library_reads is always an int, 0 when neither source could determine it, so that
    callers can store it in an integer array.
    """
    idx, sid, bam_path = task
    library_reads = 0
    library_source = "index"

    try:
        bam = pysam.AlignmentFile(bam_path, reference_filename=ref_fasta_str)
    except Exception as exc:
        # Nothing was measured for this sample: report no source instead of claiming 'index'.
        return (idx, sid, None, 0, "NA", f"OpenError: {exc}")

    try:
        if not bam.has_index():
            return (idx, sid, None, 0, "NA",
                    "NoIndex: missing .bai/.csi/.crai index (run samtools index)")

        # 1. Library size from the index (fast, exact for BAM). CRAI indexes record no
        #    per-contig counters, so this is None for CRAM: keep library_reads at 0 and
        #    let the region-based fallback below supply it instead. Never assign None to
        #    library_reads, or the early returns leak it to the caller.
        indexed_depth = index_library_depth(bam)
        if indexed_depth is not None:
            library_reads = indexed_depth
            if indexed_depth < min_library_reads:
                return (idx, sid, None, library_reads, library_source,
                        f"LowLibrary: {indexed_depth} indexed mapped reads < {min_library_reads}")

        # 2. Fetch reads across the target region
        read_positions: List[int] = []
        for aln in bam.fetch(chrom, max(0, scan_start - 1), scan_end):
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

        if indexed_depth is None:
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

        # Nothing measured means no source either: a 0-read row must not claim 'index'.
        return (idx, sid, None, library_reads, library_source if library_reads else "NA", msg)
    
    finally:
        bam.close()


def estimate_cnv_span(
    z_row: np.ndarray,
    ratio_row: np.ndarray,
    bin_starts: np.ndarray,
    bin_ends: np.ndarray,
    target_indices: np.ndarray,
    z_thresh: float,
) -> Tuple[str, int, int, int]:
    """Expand contiguous abnormal bins outward from the target region to estimate CNV boundaries."""
    target_z = z_row[target_indices]
    finite = np.isfinite(target_z)
    if finite.any():
        # Start from the strongest target bin, ignoring bins without usable signal
        center_idx = int(target_indices[finite][int(np.argmax(np.abs(target_z[finite])))])
    else:
        center_idx = int(target_indices[0])
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


def load_pyplot() -> Any:
    """Import matplotlib lazily: the screener must run without it unless --plot-dir is used."""
    try:
        import matplotlib
        matplotlib.use("Agg")  # no display needed, so it also works on a cluster
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("[ERROR] --plot-dir needs matplotlib: pip install matplotlib")
    return plt


def safe_filename(name: str) -> str:
    """Keep figure filenames portable when a sample id carries punctuation or separators."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def unique_figure_name(base: str, used: set) -> str:
    """Return a file name no other sample of this run will use.

    safe_filename() is not injective ('A/B' and 'A_B' both become 'A_B'), so a collision
    has to be resolved here instead of letting one figure overwrite another.
    """
    name = f"{base}.png"
    suffix = 1
    while name in used or name == "cohort_overview.png":
        suffix += 1
        name = f"{base}_{suffix}.png"
    used.add(name)
    return name


def plot_sample(
    plt: Any,
    out_path: Path,
    row: Dict[str, object],
    bin_starts: np.ndarray,
    bin_ends: np.ndarray,
    cr_row: np.ndarray,
    z_row: np.ndarray,
    z_thresh: float,
    chrom: str,
    target_start: int,
    target_end: int,
) -> None:
    """Write one three-panel profile: whole window, zoomed copy ratio, zoomed Z-score.

    The upper panel shows the whole normalization window (arm-level shifts), the lower two
    zoom to target +/- PLOT_FLANK_BP where individual bins are wide enough to read.
    """
    centres_mb = (bin_starts + bin_ends) / 2e6
    target_lo_mb, target_hi_mb = target_start / 1e6, target_end / 1e6

    centre = (target_start + target_end) // 2
    zoom_lo = max(int(bin_starts[0]), centre - PLOT_FLANK_BP)
    zoom_hi = min(int(bin_ends[-1]), centre + PLOT_FLANK_BP)
    in_zoom = (centres_mb * 1e6 >= zoom_lo) & (centres_mb * 1e6 <= zoom_hi)
    zoom_x = centres_mb[in_zoom]

    call = str(row["cnv_call"])
    event_color = {"DEL": "#d62728", "DUP": "#1f77b4"}.get(call, "#7f7f7f")
    has_event = call in ("DEL", "DUP")
    event_mb = None
    est_start = est_end = 0
    if has_event:
        est_start, est_end = int(row["est_cnv_start"]), int(row["est_cnv_end"])
        event_mb = (est_start / 1e6, est_end / 1e6)

    fig, (ax_all, ax_cr, ax_z) = plt.subplots(
        3, 1, figsize=(10, 8.5), gridspec_kw={"height_ratios": [1.0, 1.1, 0.9]}
    )
    fig.suptitle(
        f"{row['sample_id']}  {call}  |  {chrom}:{target_start}-{target_end}  |  "
        f"CR {row['target_mean_copy_ratio']}  mean Z {row['target_mean_zscore']}  "
        f"peak Z {row['target_peak_zscore']}  |  {row['library_reads']} reads "
        f"({row['library_norm_source']})",
        fontsize=9,
    )

    # 1) whole window: arm-level shifts are visible, the target itself is a hairline
    ax_all.plot(centres_mb, cr_row, lw=0.5, color="#444444")
    ax_all.axhline(1.0, color="black", lw=0.8)
    for gate in (DEL_COPY_RATIO_MAX, DUP_COPY_RATIO_MIN):
        ax_all.axhline(gate, color="#d62728", lw=0.6, ls="--")
    for edge in (target_lo_mb, target_hi_mb):
        ax_all.axvline(edge, color="orange", lw=0.8)
    ax_all.set_ylabel("copy ratio")
    ax_all.set_title(f"scan window {chrom}:{int(bin_starts[0])}-{int(bin_ends[-1])}",
                     fontsize=8, loc="left")

    def draw_zoom(ax: Any, values: np.ndarray, color: str) -> None:
        ax.plot(zoom_x, values[in_zoom], marker="o", ms=2.5, lw=0.9, color=color)
        ax.axvspan(target_lo_mb, target_hi_mb, color="orange", alpha=0.20, lw=0)
        if event_mb is not None:
            ax.axvspan(*event_mb, color=event_color, alpha=0.15, lw=0)
        ax.set_xlim(zoom_lo / 1e6, zoom_hi / 1e6)

    # 2) zoomed copy ratio around the target
    draw_zoom(ax_cr, cr_row, "#444444")
    ax_cr.axhline(1.0, color="black", lw=0.8)
    for gate in (DEL_COPY_RATIO_MAX, DUP_COPY_RATIO_MIN):
        ax_cr.axhline(gate, color="#d62728", lw=0.7, ls="--")
    if event_mb is not None:
        ax_cr.annotate(
            f"{call} {est_start:,}-{est_end:,} ({int(row['est_cnv_size_bp']) / 1000:.1f} kb)",
            xy=(event_mb[0], 1.0), xytext=(4, 8), textcoords="offset points",
            fontsize=8, color=event_color,
        )
    ax_cr.set_ylabel("copy ratio")
    ax_cr.set_title(f"target +/- {PLOT_FLANK_BP / 1e6:g} Mb", fontsize=8, loc="left")

    # 3) zoomed Z-score of the same bins
    draw_zoom(ax_z, z_row, "#1f77b4")
    ax_z.axhline(0.0, color="black", lw=0.8)
    for thr in (-z_thresh, z_thresh):
        ax_z.axhline(thr, color="#d62728", lw=0.7, ls="--")
    ax_z.set_ylabel("Z-score")
    ax_z.set_xlabel(f"position on {chrom} (Mb)")
    ax_z.set_title(f"threshold |Z| >= {z_thresh:g}", fontsize=8, loc="left")

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_cohort_overview(
    plt: Any,
    out_path: Path,
    sids: List[str],
    calls: List[str],
    target_cr: np.ndarray,
    flank_mad: np.ndarray,
    matrix: np.ndarray,
    bin_starts: np.ndarray,
    bin_ends: np.ndarray,
    chrom: str,
    target_start: int,
    target_end: int,
    title: str,
) -> None:
    """One figure for the whole cohort: which samples are worth a closer look, and why.

    Rows of the heatmap are sorted by target copy ratio, so a real event appears as a
    whole block of consistent colour while a noisy sample shows up as speckle. The ranked
    curve and the copy-ratio-vs-noise scatter then separate true calls from artefacts.
    """
    n = len(sids)
    target_cr = np.asarray(target_cr, dtype=np.float64)
    flank_mad = np.asarray(flank_mad, dtype=np.float64)
    # NaN target copy ratios (NO_SIGNAL) sort to the top instead of breaking the order.
    order = np.argsort(np.where(np.isfinite(target_cr), target_cr, np.inf), kind="stable")
    cr_sorted = target_cr[order]
    call_sorted = [calls[i] for i in order]
    sid_sorted = [sids[i] for i in order]

    call_colours = {"DEL": "#d62728", "DUP": "#1f77b4"}
    colours = [call_colours.get(c, "#999999") for c in calls]
    n_cand = sum(1 for c in calls if c in call_colours)

    centres_mb = (bin_starts + bin_ends) / 2e6
    n_cols = len(centres_mb)
    fig = plt.figure(figsize=(13, 8.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.0])
    ax_hm, ax_rank, ax_qc = fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    fig.suptitle(f"{title}  |  {chrom}:{target_start}-{target_end}  |  "
                 f"{n} QC-passing sample(s), {n_cand} candidate(s)", fontsize=10)

    # 1) heatmap: samples x bins, ordered by target copy ratio
    cmap = plt.get_cmap("RdBu_r")
    if hasattr(cmap, "with_extremes"):  # matplotlib >= 3.4; keep working on older clusters
        cmap = cmap.with_extremes(bad="#dddddd")
    else:
        cmap = cmap.copy()
        cmap.set_bad("#dddddd")
    # imshow spreads the columns evenly over the whole span, which silently moves a colour
    # block when --max-n-frac dropped bins and left gaps. Only even grids may use it.
    uniform = n_cols < 3 or bool(np.allclose(np.diff(centres_mb), np.diff(centres_mb)[0], rtol=1e-6))
    masked = np.ma.masked_invalid(matrix[order, :])
    if uniform:
        half_mb = 0.0 if n_cols == 1 else float(np.diff(centres_mb)[0]) / 2
        im = ax_hm.imshow(
            masked, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0.5, vmax=1.5,
            extent=(centres_mb[0] - half_mb, centres_mb[-1] + half_mb, n, 0),
        )
    else:
        # One column per surviving bin, placed at its real position; gaps stretch instead of
        # shifting every later column away from the truth.
        edges = np.empty(n_cols + 1, dtype=np.float64)
        edges[1:-1] = (centres_mb[:-1] + centres_mb[1:]) / 2
        edges[0] = centres_mb[0] - (centres_mb[1] - centres_mb[0]) / 2
        edges[-1] = centres_mb[-1] + (centres_mb[-1] - centres_mb[-2]) / 2
        im = ax_hm.pcolormesh(edges, np.arange(n + 1), masked, cmap=cmap, vmin=0.5, vmax=1.5)
        ax_hm.set_ylim(n, 0)
    for edge in (target_start / 1e6, target_end / 1e6):
        ax_hm.axvline(edge, color="black", lw=0.8, ls="--")
    ax_hm.set_xlabel(f"position on {chrom} (Mb)")
    ax_hm.set_ylabel(f"samples, sorted by target copy ratio (n={n})")
    ax_hm.set_title("copy ratio; the target region sits between the dashed lines",
                    fontsize=9, loc="left")
    cbar = fig.colorbar(im, ax=ax_hm, fraction=0.03, pad=0.01)
    cbar.set_label("copy ratio (clipped at 0.5/1.5)")

    # Small cohorts are labelled by id; for large ones the rank panel below marks the
    # candidates instead of turning the y axis into an unreadable wall of labels.
    if n <= 40:
        ax_hm.set_yticks(np.arange(n) + 0.5)
        ax_hm.set_yticklabels(sid_sorted, fontsize=6)
    else:
        ax_hm.set_yticks([])

    # Label the candidates next to their dots. Many candidates in a big cohort sit in the
    # same rank tail, so keep every label only when they stay apart; otherwise always label
    # the most extreme DEL/DUP and let the candidates table carry the full id list.
    ranked_calls = [(pos, sid, call) for pos, (sid, call) in
                    enumerate(zip(sid_sorted, call_sorted)) if call in call_colours]
    min_gap = max(3.0, 0.02 * n)
    extremes = set()
    for call in ("DEL", "DUP"):  # the most extreme candidate of each direction
        positions = [pos for pos, _sid, c in ranked_calls if c == call]
        if positions:
            extremes.add(positions[0] if call == "DEL" else positions[-1])
    labelled: set = set()
    last_pos = -min_gap - 1
    for pos, _sid, _call in ranked_calls:
        if pos in extremes or pos - last_pos >= min_gap:
            labelled.add(pos)
            last_pos = pos

    # 2) ranked target copy ratio: how far the outlier tail reaches beyond the cohort bulk
    ax_rank.plot(np.arange(1, n + 1), cr_sorted, marker="o", ms=3, lw=0.6, color="#444444")
    ax_rank.axhline(1.0, color="black", lw=0.8)
    for gate in (DEL_COPY_RATIO_MAX, DUP_COPY_RATIO_MIN):
        ax_rank.axhline(gate, color="#d62728", lw=0.7, ls="--")
    stack = {"DEL": 0, "DUP": 0}  # stagger labels when two dots sit at nearly the same height
    for rank, (call, value, sid) in enumerate(zip(call_sorted, cr_sorted, sid_sorted), start=1):
        if call in call_colours:
            ax_rank.plot(rank, value, marker="o", ms=5, color=call_colours[call])
            if rank - 1 in labelled and np.isfinite(value):
                # keep the label inside the axes at both ends of the rank axis
                ha = "left" if rank <= 2 else "right" if rank >= n - 1 else "center"
                # label towards the inside of the plot: the tails sit on the axis edges
                dy = (-11 - 10 * stack["DUP"]) if call == "DUP" else (11 + 10 * stack["DEL"])
                stack[call] += 1
                ax_rank.annotate(
                    f"{sid} ({call}, {value:.2f})", xy=(rank, value),
                    xytext=(0, dy), textcoords="offset points",
                    ha=ha, fontsize=6, color=call_colours[call], annotation_clip=True,
                )
    if len(labelled) < n_cand:
        rank_title = (f"copy ratio by rank: {n_cand} candidates, extremes labelled "
                      f"(full list in the candidates table)")
    else:
        rank_title = "copy ratio by rank: the tails are the candidates"
    ax_rank.set_xlabel("samples (sorted)")
    ax_rank.set_ylabel("target mean copy ratio")
    ax_rank.set_title(rank_title, fontsize=9, loc="left")
    # leave room above/below the gates so the labels of the outliers stay visible
    finite_cr = cr_sorted[np.isfinite(cr_sorted)]
    if finite_cr.size:
        pad = max(0.05, 0.1 * float(np.ptp(finite_cr)))
        ax_rank.set_ylim(float(finite_cr.min()) - pad, float(finite_cr.max()) + pad)

    # 3) quality context: an outlier copy ratio inside a noisy sample is not a CNV
    ax_qc.scatter(flank_mad, target_cr, s=12, c=colours, linewidths=0)
    ax_qc.axhline(1.0, color="black", lw=0.8)
    for gate in (DEL_COPY_RATIO_MAX, DUP_COPY_RATIO_MIN):
        ax_qc.axhline(gate, color="#d62728", lw=0.7, ls="--")
    finite_mad = flank_mad[np.isfinite(flank_mad) & (flank_mad > 0)]
    # A log axis only helps when the cohort's noise really spans orders of magnitude.
    if finite_mad.size and float(finite_mad.max()) > 20 * float(finite_mad.min()):
        ax_qc.set_xscale("log")
    ax_qc.set_xlabel("flanking MAD of the copy ratio")
    ax_qc.set_ylabel("target mean copy ratio")
    ax_qc.set_title("significance vs noise: high and far right = noisy, not a CNV",
                    fontsize=9, loc="left")
    for label, colour in (("NORMAL", "#999999"), ("DEL", "#d62728"), ("DUP", "#1f77b4")):
        ax_qc.scatter([], [], s=14, c=colour, label=label)
    ax_qc.legend(loc="upper right", fontsize=6, frameon=False)
    if finite_cr.size:
        ax_qc.set_ylim(float(finite_cr.min()) - pad, float(finite_cr.max()) + pad)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


ALL_COLUMNS = [
    "sample_id", "status", "library_reads", "library_norm_source", "flanking_mad",
    "target_mean_copy_ratio", "target_mean_zscore", "target_peak_zscore",
    "target_abnormal_bins", "total_target_bins", "cnv_call",
    "est_cnv_start", "est_cnv_end", "est_cnv_size_bp",
]


def qc_failed_row(
    sid: str, status: str, library_reads: int, library_source: str, n_target_bins: int, cnv_call: str = "FAILED_QC"
) -> Dict[str, object]:
    """Per-sample report row for a sample that was not evaluated against the cohort."""
    return {
        "sample_id": sid, 
        "status": status,
        "library_reads": library_reads, 
        "library_norm_source": library_source,
        "flanking_mad": "NA",
        "target_mean_copy_ratio": "NA",
        "target_mean_zscore": "NA",
        "target_peak_zscore": "NA",
        "target_abnormal_bins": "NA",
        "total_target_bins": n_target_bins,
        "cnv_call": cnv_call,
        "est_cnv_start": "NA",
        "est_cnv_end": "NA",
        "est_cnv_size_bp": "NA"
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

    ref_fasta_str = str(args.ref_fasta) if args.ref_fasta else None

    # Peek at the first readable sample once: its contig names and lengths validate --target
    # and keep a fitted --scan window inside the chromosome, so a mismatch fails up front
    # instead of turning the whole cohort into ReadError. A few samples are tried because
    # one broken first entry must not blind the checks for the whole cohort.
    contig_lengths: Dict[str, int] = {}
    for _, probe_path in sample_bams[:5]:
        try:
            with pysam.AlignmentFile(probe_path, reference_filename=ref_fasta_str) as probe:
                contig_lengths = dict(zip(probe.references, probe.lengths))
            break
        
        except (OSError, ValueError):
            # An unreadable sample is reported by its worker; try the next one for the header.
            continue

    chrom, target_start, target_end = resolve_target_region(args.target, contig_lengths)
    scan_start, scan_end, window_source = resolve_scan_window(
        chrom, target_start, target_end, args.scan, contig_lengths
    )
    if scan_end - scan_start + 1 < args.bin_size:
        sys.exit(f"[ERROR] Scan window {chrom}:{scan_start}-{scan_end} is shorter than "
                 f"--bin-size {args.bin_size}: no bin can be built.")
    
    print(f"[INFO] Target {chrom}:{target_start}-{target_end}; normalization window "
          f"{chrom}:{scan_start}-{scan_end} ({window_source}).", file=sys.stderr)

    # Figures are opt-in: import matplotlib here, before the worker pool, so a missing
    # package fails in seconds instead of after the whole cohort has been processed.
    plt_backend = load_pyplot() if args.plot_dir is not None else None
    plot_dir = Path(args.plot_dir) if args.plot_dir is not None else None
    if plot_dir is not None:
        plot_dir.mkdir(parents=True, exist_ok=True)

    # The overview figure answers "is anything wrong at all?"; per-sample figures are for
    # following up on just the samples the overview singles out. Naming samples on the
    # command line therefore REPLACES the default (candidates), it does not add to it:
    # '--plot-samples N0123,N0456' used to keep drawing every candidate as well, which on
    # a noisy cohort means hundreds of figures nobody asked for.
    plot_all = False
    draw_candidates = True  # the policy when --plot-samples is absent
    plot_requested: set = set()
    if args.plot_samples:
        tokens = [token.strip() for token in args.plot_samples.split(",") if token.strip()]
        flags = {token.lower() for token in tokens}
        plot_all = "all" in flags
        draw_candidates = "candidates" in flags
        plot_requested = {token for token in tokens if token.lower() not in ("all", "candidates")}
        if not tokens:
            draw_candidates = False
        if plt_backend is None:
            print("[WARN] --plot-samples only applies together with --plot-dir.", file=sys.stderr)
    if plt_backend is not None:
        wanted = []
        if plot_all:
            wanted.append("every QC-passing sample")
        else:
            if plot_requested:
                wanted.append("the requested sample(s): " + ", ".join(sorted(plot_requested)))
            if draw_candidates:
                wanted.append("the called candidates")
        print(f"[INFO] Per-sample figures: {'; '.join(wanted) or 'none'}"
              f" (plus the cohort overview).", file=sys.stderr)

    bin_starts, bin_ends, gc_array, target_indices = build_bins(
        chrom, target_start, target_end, scan_start, scan_end,
        args.bin_size, args.step_size, args.ref_fasta, args.max_n_frac, contig_lengths,
    )
    n_bins = len(bin_starts)
    print(f"[INFO] Constructed {n_bins} valid bins ({len(target_indices)} overlapping the target region).",
          file=sys.stderr)
    target_bp = target_end - target_start + 1
    if target_bp < args.bin_size:
        print(f"[WARN] Target region is {target_bp} bp, shorter than --bin-size {args.bin_size} bp: "
              f"each target bin also covers some flanking sequence. Lower --bin-size/--step-size "
              f"for small loci.", file=sys.stderr)
    if n_bins < 50:
        print(f"[WARN] Only {n_bins} bin(s) in the normalization window: the cohort median/MAD "
              f"are noisy. Widen --scan.", file=sys.stderr)

    # Discretize GC content into 2% strata (0..50) for robust per-sample median correction
    gc_stratum_ids = np.floor(gc_array * 50).astype(np.int32) if gc_array is not None else None

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
                chrom,
                scan_start,
                scan_end,
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
                sid, status_list[i], int(lib_reads[i]), lib_source[i], len(target_indices),
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

    # Re-center each sample by its own median copy ratio across the normalization window
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

    # Evaluate target-region metrics per valid sample
    candidate_rows: List[Dict[str, object]] = []
    report_rows: List[Dict[str, object]] = []
    n_plots = 0
    n_plot_errors = 0
    first_plot_error: Optional[Exception] = None
    # Collected for the single cohort overview figure (one row per QC-passing sample).
    ov_sids: List[str] = []
    ov_calls: List[str] = []
    ov_target_cr: List[float] = []
    ov_flank_mad: List[float] = []
    ov_rows: List[int] = []
    plotted_ids: set = set()
    used_figure_names: set = set()
    written_figures: set = set()

    v_ptr = 0
    for orig_idx, (sid, _) in enumerate(sample_bams):
        if not valid_mask[orig_idx]:
            report_rows.append(qc_failed_row(
                sid, status_list[orig_idx], int(lib_reads[orig_idx]), lib_source[orig_idx], len(target_indices)
            ))
            continue

        cr_row = copy_ratio[v_ptr, :]
        z_row = z_matrix[v_ptr, :]
        v_ptr += 1

        flank_mad = float(np.nanmedian(np.abs(cr_row - 1.0)))
        target_cr = float(np.nanmean(cr_row[target_indices]))
        target_z_vals = z_row[target_indices]
        finite_z = np.isfinite(target_z_vals)
        if finite_z.any():
            abn_bins = int(np.sum(np.abs(target_z_vals[finite_z]) >= args.z_thresh))
            target_mean_z = float(np.mean(target_z_vals[finite_z]))
            target_peak_z = float(target_z_vals[finite_z][int(np.argmax(np.abs(target_z_vals[finite_z])))])
        else:
            abn_bins = 0
            target_mean_z = target_peak_z = float("nan")

        # Call a candidate CNV when >=2 target bins exceed the Z threshold, or when the mean
        # Z does, and the mean copy ratio confirms the same direction.
        is_del = (target_mean_z <= -args.z_thresh or (abn_bins >= 2 and target_peak_z <= -args.z_thresh)) \
            and target_cr <= DEL_COPY_RATIO_MAX
        is_dup = (target_mean_z >= args.z_thresh or (abn_bins >= 2 and target_peak_z >= args.z_thresh)) \
            and target_cr >= DUP_COPY_RATIO_MIN

        if is_del or is_dup:
            ctype, c_start, c_end, c_size = estimate_cnv_span(
                z_row, cr_row, bin_starts, bin_ends, target_indices, args.z_thresh
            )
        elif not np.isfinite(target_cr):
            # No usable signal in the target bins: report it instead of calling NORMAL
            ctype, c_start, c_end, c_size = "NO_SIGNAL", 0, 0, 0
        else:
            ctype, c_start, c_end, c_size = "NORMAL", 0, 0, 0

        row_dict: Dict[str, object] = {
            "sample_id": sid,
            "status": "OK",
            "library_reads": int(lib_reads[orig_idx]),
            "library_norm_source": lib_source[orig_idx],
            "flanking_mad": fmt_num(flank_mad, 4),
            "target_mean_copy_ratio": fmt_num(target_cr, 4),
            "target_mean_zscore": fmt_num(target_mean_z, 3),
            "target_peak_zscore": fmt_num(target_peak_z, 3),
            "target_abnormal_bins": abn_bins,
            "total_target_bins": len(target_indices),
            "cnv_call": ctype,
            "est_cnv_start": c_start if ctype not in ("NORMAL", "NO_SIGNAL") else "NA",
            "est_cnv_end": c_end if ctype not in ("NORMAL", "NO_SIGNAL") else "NA",
            "est_cnv_size_bp": c_size if ctype not in ("NORMAL", "NO_SIGNAL") else "NA",
        }
        report_rows.append(row_dict)
        if ctype not in ("NORMAL", "NO_SIGNAL"):
            candidate_rows.append(row_dict)
        ov_sids.append(sid)
        ov_calls.append(ctype)
        ov_target_cr.append(target_cr)
        ov_flank_mad.append(flank_mad)
        ov_rows.append(v_ptr - 1)

        if plt_backend is not None and (
            plot_all
            or sid in plot_requested
            or (draw_candidates and ctype not in ("NORMAL", "NO_SIGNAL"))
        ):
            figure_name = unique_figure_name(safe_filename(sid), used_figure_names)
            try:
                plot_sample(
                    plt_backend, plot_dir / figure_name, row_dict,
                    bin_starts, bin_ends, cr_row, z_row, args.z_thresh,
                    chrom, target_start, target_end,
                )
                n_plots += 1
                plotted_ids.add(sid)
                written_figures.add(figure_name)
            except Exception as exc:  # a figure must never abort a finished analysis
                n_plot_errors += 1
                if first_plot_error is None:
                    first_plot_error = exc

    write_report(all_out, report_rows)
    write_report(cand_out, candidate_rows)

    if plot_dir is not None:
        overview_path = plot_dir / "cohort_overview.png"
        try:
            plot_cohort_overview(
                plt_backend, overview_path,
                ov_sids, ov_calls, np.array(ov_target_cr), np.array(ov_flank_mad),
                copy_ratio[np.array(ov_rows), :], bin_starts, bin_ends,
                chrom, target_start, target_end, str(args.out_prefix.name),
            )
            print(f"[INFO] Wrote cohort overview {overview_path}", file=sys.stderr)
            written_figures.add(overview_path.name)
        except Exception as exc:  # same rule as the per-sample figures
            n_plot_errors += 1
            if first_plot_error is None:
                first_plot_error = exc
        print(f"[INFO] Wrote {n_plots} per-sample figure(s) to {plot_dir}", file=sys.stderr)
        if n_plot_errors:
            print(f"[WARN] {n_plot_errors} figure(s) could not be written (first error: "
                  f"{type(first_plot_error).__name__}: {first_plot_error})", file=sys.stderr)
        if plot_requested:
            missing = sorted(plot_requested - plotted_ids)
            if missing:
                print(f"[WARN] No figure for {', '.join(missing)}: unknown sample id or the "
                      f"sample did not pass QC.", file=sys.stderr)
        stale = sorted(p.name for p in plot_dir.glob("*.png") if p.name not in written_figures)
        if stale:
            shown = ", ".join(stale[:5]) + (" ..." if len(stale) > 5 else "")
            print(f"[INFO] {len(stale)} PNG(s) in {plot_dir} were not written by this run "
                  f"(kept, but they may be left over from an older run): {shown}", file=sys.stderr)

    if args.save_matrix:
        import gzip
        mat_out = Path(f"{args.out_prefix}.copy_ratio_matrix.tsv.gz")
        bin_headers = [f"{chrom}:{s}-{e}" for s, e in zip(bin_starts, bin_ends)]
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
    # Make the library-size source visible: a window-based estimate is expected for CRAM
    # (.crai stores no read counters) and does not bias copy ratios, but it must never
    # look like a silently zero-length library.
    n_region = sum(1 for src in lib_source if src == "region")
    if n_region:
        print(f"[INFO] {n_region}/{total_input} sample(s) used a window-based library size "
              f"(no per-contig counters in the index, e.g. CRAM/.crai). Per-sample scale "
              f"cancels out in the cohort/sample median normalization.", file=sys.stderr)
    if n_no_signal:
        print(f"[WARN] {n_no_signal} sample(s) had no usable target-region signal (called NO_SIGNAL).", file=sys.stderr)
    if n_valid < total_input:
        print(f"[INFO] {total_input - n_valid} sample(s) skipped, by reason:", file=sys.stderr)
        print_status_summary(status_list)
        
    print(f"[SUMMARY] Candidate CNVs in {chrom}:{target_start}-{target_end}: "
          f"{len(candidate_rows)} (written to {cand_out})")


if __name__ == "__main__":
    main()



