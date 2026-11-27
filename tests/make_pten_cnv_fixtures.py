#!/usr/bin/env python3
"""
Deterministic fixture generator for tests/test_screen_pten_cnv.py.

Builds a tiny synthetic BAM cohort that exercises scripts/screen_pten_cnv.py:
  * 10 normal samples, 1 sample with a 50% deletion, 1 with a 100% duplication
    inside a 20 kb window that overlaps the (test-shrunk) PTEN core locus
  * a 600 kb three-contig reference so GC correction and autosomal library
    normalization are both engaged

The real PTEN coordinates (chr10:87.86-87.97 Mb) would need an ~88 Mb FASTA, so
the tests instead patch screen_pten_cnv.PTEN_CORE_START/END to TEST_PTEN_START/END
below and scan [SCAN_START, SCAN_END] of these short contigs.

Usage:
    python tests/make_pten_cnv_fixtures.py            # write fixtures if missing
    python tests/make_pten_cnv_fixtures.py --force    # regenerate everything
"""

import argparse
import random
from pathlib import Path

import pysam

CONTIGS = (("chr1", 200_000), ("chr2", 200_000), ("chr10", 200_000))
READ_LEN = 100
READ_STEP = 45          # ~2.2x coverage per sample
BASE_KEEP_RATE = 0.92   # per-sample downsampling, mimics library/sampling noise

SCAN_START = 40_000
SCAN_END = 160_000
BIN_SIZE = 10_000
STEP_SIZE = 5_000

# Test locus: shrunk PTEN core, and the CNV window used for DEL/DUP samples
TEST_PTEN_START = 100_000
TEST_PTEN_END = 110_000
CNV_START = 95_000
CNV_END = 115_000
DEL_FRACTION = 0.5

NORMAL_SAMPLES = tuple(f"N{i:02d}" for i in range(1, 11))
DEL_SAMPLE = "DEL1"
DUP_SAMPLE = "DUP1"

DEFAULT_OUTDIR = Path(__file__).resolve().parent / "data" / "pten_cnv"


def write_reference(path: Path) -> None:
    rng = random.Random(20261127)
    with path.open("w", encoding="utf-8") as fh:
        for name, length in CONTIGS:
            fh.write(f">{name}\n")
            seq = "".join(rng.choice("ACGT") for _ in range(length))
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    pysam.faidx(str(path))


def write_sample(outdir: Path, header: dict, sid: str, kind: str) -> None:
    contig_names = [name for name, _ in CONTIGS]
    rng = random.Random(f"{sid}:{kind}")
    out_path = outdir / f"{sid}.bam"
    read_id = 0
    with pysam.AlignmentFile(str(out_path), "wb", header=header) as out:
        for name, length in CONTIGS:
            for pos in range(200, length - 200 - READ_LEN, READ_STEP):
                if rng.random() > BASE_KEEP_RATE:
                    continue
                copies = 1
                in_window = name == "chr10" and CNV_START <= pos <= CNV_END
                if in_window and kind == "DEL":
                    if rng.random() > DEL_FRACTION:   # drop half of the reads
                        continue
                elif in_window and kind == "DUP":
                    copies = 2                        # doubles the read count
                for _ in range(copies):
                    aln = pysam.AlignedSegment()
                    aln.query_name = f"r{read_id}"
                    read_id += 1
                    aln.query_sequence = "A" * READ_LEN
                    aln.flag = 0
                    aln.reference_id = contig_names.index(name)
                    aln.reference_start = pos
                    aln.mapping_quality = 60
                    aln.cigarstring = f"{READ_LEN}M"
                    aln.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
                    out.write(aln)
    pysam.index(str(out_path))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--force", action="store_true", help="Regenerate even if the fixtures exist.")
    args = parser.parse_args()

    outdir: Path = args.outdir
    ref = outdir / "ref.fa"
    if ref.exists() and not args.force:
        print(f"[INFO] Fixtures already present in {outdir} (use --force to regenerate).")
        return

    outdir.mkdir(parents=True, exist_ok=True)
    write_reference(ref)

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": length} for name, length in CONTIGS],
    }
    samples = [(sid, "NORMAL") for sid in NORMAL_SAMPLES]
    samples += [(DEL_SAMPLE, "DEL"), (DUP_SAMPLE, "DUP")]
    for sid, kind in samples:
        write_sample(outdir, header, sid, kind)

    # Paths are relative to the fixture directory: the tests run with that as CWD
    # so the file stays machine independent.
    with (outdir / "bam_list.txt").open("w", encoding="utf-8") as fh:
        for sid, _ in samples:
            fh.write(f"{sid}\t{sid}.bam\n")

    print(f"[DONE] Wrote {len(samples)} samples to {outdir}")


if __name__ == "__main__":
    main()