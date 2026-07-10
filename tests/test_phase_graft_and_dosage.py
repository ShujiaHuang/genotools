import importlib
import sys
from pathlib import Path

import pysam
import pytest

# Ensure the scripts directory is on sys.path so that multiprocessing
# workers can re-import the module when using the 'spawn' start method.
_SCRIPTS_DIR = str(Path(__file__).resolve().parents[1] / "scripts")
if _SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIR)

phase_graft_and_dosage = importlib.import_module("phase_graft_and_dosage")


def build_test_vcf(path: Path, samples: list[str], *, phased: bool) -> None:
    header = pysam.VariantHeader()
    for sample in samples:
        header.add_sample(sample)
    header.contigs.add("chr1")

    header.formats.add("GT", "1", "String", "Genotype")
    if not phased:
        header.formats.add("PL", "G", "Integer", "Phred-scaled genotype likelihoods")

    with pysam.VariantFile(str(path), "w", header=header) as vf:
        rec = vf.new_record(contig="chr1", start=100, stop=101, alleles=["A", "C"])
        for sample in samples:
            sample_fmt = rec.samples[sample]
            if phased:
                sample_fmt["GT"] = (0, 1) if sample == "S1" else (1, 1)
                sample_fmt.phased = True
            else:
                sample_fmt["GT"] = (0, 0) if sample == "S1" else (1, 1)
                sample_fmt["PL"] = (0, 10, 20) if sample == "S1" else (20, 10, 0)
        vf.write(rec)


def test_harmonize_cohorts_accepts_different_sample_order(tmp_path: Path) -> None:
    unphased_path = tmp_path / "unphased.vcf.gz"
    phased_path = tmp_path / "phased.vcf.gz"
    output_path = tmp_path / "output.vcf.gz"

    build_test_vcf(unphased_path, ["S1", "S2"], phased=False)
    build_test_vcf(phased_path, ["S2", "S1"], phased=True)

    phase_graft_and_dosage.harmonize_cohorts(
        str(unphased_path), str(phased_path), str(output_path), threads=2
    )

    with pysam.VariantFile(str(output_path)) as vf:
        rec = next(vf)
        assert rec.samples["S1"]["GT"] == (0, 1)
        assert rec.samples["S2"]["GT"] == (1, 1)
        assert rec.samples["S1"].phased is True
        assert rec.samples["S2"].phased is True
        assert rec.samples["S1"].get("GP") is not None
        assert rec.samples["S1"].get("DS") is not None


def test_harmonize_cohorts_skips_unsupported_single_value_pl(tmp_path: Path) -> None:
    unphased_path = tmp_path / "unphased.vcf.gz"
    phased_path = tmp_path / "phased.vcf.gz"
    output_path = tmp_path / "output.vcf.gz"

    unphased_content = """##fileformat=VCFv4.2
##source=test
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT:PL\t0/0:1
"""
    phased_content = """##fileformat=VCFv4.2
##source=test
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/0
"""

    unphased_path.write_text(unphased_content)
    phased_path.write_text(phased_content)

    phase_graft_and_dosage.harmonize_cohorts(
        str(unphased_path), str(phased_path), str(output_path), threads=1
    )

    with pysam.VariantFile(str(output_path)) as vf:
        rec = next(vf)
        assert rec.samples["S1"].get("GP") is None
        assert rec.samples["S1"].get("DS") is None


def test_harmonize_cohorts_exits_on_record_mismatch(tmp_path: Path) -> None:
    unphased_path = tmp_path / "unphased.vcf.gz"
    phased_path = tmp_path / "phased.vcf.gz"
    output_path = tmp_path / "output.vcf.gz"

    header = pysam.VariantHeader()
    header.add_sample("S1")
    header.contigs.add("chr1")
    header.formats.add("GT", "1", "String", "Genotype")
    header.formats.add("PL", "G", "Integer", "Phred-scaled genotype likelihoods")

    with pysam.VariantFile(str(unphased_path), "w", header=header) as vf:
        for pos in (100, 200):
            rec = vf.new_record(contig="chr1", start=pos, stop=pos + 1, alleles=["A", "C"])
            rec.samples["S1"]["GT"] = (0, 0)
            rec.samples["S1"]["PL"] = (0, 10, 20)
            vf.write(rec)

    with pysam.VariantFile(str(phased_path), "w", header=header) as vf:
        rec = vf.new_record(contig="chr1", start=100, stop=101, alleles=["A", "C"])
        rec.samples["S1"]["GT"] = (0, 0)
        vf.write(rec)

    with pytest.raises(SystemExit):
        phase_graft_and_dosage.harmonize_cohorts(
            str(unphased_path), str(phased_path), str(output_path), threads=1
        )
