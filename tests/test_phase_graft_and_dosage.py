import importlib.util
from pathlib import Path

import pysam


MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "phase_graft_and_dosage.py"
SPEC = importlib.util.spec_from_file_location("phase_graft_and_dosage", MODULE_PATH)
phase_graft_and_dosage = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(phase_graft_and_dosage)


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
        str(unphased_path), str(phased_path), str(output_path), threads=1
    )

    with pysam.VariantFile(str(output_path)) as vf:
        rec = next(vf)
        assert rec.samples["S1"]["GT"] == (0, 1)
        assert rec.samples["S2"]["GT"] == (1, 1)
        assert rec.samples["S1"].phased is True
        assert rec.samples["S2"].phased is True
        assert rec.samples["S1"].get("GP") is not None
        assert rec.samples["S1"].get("DS") is not None
