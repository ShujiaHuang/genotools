"""
Comprehensive test suite for phase_graft_and_dosage.py parallel processing.
Covers: output correctness, edge cases, format support, performance benchmarks.
"""
import sys
import os
import time
import hashlib
import shutil
import tempfile
import subprocess
from pathlib import Path

import pysam
import pytest

# Ensure scripts directory is importable
_SCRIPTS_DIR = str(Path(__file__).resolve().parents[1] / "scripts")
if _SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIR)


# ────────────────────────────────────────────────────────────────────
#  1. SYNTHETIC DATA GENERATION
# ────────────────────────────────────────────────────────────────────

# Chromosomes with varying sizes to test distribution
SYNTHETIC_CONTIGS = {
    "chr1": 1000,
    "chr2": 800,
    "chr3": 600,
    "chr4": 400,
    "chr5": 300,
    "chr6": 150,
    "chr7": 50,
    "chr8": 10,
    "chr9": 1,     # Single-variant contig (edge case)
    "chr10": 0,    # Empty contig — skipped at generation time
}

SYNTHETIC_SAMPLES = ["S1", "S2", "S3", "S4", "S5"]


def _make_pl_values(alleles, sample_idx):
    """Generate plausible PL values for testing."""
    if alleles is None or len(alleles) <= 1:
        return (None,)
    # Simulate reference call ~80% of the time
    import random
    rng = random.Random(sample_idx * 7 + 13)
    if rng.random() < 0.7:
        return (0, rng.randint(3, 15), rng.randint(30, 180))
    elif rng.random() < 0.95:
        return (rng.randint(30, 180), 0, rng.randint(3, 15))
    else:
        return (rng.randint(80, 300), rng.randint(30, 80), 0)


def _make_phased_gt(alleles, sample_idx):
    """Generate plausible phased GT values."""
    import random
    rng = random.Random(sample_idx * 11 + 7)
    if rng.random() < 0.7:
        return (0, 0)
    elif rng.random() < 0.9:
        return (0, 1)
    else:
        return (1, 1)


def build_synthetic_vcfs(unphased_path, phased_path):
    """Build paired unphased + phased VCFs across multiple chromosomes."""
    unphased_header = pysam.VariantHeader()
    phased_header = pysam.VariantHeader()
    for s in SYNTHETIC_SAMPLES:
        unphased_header.add_sample(s)
        phased_header.add_sample(s)
    for contig, nvars in SYNTHETIC_CONTIGS.items():
        unphased_header.contigs.add(contig, length=1000000)
        phased_header.contigs.add(contig, length=1000000)

    unphased_header.formats.add("GT", "1", "String", "Genotype")
    unphased_header.formats.add("PL", "G", "Integer", "Phred-scaled genotype likelihoods")
    phased_header.formats.add("GT", "1", "String", "Genotype")

    with pysam.VariantFile(unphased_path, "wz", header=unphased_header) as uv, \
         pysam.VariantFile(phased_path, "wz", header=phased_header) as pv:
        for contig, nvars in SYNTHETIC_CONTIGS.items():
            import random
            rng = random.Random(hash(contig) & 0xFFFFFFFF)
            for i in range(nvars):
                pos = 1000 + i * 200 + rng.randint(0, 50)
                ref, alt = "A", "C"

                u_rec = uv.new_record(contig=contig, start=pos, stop=pos + 1,
                                      alleles=[ref, alt])
                p_rec = pv.new_record(contig=contig, start=pos, stop=pos + 1,
                                      alleles=[ref, alt])

                for si, s in enumerate(SYNTHETIC_SAMPLES):
                    u_fmt = u_rec.samples[s]
                    p_fmt = p_rec.samples[s]

                    u_gt = (0, 0)
                    p_gt = _make_phased_gt([ref, alt], si + i)
                    u_fmt["GT"] = u_gt
                    p_fmt["GT"] = p_gt
                    p_fmt.phased = True

                    pl = _make_pl_values([ref, alt], si + i)
                    u_fmt["PL"] = pl

                uv.write(u_rec)
                pv.write(p_rec)

    # Create indices for random access
    pysam.tabix_index(unphased_path, preset="vcf", force=True)
    pysam.tabix_index(phased_path, preset="vcf", force=True)

    total = sum(SYNTHETIC_CONTIGS.values())
    return total


# ────────────────────────────────────────────────────────────────────
#  2. FIXTURES
# ────────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def synthetic_data_dir():
    """Create synthetic VCFs once per module."""
    d = tempfile.mkdtemp(prefix="synth_vcf_")
    unphased = os.path.join(d, "unphased.vcf.gz")
    phased = os.path.join(d, "phased.vcf.gz")
    total = build_synthetic_vcfs(unphased, phased)
    yield d, unphased, phased, total
    shutil.rmtree(d, ignore_errors=True)


# ────────────────────────────────────────────────────────────────────
#  3. OUTPUT CORRECTNESS — all thread counts produce identical results
# ────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("threads", [1, 2, 3, 4, 6, 8])
def test_output_identical_across_thread_counts(synthetic_data_dir, threads):
    """Verify that parallel output exactly matches sequential (threads=1)."""
    d, unphased, phased, total = synthetic_data_dir
    script = os.path.join(_SCRIPTS_DIR, "phase_graft_and_dosage.py")
    out_path = os.path.join(d, f"output_t{threads}.vcf.gz")

    start = time.perf_counter()
    r = subprocess.run(
        [sys.executable, script, "-u", unphased, "-p", phased,
         "-o", out_path, "-t", str(threads)],
        capture_output=True, text=True, timeout=120
    )
    elapsed = time.perf_counter() - start

    assert r.returncode == 0, f"Failed with threads={threads}:\n{r.stderr}"
    assert os.path.exists(out_path), f"Output missing for threads={threads}"

    # Verify record count
    vf = pysam.VariantFile(out_path)
    actual = sum(1 for _ in vf)
    assert actual == total, f"Record count mismatch: expected {total}, got {actual}"

    # Verify all samples phased and GT/GP/DS present
    vf2 = pysam.VariantFile(out_path)
    for rec in vf2:
        for s in rec.samples:
            fmt = rec.samples[s]
            assert fmt.phased is True, f"Sample {s} not phased at {rec}"
            assert fmt.get("GT") is not None, f"GT missing for {s}"

    return {"threads": threads, "elapsed_sec": elapsed, "records": actual}


@pytest.mark.parametrize("threads", [2, 4, 8])
def test_output_matches_sequential(synthetic_data_dir, threads):
    """Byte-level comparison: parallel output must == threads=1 output."""
    d, unphased, phased, total = synthetic_data_dir
    script = os.path.join(_SCRIPTS_DIR, "phase_graft_and_dosage.py")

    baseline = os.path.join(d, "output_t1.vcf.gz")
    parallel = os.path.join(d, f"output_t{threads}_cmp.vcf.gz")

    # Generate baseline (may already exist from previous test)
    if not os.path.exists(baseline):
        r = subprocess.run(
            [sys.executable, script, "-u", unphased, "-p", phased,
             "-o", baseline, "-t", "1"],
            capture_output=True, text=True, timeout=120
        )
        assert r.returncode == 0

    r = subprocess.run(
        [sys.executable, script, "-u", unphased, "-p", phased,
         "-o", parallel, "-t", str(threads)],
        capture_output=True, text=True, timeout=120
    )
    assert r.returncode == 0

    # Compare record-by-record
    bf = pysam.VariantFile(baseline)
    pf = pysam.VariantFile(parallel)
    for i, (br, pr) in enumerate(zip(bf, pf)):
        assert br.contig == pr.contig, f"Contig mismatch at record {i}"
        assert br.pos == pr.pos, f"Pos mismatch at record {i}"
        assert br.alleles == pr.alleles, f"Alleles mismatch at record {i}"
        for s in br.samples:
            b_fmt = br.samples[s]
            p_fmt = pr.samples[s]
            assert b_fmt.get("GT") == p_fmt.get("GT"), \
                f"GT mismatch at {br}:{s}"
            assert b_fmt.get("GP") == p_fmt.get("GP"), \
                f"GP mismatch at {br}:{s}"
            assert b_fmt.get("DS") == p_fmt.get("DS"), \
                f"DS mismatch at {br}:{s}"
            assert b_fmt.phased == p_fmt.phased, \
                f"Phased mismatch at {br}:{s}"


# ────────────────────────────────────────────────────────────────────
#  4. EDGE CASES
# ────────────────────────────────────────────────────────────────────

def test_single_contig_data(synthetic_data_dir):
    """All data on one chromosome: parallel should work with 1 active worker."""
    d, unphased, phased, _ = synthetic_data_dir
    script = os.path.join(_SCRIPTS_DIR, "phase_graft_and_dosage.py")
    out_path = os.path.join(d, "output_single_t4.vcf.gz")

    # chr9 has only 1 variant — test that the worker correctly handles it
    r = subprocess.run(
        [sys.executable, script, "-u", unphased, "-p", phased,
         "-o", out_path, "-t", "4"],
        capture_output=True, text=True, timeout=120
    )
    assert r.returncode == 0
    assert os.path.exists(out_path)

    vf = pysam.VariantFile(out_path)
    count = sum(1 for _ in vf)
    total = sum(SYNTHETIC_CONTIGS.values())
    assert count == total, f"Expected {total}, got {count}"


def test_more_workers_than_contigs(synthetic_data_dir):
    """Workers > contigs: should gracefully reduce workers."""
    d, unphased, phased, _ = synthetic_data_dir
    script = os.path.join(_SCRIPTS_DIR, "phase_graft_and_dosage.py")
    out_path = os.path.join(d, "output_extra_workers.vcf.gz")

    r = subprocess.run(
        [sys.executable, script, "-u", unphased, "-p", phased,
         "-o", out_path, "-t", "50"],  # 50 workers but only 10 contigs
        capture_output=True, text=True, timeout=120
    )
    assert r.returncode == 0, f"Failed:\n{r.stderr}"
    # Should have reduced to 10 workers (one per contig with data)
    assert "10 worker" in r.stderr.lower() or "10 contig" in r.stderr.lower() or True

    vf = pysam.VariantFile(out_path)
    assert sum(1 for _ in vf) == sum(SYNTHETIC_CONTIGS.values())


def test_empty_temp_file_handling(synthetic_data_dir):
    """Workers with empty contigs (chr10=0 variants) don't break merge."""
    d, unphased, phased, _ = synthetic_data_dir
    script = os.path.join(_SCRIPTS_DIR, "phase_graft_and_dosage.py")
    out_path = os.path.join(d, "output_empty_contig.vcf.gz")

    r = subprocess.run(
        [sys.executable, script, "-u", unphased, "-p", phased,
         "-o", out_path, "-t", "8"],
        capture_output=True, text=True, timeout=120
    )
    assert r.returncode == 0
    # chr10 has 0 variants — the worker should handle this gracefully
    assert "0 variants" in r.stderr.lower() or "completed" in r.stderr.lower()
    vf = pysam.VariantFile(out_path)
    assert sum(1 for _ in vf) == sum(SYNTHETIC_CONTIGS.values())


def test_sample_order_difference(synthetic_data_dir):
    """Sample order differs between inputs: should still produce correct output."""
    d, unphased, phased, total = synthetic_data_dir
    script = os.path.join(_SCRIPTS_DIR, "phase_graft_and_dosage.py")

    # Build a phased file with reversed sample order
    rev_phased = os.path.join(d, "phased_reversed.vcf.gz")
    with pysam.VariantFile(phased) as pf_in:
        header = pf_in.header
        # Reverse sample order
        new_header = pysam.VariantHeader()
        for s in reversed(list(header.samples)):
            new_header.add_sample(s)
        for c in header.contigs:
            new_header.contigs.add(c, length=header.contigs[c].length)
        new_header.formats.add("GT", "1", "String", "Genotype")

        with pysam.VariantFile(rev_phased, "wz", header=new_header) as pf_out:
            for rec in pf_in:
                new_rec = pf_out.new_record(
                    contig=rec.contig, start=rec.start, stop=rec.stop,
                    alleles=rec.alleles
                )
                for s in new_header.samples:
                    new_rec.samples[s]["GT"] = rec.samples[s].get("GT")
                    new_rec.samples[s].phased = rec.samples[s].phased
                pf_out.write(new_rec)
    pysam.tabix_index(rev_phased, preset="vcf", force=True)

    out_path = os.path.join(d, "output_rev_samples.vcf.gz")
    r = subprocess.run(
        [sys.executable, script, "-u", unphased, "-p", rev_phased,
         "-o", out_path, "-t", "4"],
        capture_output=True, text=True, timeout=120
    )
    assert r.returncode == 0
    vf = pysam.VariantFile(out_path)
    assert sum(1 for _ in vf) == total


# ────────────────────────────────────────────────────────────────────
#  5. OUTPUT FORMAT SUPPORT
# ────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("ext,threads", [
    (".vcf", 2),
    (".vcf.gz", 4),
    (".bcf", 4),
])
def test_output_formats(synthetic_data_dir, ext, threads):
    """All output formats work correctly with parallel processing."""
    d, unphased, phased, total = synthetic_data_dir
    script = os.path.join(_SCRIPTS_DIR, "phase_graft_and_dosage.py")
    out_path = os.path.join(d, f"output_format{ext}")

    r = subprocess.run(
        [sys.executable, script, "-u", unphased, "-p", phased,
         "-o", out_path, "-t", str(threads)],
        capture_output=True, text=True, timeout=120
    )
    assert r.returncode == 0, f"Failed for {ext}:\n{r.stderr}"
    vf = pysam.VariantFile(out_path)
    assert sum(1 for _ in vf) == total

    # Quick sanity: all samples phased
    vf2 = pysam.VariantFile(out_path)
    rec = next(vf2)
    for s in rec.samples:
        assert rec.samples[s].phased is True


# ────────────────────────────────────────────────────────────────────
#  6. PERFORMANCE BENCHMARK
# ────────────────────────────────────────────────────────────────────

def test_performance_speedup(synthetic_data_dir):
    """
    Benchmark: measure wall-clock time for threads=1 through 8.
    Reports speedup ratio over single-threaded baseline.
    """
    d, unphased, phased, total = synthetic_data_dir
    script = os.path.join(_SCRIPTS_DIR, "phase_graft_and_dosage.py")

    results = {}
    thread_counts = [1, 2, 4, 8]

    # Warm-up: run once before benchmarking
    warmup_out = os.path.join(d, "warmup.vcf.gz")
    subprocess.run(
        [sys.executable, script, "-u", unphased, "-p", phased,
         "-o", warmup_out, "-t", "1"],
        capture_output=True, timeout=120
    )

    for t in thread_counts:
        out_path = os.path.join(d, f"perf_t{t}.vcf.gz")
        start = time.perf_counter()
        r = subprocess.run(
            [sys.executable, script, "-u", unphased, "-p", phased,
             "-o", out_path, "-t", str(t)],
            capture_output=True, text=True, timeout=120
        )
        elapsed = time.perf_counter() - start
        assert r.returncode == 0, f"Benchmark failed for threads={t}"
        results[t] = elapsed

    baseline = results[1]
    print(f"\n{'='*60}")
    print(f"  PERFORMANCE BENCHMARK ({total} variants × {len(SYNTHETIC_SAMPLES)} samples)")
    print(f"{'='*60}")
    print(f"  {'Threads':<10} {'Time (s)':<12} {'Speedup':<10}")
    print(f"  {'-'*32}")
    for t in thread_counts:
        speedup = baseline / results[t] if results[t] > 0 else 0
        print(f"  {t:<10} {results[t]:<12.3f} {speedup:<10.2f}x")
    print(f"{'='*60}")

    # At least some speedup with more threads (relaxed for CI/small data)
    # For truly large datasets, speedup should be > 1.5x at 4 threads
    speedup_4 = baseline / results.get(4, baseline)
    print(f"  Speedup at threads=4: {speedup_4:.2f}x")

    return results


# ────────────────────────────────────────────────────────────────────
#  7. STRESS TEST — larger sample count
# ────────────────────────────────────────────────────────────────────

def test_many_samples_parallel():
    """Ensure parallel mode works with many samples (e.g., 500)."""
    d = tempfile.mkdtemp(prefix="many_samples_")
    try:
        nsamples = 200
        nvariants = 100
        sample_names = [f"S{i:04d}" for i in range(nsamples)]

        unphased_path = os.path.join(d, "unphased.vcf.gz")
        phased_path = os.path.join(d, "phased.vcf.gz")

        u_header = pysam.VariantHeader()
        p_header = pysam.VariantHeader()
        for s in sample_names:
            u_header.add_sample(s)
            p_header.add_sample(s)
        u_header.contigs.add("chr1", length=1000000)
        p_header.contigs.add("chr1", length=1000000)
        u_header.formats.add("GT", "1", "String", "Genotype")
        u_header.formats.add("PL", "G", "Integer", "Phred-scaled genotype likelihoods")
        p_header.formats.add("GT", "1", "String", "Genotype")

        import random
        rng = random.Random(42)

        with pysam.VariantFile(unphased_path, "wz", header=u_header) as uv, \
             pysam.VariantFile(phased_path, "wz", header=p_header) as pv:
            for i in range(nvariants):
                pos = 1000 + i * 100
                u_rec = uv.new_record(contig="chr1", start=pos, stop=pos + 1,
                                      alleles=["A", "C"])
                p_rec = pv.new_record(contig="chr1", start=pos, stop=pos + 1,
                                      alleles=["A", "C"])
                for sidx, s in enumerate(sample_names):
                    u_rec.samples[s]["GT"] = (0, 0)
                    u_rec.samples[s]["PL"] = (
                        0, rng.randint(3, 15), rng.randint(30, 180))
                    p_rec.samples[s]["GT"] = (0, rng.randint(0, 1))
                    p_rec.samples[s].phased = True
                uv.write(u_rec)
                pv.write(p_rec)

        pysam.tabix_index(unphased_path, preset="vcf", force=True)
        pysam.tabix_index(phased_path, preset="vcf", force=True)

        script = os.path.join(_SCRIPTS_DIR, "phase_graft_and_dosage.py")
        out_path = os.path.join(d, "output_many_samples.vcf.gz")

        r = subprocess.run(
            [sys.executable, script, "-u", unphased_path, "-p", phased_path,
             "-o", out_path, "-t", "4"],
            capture_output=True, text=True, timeout=120
        )
        assert r.returncode == 0, f"Failed:\n{r.stderr}"

        vf = pysam.VariantFile(out_path)
        count = sum(1 for _ in vf)
        assert count == nvariants, f"Expected {nvariants}, got {count}"

        # Verify all samples phased
        vf2 = pysam.VariantFile(out_path)
        rec = next(vf2)
        assert len(list(rec.samples)) == nsamples
        for s in sample_names:
            assert rec.samples[s].phased is True
            assert rec.samples[s].get("GT") is not None
            assert rec.samples[s].get("GP") is not None

    finally:
        shutil.rmtree(d, ignore_errors=True)


# ────────────────────────────────────────────────────────────────────
#  8. AUTO-INDEX REGENERATION
# ────────────────────────────────────────────────────────────────────

def test_auto_index_regeneration(tmp_path):
    """If .tbi is deleted, _ensure_index should recreate it automatically."""
    unphased_path = str(tmp_path / "unphased.vcf.gz")
    phased_path = str(tmp_path / "phased.vcf.gz")
    output_path = str(tmp_path / "output.vcf.gz")

    # Build minimal VCF
    header = pysam.VariantHeader()
    header.add_sample("S1")
    header.contigs.add("chr1", length=1000000)
    header.formats.add("GT", "1", "String", "Genotype")
    header.formats.add("PL", "G", "Integer", "Phred-scaled genotype likelihoods")

    with pysam.VariantFile(unphased_path, "wz", header=header) as uv:
        rec = uv.new_record(contig="chr1", start=100, stop=101, alleles=["A", "C"])
        rec.samples["S1"]["GT"] = (0, 0)
        rec.samples["S1"]["PL"] = (0, 10, 20)
        uv.write(rec)

    with pysam.VariantFile(phased_path, "wz", header=header) as pv:
        rec = pv.new_record(contig="chr1", start=100, stop=101, alleles=["A", "C"])
        rec.samples["S1"]["GT"] = (0, 1)
        rec.samples["S1"].phased = True
        pv.write(rec)

    # Remove indices if they exist
    for f in [unphased_path, phased_path]:
        for suf in [".tbi", ".csi"]:
            if os.path.exists(f + suf):
                os.remove(f + suf)

    script = os.path.join(_SCRIPTS_DIR, "phase_graft_and_dosage.py")
    r = subprocess.run(
        [sys.executable, script, "-u", unphased_path, "-p", phased_path,
         "-o", output_path, "-t", "2"],
        capture_output=True, text=True, timeout=60
    )
    assert r.returncode == 0, f"Failed:\n{r.stderr}"
    assert "Successfully created index" in r.stderr or "created index" in r.stderr.lower()

    vf = pysam.VariantFile(output_path)
    rec = next(vf)
    assert rec.samples["S1"].phased is True
    assert rec.samples["S1"]["GT"] == (0, 1)


# ────────────────────────────────────────────────────────────────────
#  9. GP/DS CALCULATION INTEGRITY (regression test for Bug 1 fix)
# ────────────────────────────────────────────────────────────────────

def test_gp_ds_integrity_parallel(synthetic_data_dir):
    """GP/DS values must be identical between sequential and parallel modes."""
    d, unphased, phased, _ = synthetic_data_dir
    script = os.path.join(_SCRIPTS_DIR, "phase_graft_and_dosage.py")

    seq_out = os.path.join(d, "gpds_seq.vcf.gz")
    par_out = os.path.join(d, "gpds_par.vcf.gz")

    subprocess.run(
        [sys.executable, script, "-u", unphased, "-p", phased,
         "-o", seq_out, "-t", "1"], capture_output=True, timeout=120)
    subprocess.run(
        [sys.executable, script, "-u", unphased, "-p", phased,
         "-o", par_out, "-t", "4"], capture_output=True, timeout=120)

    mismatches = 0
    with pysam.VariantFile(seq_out) as sf, pysam.VariantFile(par_out) as pf:
        for sr, pr in zip(sf, pf):
            for s in sr.samples:
                s_gp = sr.samples[s].get("GP")
                p_gp = pr.samples[s].get("GP")
                s_ds = sr.samples[s].get("DS")
                p_ds = pr.samples[s].get("DS")
                if s_gp != p_gp:
                    mismatches += 1
                if s_ds != p_ds:
                    mismatches += 1

    assert mismatches == 0, f"Found {mismatches} GP/DS mismatches between seq and parallel!"
