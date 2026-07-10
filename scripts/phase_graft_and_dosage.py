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

import os
import sys
import math
import logging
import tempfile
import argparse
import multiprocessing
from typing import Tuple, Optional

import pysam

try:
    import numpy as np
except ImportError:  # pragma: no cover - optional dependency
    np = None


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


def calculate_gp_ds(pl: Tuple[int, ...]) -> Tuple[Optional[Tuple[float, ...]], Optional[float]]:
    """
    Convert Phred-scaled likelihoods (PL) to Genotype Probabilities (GP) and Dosage (DS).
    """
    if not pl:
        return None, None

    if len(pl) not in (2, 3):
        logger.debug("Skipping GP/DS generation for unsupported PL length %d: %r", len(pl), pl)
        return None, None

    try:
        if np is not None:
            arr = np.asarray(pl, dtype=np.float64)
            likelihoods = np.power(10.0, -arr / 10.0)
            sum_l = float(likelihoods.sum())
            if sum_l == 0:
                return None, None
            gps = likelihoods / sum_l
            if len(pl) == 3:
                ds = float(gps[1] + 2.0 * gps[2])
                return tuple(round(float(p), 4) for p in gps), round(ds, 4)
            elif len(pl) == 2:
                return tuple(round(float(p), 4) for p in gps), round(float(gps[1]), 4)
            return tuple(round(float(p), 4) for p in gps), None

        likelihoods = [math.pow(10, -val / 10.0) for val in pl]
        sum_l = sum(likelihoods)

        if sum_l == 0:
            return None, None

        gps = tuple(l / sum_l for l in likelihoods)

        if len(pl) == 3:
            ds = gps[1] + 2.0 * gps[2]
            return tuple(round(p, 4) for p in gps), round(ds, 4)
        elif len(pl) == 2:
            return tuple(round(p, 4) for p in gps), round(gps[1], 4)
        return tuple(round(p, 4) for p in gps), None

    except Exception as e:
        logger.debug(f"Numerical error parsing PL {pl}: {e}")
        return None, None


def _process_contig_chunk(args: Tuple[str, str, str, str, int, Tuple[str, ...], Tuple[str, ...], str]) -> Tuple[str, str, int]:
    contig, unphased_in, phased_in, temp_output, threads, unphased_samples, phased_samples, write_mode = args
    with pysam.VariantFile(unphased_in, "r", threads=threads) as unphased_vcf, \
         pysam.VariantFile(phased_in, "r", threads=threads) as phased_vcf:

        header = unphased_vcf.header.copy()
        if "GT" not in header.formats:
            header.formats.add("GT", "1", "String", "Genotype")
        if "GP" not in header.formats:
            header.formats.add("GP", "G", "Float", "Genotype Probabilities derived strictly from PL")
        if "DS" not in header.formats:
            header.formats.add("DS", "1", "Float", "Alternate Allele Dosage calculated from GP")

        for target_header in (unphased_vcf.header, phased_vcf.header):
            if "GT" not in target_header.formats:
                target_header.formats.add("GT", "1", "String", "Genotype")
            if "GP" not in target_header.formats:
                target_header.formats.add("GP", "G", "Float", "Genotype Probabilities derived strictly from PL")
            if "DS" not in target_header.formats:
                target_header.formats.add("DS", "1", "Float", "Alternate Allele Dosage calculated from GP")

        with pysam.VariantFile(temp_output, write_mode, header=header, threads=threads) as out_vcf:
            u_iter = unphased_vcf.fetch(contig=contig)
            p_iter = phased_vcf.fetch(contig=contig)
            u_rec = next(u_iter, None)
            p_rec = next(p_iter, None)
            processed_count = 0

            while u_rec is not None and p_rec is not None:
                if u_rec.pos != p_rec.pos:
                    logger.critical(
                        f"CRITICAL DESYNCHRONIZATION DETECTED FOR {contig}!\n"
                        f"Unphased iterator at -> {u_rec.contig}:{u_rec.pos}\n"
                        f"Phased iterator at   -> {p_rec.contig}:{p_rec.pos}\n"
                        f"Please ensure both input files contain identical variants and sorting orders."
                    )
                    raise ValueError("Record desynchronization detected")

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
                        try:
                            gp, ds = calculate_gp_ds(pl)
                            if gp is not None:
                                u_fmt["GP"] = gp
                            if ds is not None:
                                u_fmt["DS"] = ds
                        except TypeError:
                            logger.debug("Skipping GP/DS assignment for unsupported PL tuple %r", pl)

                out_vcf.write(u_rec)
                processed_count += 1

                u_rec = next(u_iter, None)
                p_rec = next(p_iter, None)

            if u_rec is not None or p_rec is not None:
                raise ValueError(f"Unequal record counts detected for contig {contig}")

    return contig, temp_output, processed_count


def harmonize_cohorts(unphased_in: str, phased_in: str, output_out: str, threads: int = 1) -> None:
    """
    Executes a high-throughput harmonization between the unphased baseline and phased reference data.
    Processes chromosomes independently when possible, which allows the workload to be distributed
    across multiple CPU workers while keeping the output order stable.
    """
    write_mode = get_pysam_write_mode(output_out)
    logger.info(f"Input Unphased Source: {unphased_in}")
    logger.info(f"Input Phased Reference: {phased_in}")
    logger.info(f"Output Destination: {output_out} (Mode: '{write_mode}')")

    with pysam.VariantFile(unphased_in, "r", threads=threads) as unphased_vcf, \
         pysam.VariantFile(phased_in, "r", threads=threads) as phased_vcf:

        header = unphased_vcf.header.copy()
        if "GT" not in header.formats:
            header.formats.add("GT", "1", "String", "Genotype")
        if "GP" not in header.formats:
            header.formats.add("GP", "G", "Float", "Genotype Probabilities derived strictly from PL")
        if "DS" not in header.formats:
            header.formats.add("DS", "1", "Float", "Alternate Allele Dosage calculated from GP")

        for target_header in (unphased_vcf.header, phased_vcf.header):
            if "GT" not in target_header.formats:
                target_header.formats.add("GT", "1", "String", "Genotype")
            if "GP" not in target_header.formats:
                target_header.formats.add("GP", "G", "Float", "Genotype Probabilities derived strictly from PL")
            if "DS" not in target_header.formats:
                target_header.formats.add("DS", "1", "Float", "Alternate Allele Dosage calculated from GP")

        unphased_samples = list(unphased_vcf.header.samples)
        phased_samples = list(phased_vcf.header.samples)
        if set(unphased_samples) != set(phased_samples):
            logger.error(
                f"Sample mismatch error! Unphased samples: {unphased_samples[:5]}..., "
                f"Phased samples: {phased_samples[:5]}..."
            )
            sys.exit(1)

        contigs = [contig for contig in unphased_vcf.header.contigs if contig in phased_vcf.header.contigs]
        if not contigs:
            logger.warning("No shared contigs found between the input files.")
            contigs = []

        logger.info(f"Cohort validation complete ({len(unphased_samples)} samples matched). Executing stream...")

        if threads <= 1 or len(contigs) <= 1:
            processed_count = 0
            with pysam.VariantFile(output_out, write_mode, header=header, threads=threads) as out_vcf:
                while True:
                    u_rec = next(unphased_vcf, None)
                    p_rec = next(phased_vcf, None)

                    if u_rec is None and p_rec is None:
                        break
                    if u_rec is None or p_rec is None:
                        logger.critical(
                            "CRITICAL DESYNCHRONIZATION DETECTED!\n"
                            f"Unphased iterator exhausted: {u_rec is None}\n"
                            f"Phased iterator exhausted: {p_rec is None}"
                        )
                        sys.exit(1)

                    if u_rec.contig != p_rec.contig or u_rec.pos != p_rec.pos:
                        logger.critical(
                            f"CRITICAL DESYNCHRONIZATION DETECTED!\n"
                            f"Unphased iterator at -> {u_rec.contig}:{u_rec.pos}\n"
                            f"Phased iterator at   -> {p_rec.contig}:{p_rec.pos}\n"
                            f"Please ensure both input files contain identical variants and sorting orders."
                        )
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
                            try:
                                gp, ds = calculate_gp_ds(pl)
                                if gp is not None:
                                    u_fmt["GP"] = gp
                                if ds is not None:
                                    u_fmt["DS"] = ds
                            except TypeError:
                                logger.debug("Skipping GP/DS assignment for unsupported PL tuple %r", pl)

                    out_vcf.write(u_rec)
                    processed_count += 1

                logger.info(f"Execution finished. Total processed variants: {processed_count}.")
            return

        temp_dir = tempfile.mkdtemp(prefix="phase_graft_", dir="/tmp")
        try:
            args = []
            for contig in contigs:
                temp_path = os.path.join(temp_dir, f"{contig}.tmp")
                args.append((contig, unphased_in, phased_in, temp_path, max(1, threads // len(contigs)), unphased_samples, phased_samples, write_mode))

            with multiprocessing.get_context("spawn").Pool(processes=min(threads, len(contigs))) as pool:
                results = pool.map(_process_contig_chunk, args)

            processed_count = 0
            with pysam.VariantFile(output_out, write_mode, header=header, threads=threads) as out_vcf:
                for contig, temp_path, count in sorted(results, key=lambda item: contigs.index(item[0])):
                    with pysam.VariantFile(temp_path, "r", threads=threads) as chunk_vcf:
                        for rec in chunk_vcf:
                            out_vcf.write(rec)
                    processed_count += count
                    os.remove(temp_path)

            logger.info(f"Execution finished. Total processed variants: {processed_count}.")
        finally:
            for filename in os.listdir(temp_dir):
                os.remove(os.path.join(temp_dir, filename))
            os.rmdir(temp_dir)


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
                        help="Number of background I/O compression/decompression threads allocation.")
    
    args = parser.parse_args()
    harmonize_cohorts(args.unphased, args.phased, args.output, args.threads)


if __name__ == "__main__":
    main()



