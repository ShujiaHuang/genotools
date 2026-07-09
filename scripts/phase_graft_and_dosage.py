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
import math
import logging
import argparse
from typing import Tuple, Optional, Dict

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


def calculate_gp_ds(pl: Tuple[int, ...]) -> Tuple[Optional[Tuple[float, ...]], Optional[float]]:
    """
    Convert Phred-scaled likelihoods (PL) to Genotype Probabilities (GP) and Dosage (DS).
    Optimized for numerical stability in deep sequencing data.
    """
    if not pl:
        return None, None

    try:
        # Prevent math overflow/underflow by converting PL to standard likelihoods (10^(-PL/10))
        likelihoods = [math.pow(10, -val / 10.0) for val in pl]
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
            
        # Multi-allelic site fallback: calculate GP but omit DS to prevent model bias
        else:
            return tuple(round(p, 4) for p in gps), None
            
    except Exception as e:
        logger.debug(f"Numerical error parsing PL {pl}: {e}")
        return None, None


def harmonize_cohorts(unphased_in: str, phased_in: str, output_out: str, threads: int = 1) -> None:
    """
    Executes a high-throughput synchronous streaming merge between the unphased baseline
    and phased reference data. Injects exact phase states, computes accurate GP/DS values,
    and dynamically stream-compresses to the target output format.
    """
    write_mode = get_pysam_write_mode(output_out)
    logger.info(f"Input Unphased Source: {unphased_in}")
    logger.info(f"Input Phased Reference: {phased_in}")
    logger.info(f"Output Destination: {output_out} (Mode: '{write_mode}')")

    # Utilize Context Managers to guarantee file descriptors are freed, even upon pipeline failure
    # Note: pysam automatically auto-detects input formats (.vcf, .vcf.gz, .bcf) using standard mode 'r'
    with pysam.VariantFile(unphased_in, "r", threads=threads) as unphased_vcf, \
         pysam.VariantFile(phased_in, "r", threads=threads) as phased_vcf:
         
        # Harmonize metadata headers
        header = unphased_vcf.header
        header.formats.add("GP", "G", "Float", "Genotype Probabilities derived strictly from PL")
        header.formats.add("DS", "1", "Float", "Alternate Allele Dosage calculated from GP")
        
        # Verify cohort sample alignment across files
        unphased_samples = list(unphased_vcf.header.samples)
        phased_samples = list(phased_vcf.header.samples)
        if unphased_samples != phased_samples:
            logger.error(f"Sample mismatch error! Unphased cohort size: {len(unphased_samples)}, "
                         f"Phased cohort size: {len(phased_samples)}.")
            sys.exit(1)
            
        logger.info(f"Cohort validation complete ({len(unphased_samples)} samples matched). Executing stream...")
        
        processed_count = 0
        
        with pysam.VariantFile(output_out, write_mode, header=header, threads=threads) as out_vcf:
            # Sync streams variant-by-variant to maintain O(1) memory footprint
            for u_rec, p_rec in zip(unphased_vcf, phased_vcf):
                
                # Strict genomic coordinate validation to trap upstream sorting or filtering discrepancies
                if u_rec.contig != p_rec.contig or u_rec.pos != p_rec.pos:
                    logger.critical(
                        f"CRITICAL DESYNCHRONIZATION DETECTED!\n"
                        f"Unphased iterator at -> {u_rec.contig}:{u_rec.pos}\n"
                        f"Phased iterator at   -> {p_rec.contig}:{p_rec.pos}\n"
                        f"Please ensure both input files contain identical variants and sorting orders."
                    )
                    sys.exit(1)
                
                # Batch processing of samples within the current variant row
                for sample in unphased_samples:
                    u_fmt = u_rec.samples[sample]
                    p_fmt = p_rec.samples[sample]
                    
                    # 1. Physical Phase Grafting
                    u_fmt['GT'] = p_fmt['GT']
                    u_fmt.phased = p_fmt.phased
                    
                    # 2. Mathematical Genotype Dosage Reconstruction
                    pl = u_fmt.get('PL')
                    if pl is not None:
                        gp, ds = calculate_gp_ds(pl)
                        if gp is not None:
                            u_fmt['GP'] = gp
                        if ds is not None:
                            u_fmt['DS'] = ds
                            
                out_vcf.write(u_rec)
                processed_count += 1
                
                if processed_count % 50000 == 0:
                    logger.info(f"Successfully processed {processed_count} genomic variants.")
                    
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
                        help="Number of background I/O compression/decompression threads allocation.")
    
    args = parser.parse_args()
    harmonize_cohorts(args.unphased, args.phased, args.output, args.threads)


if __name__ == "__main__":
    main()



