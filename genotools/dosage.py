"""Calculating genotype dosage.

Author: Shujia Huang
Date: 2020-05-09 15:37:02
"""
import sys
import gzip
import argparse
import numpy as np
from datetime import datetime


def calculate_dosage_PL(pl, is_normalized=False):
    w = np.array([0, 1, 2])  # Weight of Genotype for RR, R/A and A/A in dosage calculation
    if len(pl) < 3:
        # PARs region of X Chromosome and there just 2 element in PL field.
        w = np.array([0, 2])

    pr = 10 ** (-0.1 * np.array(pl, dtype=float))  # Transformation Phred-scale to be Probability
    if is_normalized:
        pr = recover_PL(pr)

    return round(sum(w * pr), 3)


def recover_PL(pr):
    # Recover the normalized PL
    return pr / pr.sum()


if __name__ == "__main__":
    START_TIME = datetime.now()

    cmdparser = argparse.ArgumentParser(description="Print genotype dosage")
    cmdparser.add_argument("-I", "--input", dest="input", type=str, required=True,
                           help="Input VCF. Required.")
    cmdparser.add_argument("--vcf", dest="vcf", action="store_true", help="Output in VCF format.")
    cmdparser.add_argument("-N", "--is-gatk-normalized-PL", dest="normalizePL", action="store_true",
                           help="PL value in the VCF is been normalized in the GATK process, which "
                                "you can find more detail here: https://gatk.broadinstitute.org/hc/en-us/articles/360035890451-Calculation-of-PL-and-GQ-by-HaplotypeCaller-and-GenotypeGVCFs")
    args = cmdparser.parse_args()

    DS = "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"estimated ALT dose [P(RA) + 2*P(AA)]\">"
    with gzip.open(args.input, "rt") if args.input.endswith(".gz") else open(args.input, "rt") as IN:
        for line in IN:
            if line.startswith("##"):
                if args.vcf:
                    # Output VCF header
                    if line.startswith("##FORMAT=<ID=DP"):
                        print(line.strip())
                        print(DS)
                    elif line.startswith("##FORMAT=<ID=GT"):
                        print(DS)
                        print(line.strip())
                    else:
                        print(line.strip())
                continue

            col = line.strip().split()
            if line.startswith("#CHROM"):
                if args.vcf:
                    print(line.strip())
                else:
                    print("\t".join(col[:5] + col[9:]))
                continue

            # GT:AD:DP:GQ:PGT:PID:PL
            format_idx = {s: i for i, s in enumerate(col[8].split(":"))}
            allelic_num = len(col[4].split(","))  # could be multi-allelic. e.g: A,G

            if allelic_num > 2:
                # ignore any positions which it's trio-allelic
                sys.stderr.write("[WARNING] Trio-allelic: %s" % line.strip())
                continue

            if args.vcf:
                col[8] += ":DS"

            # set output field
            fields = col[:9] if args.vcf else col[:5]

            # calculate the samples' dosage one by one
            for sample in col[9:]:
                ds = []
                pl = sample.split(":")[format_idx["PL"]]
                if pl != ".":
                    pls = np.array(pl.split(","), dtype=float)
                    is_PARs_region_in_X_chromosome = True if len(pls) % 3 else False
                    for i in range(allelic_num):
                        ds.append(
                            calculate_dosage_PL(pls[i:i + 3] if not is_PARs_region_in_X_chromosome else pls[i:i + 2],
                                                is_normalized=args.normalizePL))
                else:
                    ds = [-1]

                fields.append(sample + ":%s" % ",".join(map(str, ds)) if args.vcf else ",".join(map(str, ds)))

            print("\t".join(fields))

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write("\n** process done, %d seconds elapsed **\n" % elapsed_time.seconds)
