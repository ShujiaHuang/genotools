"""Number of variant per sample in vcf file

Author: Shujia Huang
Date: 2021-04-19 15:08:42
"""
import argparse
import gzip
import sys
import re

from datetime import datetime

if __name__ == "__main__":
    START_TIME = datetime.now()

    cmdparser = argparse.ArgumentParser(description="Number of variant per sample in vcf file.")
    cmdparser.add_argument("-I", "--input", dest="input", type=str, required=True,
                           help="Input VCF. Required.")
    args = cmdparser.parse_args()

    n = 0
    data = {"All": [0, 0]}
    with gzip.open(args.input, "rt") if args.input.endswith(".gz") else open(args.input, "rt") as IN:
        for line in IN:
            if line.startswith("#"):
                continue

            n += 1
            if n % 100000 == 0:
                elapsed_time = datetime.now() - START_TIME
                sys.stderr.write("[INFO] Processing %d records done, "
                                 "%d seconds elapsed\n" % (n, elapsed_time.seconds))

            col = line.strip().split()
            chrom = col[0]
            if chrom in ["chrX", "chrY", "chrM"]:
                continue

            pos, ref, alt = int(col[1]), col[3], col[4]
            if (len(ref) != 1) or (len(alt) != 1):  # Ignore non SNP site and multiple variants
                continue

            ac = int(re.search(";?AC=([^;]+)", col[7]).group(1))
            af = float(re.search(";?AF=([^;]+)", col[7]).group(1))

            af_bin = 0
            if ac == 1:
                af_bin = "AC=1"
            elif ac == 2:
                af_bin = "AC=2"
            elif af <= 0.001:
                af_bin = "0.1%"  # 0.001
            elif af <= 0.01:
                af_bin = "1%"    # 0.01
            elif af <= 0.02:    
                af_bin = "2%"    # 0.02
            elif af <= 0.05:
                af_bin = "5%"    # 0.05
            else:
                af_bin = ">5%"

            if af_bin not in data:
                data[af_bin] = [0, 0]  # [nHom, nHet]

            for i, s in enumerate(col[9:]):
                gt = s.split(":")[0]  # get sample GT
                
                if (gt == "1|1") or (gt == "1/1"):
                    data[af_bin][0] += 1
                    data["All"][0] += 1
                elif (gt == "1|0") or (gt == "0|1") or (gt == "0/1"):
                    data[af_bin][1] += 1
                    data["All"][1] += 1
                else:
                    pass

    print ("%s" % "\t".join(["#AF", "nHom", "nHet", "Het/Hom"]))
    for k in ["AC=1", "AC=2", "0.1%", "1%", "2%", "5%", ">5%", "All"]:
        v = data[k]
        print ("%s\t%d\t%d\t%.2f" % (k, v[0], v[1], v[1]/v[0] if v[0] > 0 else -9))

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write("\n** process done, %d seconds elapsed **\n" % elapsed_time.seconds)
