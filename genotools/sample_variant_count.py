"""Count variants per-sample in a vcf file

Author: Shujia Huang
Date: 2020-10-09 15:24:37
"""
import argparse
import gzip
import sys

from datetime import datetime

if __name__ == "__main__":
    START_TIME = datetime.now()

    cmdparser = argparse.ArgumentParser(description="Count variants per sample in a vcf file.")
    cmdparser.add_argument("-I", "--input", dest="input", type=str, required=True,
                           help="Input VCF. Required.")
    args = cmdparser.parse_args()

    sample_varaints_count = {}
    samples = []
    n = 0
    with gzip.open(args.input, "rt") if args.input.endswith(".gz") else open(args.input, "rt") as IN:
        for line in IN:
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                col = line.strip().split()
                for s in col[9:]:
                    samples.append(s)
                    sample_varaints_count[s] = [0, 0, 0, 0]
                continue

            n += 1
            if n % 100000 == 0:
                elapsed_time = datetime.now() - START_TIME
                sys.stderr.write("[INFO] Processing %d records done, " 
                                 "%d seconds elapsed\n" % (n, elapsed_time.seconds))

            col = line.strip().split()
            ref = col[3]
            alt = col[4].split(",")  # may be multiple variants
            variant_type = {str(i+1): len(a)-len(ref) for i, a in enumerate(alt)}
            for i, s in enumerate(col[9:]):
                gt = s.split(":")[0]  # get sample GT

                if gt != "0/0" and gt != "0|0" and gt != "0" and ("." not in gt):
                    sample_varaints_count[samples[i]][0] += 1

                    vt = variant_type[gt[-1]] if gt[-1] in variant_type else variant_type[gt[0]]
                    if vt == 0:  # SNP
                        sample_varaints_count[samples[i]][1] += 1
                    elif vt > 0:  # Insertion
                        sample_varaints_count[samples[i]][2] += 1
                    else:  # Deletion
                        sample_varaints_count[samples[i]][3] += 1

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write("[INFO] All %d records loaded, %d seconds elapsed\n" % (n, elapsed_time.seconds))

    print("#SAMPLE_ID\tVariantCount\tSNP\tInsertion\tDeletion")
    for s in samples:
        print("%s\t%s" % (s, "\t".join(map(str, sample_varaints_count[s]))))

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write("\n** process done, %d seconds elapsed **\n" % elapsed_time.seconds)
