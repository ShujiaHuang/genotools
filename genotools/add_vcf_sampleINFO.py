"""

Author: Shujia Huang
Date: 2020-10-24 19:10:31
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
    cmdparser.add_argument("-S", "--sample", dest="sample", type=str, required=True,
                           help="SAMPLE. Required.")
    cmdparser.add_argument("-H", "--header", dest="header", type=str, required=True,
                           help="VCF Header info. Required.")
    args = cmdparser.parse_args()

    with open(args.header) as I:
        for line in I:
            print("%s" % line.strip())


    header_9col = {}
    with gzip.open(args.input, "rt") if args.input.endswith(".gz") else open(args.input, "rt") as IN:
        for line in IN:
            if line.startswith("#"):
                continue
            col = line.strip().split()
            k = col[0] + ":" + col[1] + ":" + col[2]
            header_9col[k] = col

    n = 0
    with gzip.open(args.sample, "rt") if args.sample.endswith(".gz") else open(args.sample, "rt") as IN:
        for line in IN:
            if line.startswith("#"):
                continue

            n += 1
            if n % 100000 == 0:
                elapsed_time = datetime.now() - START_TIME
                sys.stderr.write("[INFO] Processing %d records done, %d seconds elapsed\n" % (n, elapsed_time.seconds))

            col = line.strip().split()
            k = col[0] + ":" + col[1] + ":" + col[2]

            if k in header_9col:
                print("%s" % "\t".join(header_9col[k] + col[8:]))

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write("\n** process done, %d seconds elapsed **\n" % elapsed_time.seconds)