"""Random sampling samples

Author: Shujia Huang
Date: 2021-12-03 11:59:42
"""
import argparse
import sys

import numpy as np

from datetime import datetime



if __name__ == "__main__":
    START_TIME = datetime.now()

    cmdparser = argparse.ArgumentParser(description="Count variants per sample in a vcf file.")
    cmdparser.add_argument("-I", "--input", dest="input", type=str, required=True,
                           help="Input file. Required.")
    cmdparser.add_argument("-n", "--num", dest="num", type=int, required=True,
                           help="The max sample size for selecting randomly. Required.")
    args = cmdparser.parse_args()

    data = {}
    with open(args.input) as I:
        """
        #Sample	SuperReg	Province	Ethnicity	Dialect_Area	Dialect	Sex
        00113051204M47BFF2	SouthChina	Guangdong	Han	Cantonese	Cantonese	Female
        00113061190M24BFF2		Guangdong	Han	Cantonese	Cantonese	Female
        00113071013M25BFF2	SouthChina	Guangdong	Han	Cantonese	Cantonese	Female
        """
        for line in I:
            if line.startswith("#"):
                continue

            col = line.strip().split("\t")
            region = col[1]
            if region == "":
                continue

            if region not in data:
                data[region] = []

            data[region].append(col[0])

    for r, samples in data.items():
        np.random.shuffle(samples)

        for sample in samples[:args.num]:
            print ("%s\t%s" % (r, sample))

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write("\n** process done, %d seconds elapsed **\n" % elapsed_time.seconds)
