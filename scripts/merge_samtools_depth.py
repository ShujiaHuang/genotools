"""Merge

Author: Shujia Huang
Date: 2021-12-17 09:59:04
"""
import sys
from datetime import datetime


class BedGenerator(object):
    def __init__(self):
        self.chromosome = ""
        self.start = -1
        self.last = -1
        self.depth = -1

    def add_position(self, chromo, pos, depth):

        if self.chromosome == "":
            self.chromosome = chromo
            self.start = pos
            self.last = pos
            self.depth = depth

        elif (chromo == self.chromosome) and (pos == self.last + 1) and (depth == self.depth):
            # Merge
            self.last = pos
        else:
            # Output: bed format with 1-base
            print(self.chromosome, self.start, self.last, self.depth, sep="\t")

            # reset the value
            self.chromosome = chromo
            self.start = pos
            self.last = pos
            self.depth = depth

    def out_last_pos(self):
        if self.last != -1:
            print(self.chromosome, self.start, self.last, self.depth, sep="\t")


if __name__ == "__main__":
    START_TIME = datetime.now()

    bed = BedGenerator()
    n = 0
    for line in sys.stdin:  # $ samtools depth --reference ref.fa in.cram | python merge_samtools_depth.py > out.bed
        """
        Input:
        chr1    10004   1
        chr1    10005   1
        chr1    10006   1
        chr1    10007   2
        chr1    10008   2
        chr1    10009   2
        chr1    10010   3

        Output:
        chr1    10004   10006   1
        chr1    10007   10009   2
        chr1    10010   10010   3

        """

        col = line.strip().split()
        chrom, pos, depth = col[0], int(col[1]), int(col[2])
        bed.add_position(chrom, pos, depth)

        n += 1
        if n % 1000000 == 0:
            sys.stderr.write(f"[INFO] Parsing {n} lines and hit: {chrom} - {pos}\n")

    sys.stderr.write(f"[INFO] Parsed total {n} lines and hit the end position: {chrom} - {pos}\n")
    bed.out_last_pos()

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write(f"\n** process done, {elapsed_time.seconds} seconds elapsed **\n")






