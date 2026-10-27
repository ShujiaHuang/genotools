"""Merge

Author: Shujia Huang
Date: 2021-12-17 09:59:04
"""
import argparse
import gzip
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


def open_depth_file(path):
    """Open a 'samtools depth' output file, transparently handling gzip/bgzip and '-' (stdin)."""
    if path == "-":
        return sys.stdin
    if path.endswith((".gz", ".bgz", ".bgzip")):
        return gzip.open(path, mode="rt")
    return open(path, mode="rt")


def main():
    START_TIME = datetime.now()

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
    parser = argparse.ArgumentParser(
        description="Merge consecutive same-depth positions of 'samtools depth' output "
                    "into BED intervals. Without input files the data is read from stdin.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Example: $ samtools depth --reference ref.fa in.cram | %(prog)s > out.bed\n"
               "         $ %(prog)s out.depth out.depth.gz > out.bed"
    )
    parser.add_argument("files", nargs="*", default=["-"], type=str,
                        help="Input 'samtools depth' file(s); '.gz/.bgz/.bgzip' are "
                             "decompressed on the fly. Use '-' for stdin. When omitted, "
                             "stdin is read.")
    args = parser.parse_args()

    bed = BedGenerator()
    n = 0
    chrom, pos = "", -1

    for file_path in args.files:
        with open_depth_file(file_path) as fh:
            for line in fh:
                col = line.strip().split()
                if not col:
                    continue

                chrom, pos, depth = col[0], int(col[1]), int(col[2])
                bed.add_position(chrom, pos, depth)

                n += 1
                if n % 1000000 == 0:
                    sys.stderr.write(f"[INFO] Parsing {n} lines from '{file_path}' and hit: {chrom} - {pos}\n")

    sys.stderr.write(f"[INFO] Parsed total {n} lines and hit the end position: {chrom} - {pos}\n")
    bed.out_last_pos()

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write(f"\n** process done, {elapsed_time.seconds} seconds elapsed **\n")


if __name__ == "__main__":
    main()
