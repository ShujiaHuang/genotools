"""Based on the ratio of X chromosome-derived shotgun sequencing data to the autosomal coverage to establish the probability of an XX or XY karyotype for samples.

Author: Shujia Huang
Date: 2022-08-02
"""
import argparse
import sys
import gzip

from datetime import datetime
import numpy as np

START_TIME = datetime.now()


def infer_sex(ci1, ci2):

    sex = "-" 
    if (ci1 > 0.8) {  # The sample should be assigned as Female
        sex = "Female"
    } elif (ci2 < 0.6) {  # The sample should be assigned as Male
        sex = "Male"
    } elif (ci1 > 0.6 and ci2 > 0.8) {
        # The sample is consistent with XX but not XY
        sex = "XX"
    } elif (ci1 < 0.6 and ci2 < 0.8) {
        # The sample is consistent with XY but not XX
        sex = "XY"
    } else { # The sample could not be assigned
        sex = "-"
    }

    return sex


def main(fname, samplename):
    autosome_id = set(["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                        "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
                        "chr18","chr19","chr20","chr21","chr22"])

    chromosome_id = set(["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                         "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
                         "chr18","chr19","chr20","chr21","chr22","chrX","chrY"])


    chromo_len = {}
    chromo_map = {}
    total_len, total_map = 0, 0
    with gzip.open(fname, "rt") if fname.endswith(".gz") else open(fname, "rt") as IN:
        """ Zero-based in bed format by bedtools.

        chr1    0       9999    0
        chr1    9999    10000   1
        chr1    10000   10001   3
        chr1    10001   10002   6
        chr1    10002   10003   11
        chr1    10003   10004   18
        chr1    10004   10005   21
        """
        for line in IN:
            if line.startswith("#"):
                continue

            col = line.strip().split()
            chr_name, start, end, depth = col[0], int(col[1]), int(col[2]), int(col[3])

            if chr_name not in chromosome_id:
                continue

            total_len += (end - start)
            total_map += ((end - start) * depth)  # total mapped bases

            if chr_name not in chromo_len:
                chromo_len[chr_name] = 0
                chromo_map[chr_name] = 0
                sys.stderr.write("*** processing %s ***\n" % chr_name)

            chromo_len[chr_name] += (end - start)
            chromo_map[chr_name] += ((end - start) * depth)


    rt = {k:(chromo_map[k]/total_map) / (chromo_len[k]/total_len) for k in chromosome_id}
    xa = [rt["chrX"]/rt[k] for k in autosome_id]

    rx = np.mean(xa)
    ci = 1.96 * np.std(xa) / np.sqrt(len(xa))
    ci1, ci2 = rx - ci, rx + ci

    sex = infer_sex(ci1, ci2)
    print("#Sample_id\tRx\t95%-CI-Lower\t95%-CI-Upper\tInferred-Sex")
    print("%s" % "\t".join(map(str, [samplename, rx, ci1, ci2, sex])))



if __name__ == '__main__':

    cmdparser = argparse.ArgumentParser(description="Based on the ratio of X chromosome-derived "
        "shotgun sequencing data to the autosomal coverage to establish the probability of an XX "
        "or XY karyotype for samples.")
    cmdparser.add_argument("-b", "--bed", dest="bed", type=str, required=True,
                           help="Input the coverage file in .bed format. Required.")
    cmdparser.add_argument("-s", "--samplename", dest="samplename", type=str, required=True,
                           help="Sample name. Required.")
    args = cmdparser.parse_args()

    main(args.bed, args.samplename)

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write("\n** process done, %d seconds elapsed **\n" % elapsed_time.seconds)



