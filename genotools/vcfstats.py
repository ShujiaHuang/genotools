"""Distribution of the variants's length

Author: Shujia Huang
Date: 2021-04-19 15:08:42
"""
import argparse
import gzip
import sys
import re

from datetime import datetime


def is_overlap(chrom, start_pos, end_pos, bed_region, index):
    
    if chrom not in bed_region:
        return False

    is_overlap = False
    start_idx = index[chrom]
    for i in range(start_idx, len(bed_region[chrom])):

        if start_pos > bed_region[chrom][i][1]: continue
        if end_pos   < bed_region[chrom][i][0]: break

        index[chrom] = i
        is_overlap = True

        break

    return is_overlap


if __name__ == "__main__":
    START_TIME = datetime.now()

    cmdparser = argparse.ArgumentParser(description="Count variants per sample in a vcf file.")
    cmdparser.add_argument("-I", "--input", dest="input", type=str, required=True,
                           help="Input VCF. Required.")
    cmdparser.add_argument("-b", "--bed", dest="bed", type=str, required=True,
                           help="Input CDS region. bed format. Required.")
    args = cmdparser.parse_args()

    bed_region = {}
    index = {}
    with open(args.bed) as I:
        for line in I:
            col = line.strip().split()
            chr_id, start, end = col[0], int(col[1]), int(col[2])
            bed_region.setdefault(chr_id, []).append([start, end])  # sorted
            index[chr_id] = 0

    n = 0
    print ("%s" % "\t".join(["CHROM", "POS", "REF", "ALT", "AC", "Raw_AF", "AF", "Variant_length", "Variant_type", "Known_or_novel", "Genome_region"]))
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
            chrom = col[0]
            pos = int(col[1])
            variants_kn = "Known" if col[2] != "." else "Novel"
            ref = col[3]
            alt = col[4].split(",")  # may be multiple variants
            ac =   int(re.search(";?AC=([^;]+)", col[7]).group(1))
            af = float(re.search(";?AF=([^;]+)", col[7]).group(1))
            variant_length = [len(a)-len(ref) for i, a in enumerate(alt)]
            variant_type = ["SNP" if s == 0 else ("Deletion" if s<0 else "Insertion") for s in variant_length]

            if len(variant_type) > 1:  # ignore multiple variants
                continue

            is_in_coding_region = is_overlap(chrom, 
                                             pos, 
                                             pos+abs(variant_length[0]) if variant_type[0] == "Deletion" else pos, 
                                             bed_region,
                                             index)
            g_region = "CDS" if is_in_coding_region else "Other"

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

            print ("\t".join(map(str, [chrom, pos, ref] + alt + [ac, af, af_bin] + variant_length + variant_type + [variants_kn, g_region])))
            

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write("\n** process done, %d seconds elapsed **\n" % elapsed_time.seconds)
