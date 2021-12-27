"""
Author: Shujia Huang
Date: 2021-05-18 08:32:34
"""
import argparse
import sys
import gzip

from datetime import datetime
import numpy as np
import pandas as pd
import statsmodels.api as sm

START_TIME = datetime.now()


def get_beta_value(fname):
    beta = {}
    with gzip.open(fname, "rt") if fname.endswith(".gz") else open(fname, "rt") as IN:
        """
        snp	chr	pos	A	B	beta	se	pval	freq	REF	ALT
        rs340874	chr1	213985913	C	T	0.0174	0.0032	6.80E-08	0.5175	T	C
        rs1371614	chr2	26930006	T	C	0.0191	0.0045	2.36E-05	0.2513	C	T
        rs780094	chr2	27518370	C	T	0.0325	0.0032	3.30E-24	0.6044	NA	NA
        rs560887	chr2	168906638	C	T	0.0731	0.0034	4.68E-100	0.7019	T	C
        """
        for line in IN:
            if line.startswith("snp") or line.startswith("SNP"):
                continue

            col = line.strip().split()
            col[1] = col[1] if col[1].startswith("chr") else "chr" + col[1]
            pos = col[1] + ":" + col[2]

            # [Effect allele, non-Effect allele, GWAS beta value]
            beta[pos] = [col[3].upper(), col[4].upper(), float(col[5])]

    return beta


def load_fam_data(fname):
    fam = {}
    with gzip.open(fname, "rt") if fname.endswith(".gz") else open(fname, "rt") as IN:
        """
        #Family IID FID MID Sex Phenotype
        6596178	16101233BFF2	0	0	2	-9
        4459907	17200664BFF2	0	00116101038M15BFF2	1	-9
        2052894 16100773BFF2    00115111159F00BFF2      00115111159M22BFF2      2       -9
        """
        for line in IN:
            if line.startswith("#"):
                continue

            col = line.strip().split()
            sid, fid, mid = col[1], col[2], col[3]
            fam[sid] = [sid, fid, mid]

    return fam


# def distinguish_allele(maternal_GT, child_GT):
#     if (sum(maternal_GT) == 2 and sum(child_GT) == 0) or (sum(maternal_GT) == 0 and sum(child_GT) == 2):
#         # Mendelian error
#         return None, None, None
#
#     # 0,1 => 0,1
#     if (sum(maternal_GT) == 1) and (sum(child_GT) == 1):
#         # can not distinguish the maternal allele
#         return None, None, None
#
#     if sum(maternal_GT) == 0:
#         h1, h2 = 0, 0
#         h3 = 0 if sum(child_GT) == 0 else 1
#
#     elif sum(maternal_GT) == 1:  # 01 or 10
#         if sum(child_GT) == 0:  # could only be 0 or 2
#             h1, h2, h3 = 0, 1, 0
#         elif sum(child_GT) == 2:
#             h1, h2, h3 = 1, 0, 1
#         else:
#             raise ValueError("[ERROR] Child Genotype error.")
#
#     elif sum(maternal_GT) == 2:
#         h1, h2 = 1, 1
#         h3 = sum(child_GT) - h2  # child_GT) could only be [0,1]/[1,0] or [1,1]
#
#     else:
#         raise ValueError("[ERROR] Maternal genotype error!")
#
#     return h1, h2, h3


def paternal_allele_origin_by_duo(sample_gt, parent_gt, is_paternal_gt=False):
    """Determine the paternal allele index in `sample_gt` by parent-offspring duos.

    :param sample_gt: The genotype of sample
    :param parent_gt: The genotype of paternal or maternal
    :param is_paternal_gt: The `parent_gt` is from paternal. default: False
    :return:
    """
    # Default value
    is_error_genotype_match = False
    paternal_allele_origin = [0, False]  # [paternal_genotype_index, is_clear_origin]

    # Genotype should be: ["0", "0"], ["0", "1"], ["1", "0"] or ["1", "1"]
    s_gt = sample_gt.split("|")
    p_gt = parent_gt.split("|")

    s_gt_sum = sum(map(int, s_gt))  # should be: 0, 1, or 2
    p_gt_sum = sum(map(int, p_gt))  # should be: 0, 1, or 2

    if s_gt_sum == 1 and p_gt_sum == 1:
        return is_error_genotype_match, paternal_allele_origin

    if p_gt_sum == 0 and s_gt_sum == 0:  # MOM/DAD: 0|0, KID: 0|0
        paternal_allele_origin = [0, False]

    elif p_gt_sum == 0 and s_gt_sum == 1:  # MOM/DAD: 0|0, KID: 0|1 or 1|0
        if is_paternal_gt:
            paternal_allele_origin = [0, True] if s_gt[0] == "0" else [1, True]
        else:
            # If first allele is the maternal allele, the second one could only be paternal
            paternal_allele_origin = [0, True] if s_gt[0] == "1" else [1, True]

    elif p_gt_sum == 1 and (s_gt_sum == 0 or s_gt_sum == 2):  # MOM/DAD: 0|1 or 1|0, KID: 0|0, 1|1
        paternal_allele_origin = [0, False]

    elif p_gt_sum == 2 and s_gt_sum == 1:  # MOM/DAD: 1|1, KID: (0|1 or 1|0)
        if is_paternal_gt:
            paternal_allele_origin = [0, True] if s_gt[0] == "1" else [1, True]
        else:
            paternal_allele_origin = [0, True] if s_gt[0] == "0" else [1, True]

    elif p_gt_sum == 2 and s_gt_sum == 2:  # MOM/DAD: 1|1, KID: 1|1
        paternal_allele_origin = [0, False]

    else:
        is_error_genotype_match = True  # probably hit HWE

    return is_error_genotype_match, paternal_allele_origin


def paternal_allele_origin_by_trio(sample_gt, father_gt, mother_gt):
    """Determine the paternal allele index in `sample_gt` by parent-offspring duos.

    :param sample_gt: The genotype of sample
    :param father_gt: The genotype of paternal
    :param mother_gt: The genotype of maternal
    """
    # Default value
    is_error_genotype_match = False
    paternal_allele_origin = [0, False]  # [paternal_genotype_index, is_clear_origin]

    # Genotype should be: ["0", "0"], ["0", "1"], ["1", "0"] or ["1", "1"]
    s_gt = sample_gt.split("|")
    f_gt = father_gt.split("|")
    m_gt = mother_gt.split("|")

    s_gt_sum = sum(map(int, s_gt))  # should be: 0, 1, or 2
    f_gt_sum = sum(map(int, f_gt))  # should be: 0, 1, or 2
    m_gt_sum = sum(map(int, m_gt))  # should be: 0, 1, or 2

    if s_gt_sum == 1 and f_gt_sum == 1 and m_gt_sum == 1:  # 0|1, 0|1, 0|1
        return is_error_genotype_match, paternal_allele_origin

    if s_gt_sum == 0:
        if f_gt_sum != 2 or m_gt_sum != 2:
            paternal_allele_origin = [0, False]
        else:
            # DAD: 1|1 or MOM: 1|1 => impossible
            is_error_genotype_match = True

    elif s_gt_sum == 1:  # KID: 0|1 or 1|0
        if f_gt_sum == 0 and (m_gt_sum == 1 or m_gt_sum == 2):  # DAD: 0|0, MOM: 0|1 or 1|0 or 1|1
            paternal_allele_origin = [0, True] if s_gt[0] == "0" else [1, True]

        elif f_gt_sum == 1 and m_gt_sum == 0:  # DAD: 0|1 or 1|0, MOM: 0|0
            paternal_allele_origin = [0, True] if s_gt[0] == "1" else [1, True]

        elif f_gt_sum == 1 and m_gt_sum == 1:  # DAD: 0|1 or 1|0, MOM: 0|1 or 1|0
            pass  # has returned the default value

        elif f_gt_sum == 1 and m_gt_sum == 2:  # DAD: 0|1 or 1|0, MOM: 1|1
            paternal_allele_origin = [0, True] if s_gt[0] == "0" else [1, True]

        elif f_gt_sum == 2 and (m_gt_sum == 0 or m_gt_sum == 1):  # DAD: 1|1, MOM: 0|0 or 0|1 or 1|0
            paternal_allele_origin = [0, True] if s_gt[0] == "1" else [1, True]

        else:  # (0|0, 0|0), (1|1, 1|1)
            is_error_genotype_match = True

    else:  # KID: 1|1
        if f_gt_sum != 0 and m_gt_sum != 0:  # DAD: 0|1 or 1|1, MOM: 0|1 or 1|1
            paternal_allele_origin = [0, False]

        else:  # DAD == 0|0 or MOM == 0|0
            is_error_genotype_match = True

    return is_error_genotype_match, paternal_allele_origin


def offspring_genotype_origin(data, fam_idx, index2sample):
    ind_format = {name: i for i, name in enumerate(data[0][8].split(":"))}
    if "GT" not in ind_format:
        raise ValueError("[ERROR] VCF ERROR: GT is not in FORMAT.")

    paternal_allele_origin = {}  # key value is the array index of child in VCF
    for d in data:
        for c, f, m in fam_idx:  # [KID, DAD, MOM]
            father = d[f].split(":") if f is not None else None
            mother = d[m].split(":") if m is not None else None
            child = d[c].split(":")

            if (("." in child[ind_format["GT"]]) or
                    (father and "." in father[ind_format["GT"]]) or
                    (mother and "." in mother[ind_format["GT"]])):
                # contain missing genotype, ignore
                continue

            if (("/" in child[ind_format["GT"]]) or
                    (father and "/" in father[ind_format["GT"]]) or
                    (mother and "/" in mother[ind_format["GT"]])):
                raise ValueError("[ERORR] Unphased sample: %s, %s or %s. Detail: %s" % (
                    index2sample[f] if f is not None else "",
                    index2sample[m] if m is not None else "",
                    index2sample[c],
                    "\t".join(d[:7] + [d[f] if f is not None else "-",
                                       d[m] if m is not None else "-",
                                       d[c]])))

            # Genotype should be: "0|0", "0|1", "1|0" or "1|1"
            father_gt = father[ind_format["GT"]] if father else None
            mother_gt = mother[ind_format["GT"]] if mother else None
            child_gt = child[ind_format["GT"]]

            if c not in paternal_allele_origin:
                # Key value is the index of child individual in VCF line. [genotype_index, is_clear_origin]
                is_error_genotype_match, paternal_allele_origin[c] = False, [0, False]

            if paternal_allele_origin[c][1]:
                continue

            # Key value is the index of child individual in VCF line. [genotype_index, is_clear_origin]
            if father_gt is None or mother_gt is None:
                # duo
                if mother_gt:
                    is_error_genotype_match, paternal_allele_origin[c] = paternal_allele_origin_by_duo(
                        child_gt, mother_gt, is_paternal_gt=False)
                elif father_gt:
                    is_error_genotype_match, paternal_allele_origin[c] = paternal_allele_origin_by_duo(
                        child_gt, father_gt, is_paternal_gt=True)
                else:
                    # Single individual
                    is_error_genotype_match, paternal_allele_origin[c] = False, [0, False]
            else:
                # Trio
                is_error_genotype_match, paternal_allele_origin[c] = paternal_allele_origin_by_trio(
                    child_gt, father_gt, mother_gt)

            if is_error_genotype_match:
                paternal_allele_origin[c] = [0, False]
                sys.stderr.write("[ERROR] Genotype match failed but still set original and "
                                 "continue: %s \t(father: %s, %s), (mother: %s, %s) and "
                                 "(child: %s, %s).\n" % (
                                     "\t".join(d[0:5]),
                                     index2sample[f] if f is not None else "-",
                                     d[f] if f is not None else "-",
                                     index2sample[m] if m is not None else "-",
                                     d[m] if m is not None else "-",
                                     index2sample[c],
                                     d[c]))
    return paternal_allele_origin


def output_origin_phased(data, paternal_allele_origin):
    ind_format = {name: i for i, name in enumerate(data[0][8].split(":"))}
    for d in data:
        for k, c in paternal_allele_origin.items():
            ind_info = d[k].split(":")  # 0|0:0:1,0,0
            gt = ind_info[ind_format["GT"]].split("|")

            # new GT
            ind_info[ind_format["GT"]] = "|".join([gt[c[0]], gt[1 - c[0]]])
            d[k] = ":".join(ind_info)

        print("%s" % "\t".join(d))

    return


def determine_variant_parent_origin(in_vcf_fn, fam, window=10000):
    """Transform the phased Child genotype from 'Haplotype-Block-A|Haplotype-Block-B'
    to 'Paternal-Haplotype|Maternal-Haplotype'. In a case:

    KID DAD MOM
    0|1 1|1 0|0

    Transform to be:

    KID DAD MOM
    1|0 1|1 0|0

    :param in_vcf_fn:
    :param child_mother_pairs:
    :param window: The size of phasing block in Beagle.
    :return:
    """
    sample2index, index2sample = {}, {}
    fam_idx = []
    n = 0
    data_buffer = []
    prewindow = {}
    with gzip.open(in_vcf_fn, "rt") if in_vcf_fn.endswith(".gz") else open(in_vcf_fn, "rt") as IN:
        # VCF file
        for line in IN:
            if line.startswith("##"):
                print(line.strip())
                continue

            col = line.strip().split()
            if line.startswith("#CHROM"):
                print(line.strip())
                for i in range(9, len(col)):  # load sample ID and the index of sample
                    sample2index[col[i]] = i
                    index2sample[i] = col[i]

                for i in range(9, len(col)):
                    sample_id = col[i]
                    if sample_id in fam:
                        sid, fid, mid = fam[sample_id]
                        if (fid != "0" and fid not in sample2index) or (mid != "0" and mid not in sample2index):
                            raise ValueError("[ERROR] %s or %s not in VCF" % (fid, mid))

                        fam_idx.append([sample2index[sid],
                                        sample2index[fid] if fid != "0" else None,
                                        sample2index[mid] if mid != "0" else None])
                continue

            n += 1
            if n % 100000 == 0:
                elapse_time = datetime.now() - START_TIME
                sys.stderr.write("[INFO] Processing %d records done, %d seconds elapsed\n" % (n, elapse_time.seconds))

            if "," in col[4]:  # ignore multi-allelic
                continue

            chrom, pos = col[0], int(col[1])
            window_num = pos // window  # The Phased window or block
            if chrom not in prewindow:
                prewindow[chrom] = window_num

            if len(data_buffer) and (prewindow[chrom] != window_num or data_buffer[0][0] != chrom):
                paternal_allele_origin_idx = offspring_genotype_origin(data_buffer, fam_idx, index2sample)
                output_origin_phased(data_buffer, paternal_allele_origin_idx)

                data_buffer = []  # clear record
                prewindow[chrom] = window_num  # reset window

            data_buffer.append(col)  # each row is one line of VCF record

        if len(data_buffer):
            paternal_allele_origin_idx = offspring_genotype_origin(data_buffer, fam_idx, index2sample)
            output_origin_phased(data_buffer, paternal_allele_origin_idx)

    elapse_time = datetime.now() - START_TIME
    sys.stderr.write("[INFO] All %d records loaded, %d seconds elapsed.\n" % (n, elapse_time.seconds))

    return


def distinguish_origin(in_vcf_fn, fam, is_dosage=False):
    sample2index, index2sample = {}, {}
    child_mother_idx = []
    data = {}
    with gzip.open(in_vcf_fn, "rt") if in_vcf_fn.endswith(".gz") else open(in_vcf_fn, "rt") as IN:
        # VCF file
        for line in IN:
            if line.startswith("##"):
                continue

            col = line.strip().split()
            if line.startswith("#CHROM"):
                for i in range(9, len(col)):  # load sample ID and the index of sample
                    sample2index[col[i]] = i
                    index2sample[i] = col[i]

                for i in range(9, len(col)):
                    sample_id = col[i]
                    if sample_id in fam:
                        sid, _, mid = fam[sample_id]
                        if mid == "0":  # only use mother-pair samples
                            continue

                        if (mid not in sample2index) or (sid not in sample2index):
                            raise ValueError("[ERROR] %s or %s not in VCF" % (mid, sid))

                        child_mother_idx.append([sample2index[sid], sample2index[mid]])

                continue

            """
            #CHROM  POS     ID      REF     ALT
            chr7   44184122    rs730497        G       A
            """
            snp = col[2] if col[2] != "." else "-".join([col[0], col[1], col[3], col[4]])

            if "," in col[4]:  # ignore multi-allelic
                continue

            ind_format = {name: i for i, name in enumerate(col[8].split(":"))}
            if ("GT" not in ind_format) or ("GP" not in ind_format):
                raise ValueError("[ERROR]VCF ERROR: GT or GP not in FORMAT.")

            for c, m in child_mother_idx:
                # Genotype should be: [0, 0], [0, 1], [1, 0] or [1, 1]
                # `child_gt` is in "paternal_hap|maternal_hap" format after "TTC" process, so
                # `child_gt[0]` is paternal allele and `child_gt[1]` is maternal allele.
                child_gt = list(map(int, col[c].split(":")[ind_format["GT"]].split("|")))  #
                mother_gt = list(map(int, col[m].split(":")[ind_format["GT"]].split("|")))

                # Should be the probability of genotype: [0.99, 0.01, 0.00] for [Hom_Ref, Het_Var, Hom_Var]
                child_gp = list(map(float, col[c].split(":")[ind_format["GP"]].split(",")))
                mother_gp = list(map(float, col[m].split(":")[ind_format["GP"]].split(",")))

                k = index2sample[m] + "_" + index2sample[c]
                if k not in data:
                    data[k] = []

                # `h1`: maternal transmitted allele/dosage
                # `h2`: maternal non-transmitted allele/dosage
                # `h3`: paternal transmitted allele/dosage (fetal only allele)
                if is_dosage:
                    mat = mother_gp[1] + 2 * mother_gp[2]
                    fet = child_gp[1] + 2 * child_gp[2]
                    h1 = child_gt[1] * fet / max(1, sum(child_gt))  # Het: divide 1; Hom: divide 2
                    h2 = (mother_gt[0] if (mother_gt[0] != child_gt[1]) else
                          mother_gt[1]) * mat / max(1, sum(mother_gt))
                    h3 = child_gt[0] * fet / max(1, sum(child_gt))
                else:
                    mat = sum(mother_gt)
                    fet = sum(child_gt)
                    h1 = child_gt[1]  # The parent origin of variants of child has been distinguish in `TTC` process.
                    h2 = mother_gt[0] if mother_gt[0] != child_gt[1] else mother_gt[1]
                    h3 = child_gt[0]

                data[k].append([snp, mat, fet, h1, h2, h3])

    is_first_line = True
    for c, m in child_mother_idx:

        k = index2sample[m] + "_" + index2sample[c]
        if is_first_line:  # output header
            is_first_line = False
            header = ["Mother", "Child"]
            for p in data[k]:
                header.append(p[0] + "_maternal_gt")
                header.append(p[0] + "_child_gt")
                header.append(p[0] + "_h1")
                header.append(p[0] + "_h2")
                header.append(p[0] + "_h3")

            print("%s" % "\t".join(header))

        record = [index2sample[m], index2sample[c]]
        for p in data[k]:
            record.append(str(p[1]))
            record.append(str(p[2]))
            record.append(str(p[3]))
            record.append(str(p[4]))
            record.append(str(p[5]))

        print("%s" % "\t".join(record))

    return


def calculate_genotype_and_haplotype_score(in_vcf_fn, pos_beta_value, fam, is_dosage=False):
    sample2index, index2sample = {}, {}
    gs, af_beta = {}, {}
    mother_child_idx = []
    n = 0
    with gzip.open(in_vcf_fn, "rt") if in_vcf_fn.endswith(".gz") else open(in_vcf_fn, "rt") as IN:
        # VCF file
        for line in IN:
            if line.startswith("##"):
                continue

            col = line.strip().split()
            if line.startswith("#CHROM"):
                for i in range(9, len(col)):  # load sample ID and the index of sample
                    sample2index[col[i]] = i
                    index2sample[i] = col[i]

                for i in range(9, len(col)):
                    sample_id = col[i]
                    if sample_id in fam:
                        sid, fid, mid = fam[sample_id]
                        if mid == "0":  # only use mother-pair samples
                            continue

                        if (mid not in sample2index) or (sid not in sample2index):
                            raise ValueError("[ERROR] %s or %s not in VCF" % (mid, sid))

                        k = mid + "_" + sid  # mother-child
                        gs[k], af_beta[k] = [], []

                        is_trio = True if (fid != "0" and mid != "0") else False
                        mother_child_idx.append([sample2index[mid], sample2index[sid], is_trio])

                continue

            n += 1
            if n % 100000 == 0:
                elapse_time = datetime.now() - START_TIME
                sys.stderr.write("[INFO] Processing %d records done, %d seconds elapsed\n" % (n, elapse_time.seconds))

            """
            #CHROM  POS     ID      REF     ALT
            chr7   44184122    rs730497        G       A
            """
            pos = col[0] + ":" + col[1]
            ref_allele = col[3].upper()
            alt_allele = col[4].upper()

            if "," in alt_allele:  # ignore multi-allelic
                continue

            if pos not in pos_beta_value:
                continue

            # a1 is the effective allele, a2 is the non-effective allele.
            a1, a2, beta = pos_beta_value[pos]
            if (ref_allele + alt_allele != a1 + a2) and (alt_allele + ref_allele != a1 + a2):
                raise ValueError("[ERROR] Alleles not matched: "
                                 "[%s, %s] != [%s, %s]" % (ref_allele, alt_allele, a1, a2))

            if ref_allele == a1:
                beta = -1.0 * beta

            info = {c.split("=")[0]: c.split("=")[-1] for c in col[7].split(";") if "=" in c}
            af = float(info["AF"])  # ALT allele frequency

            ind_format = {name: i for i, name in enumerate(col[8].split(":"))}
            if "GT" not in ind_format:
                raise ValueError("[ERROR] 'GT' filed is required in VCF for each individual.")

            if is_dosage and (("GP" not in ind_format) and ("DS" not in ind_format)):
                raise ValueError("[ERROR] 'GP' or 'DS' field is required for dosage "
                                 "for each individual.")

            for m, c, is_duo in mother_child_idx:
                # Genotype should be: [0, 0], [0, 1], [1, 0] or [1, 1]
                # `child_gt` is in "paternal_hap|maternal_hap" format after "TTC" process, so
                # `child_gt[0]` is paternal allele and `child_gt[1]` is maternal allele.
                child_gt = list(map(int, col[c].split(":")[ind_format["GT"]].split("|")))  #
                mother_gt = list(map(int, col[m].split(":")[ind_format["GT"]].split("|")))

                if (sum(mother_gt) == 0 and sum(child_gt) == 2) or (sum(mother_gt) == 2 and sum(child_gt) == 0):
                    # Mendelian error
                    continue

                is_duo = not is_trio
                if is_duo and (sum(mother_gt) == 1) and (sum(child_gt) == 1):  # 0,1 => 0,1
                    # can not distinguish the maternal allele
                    continue

                # `h1`: maternal transmitted allele/dosage
                # `h2`: maternal non-transmitted allele/dosage
                # `h3`: paternal transmitted allele/dosage (fetal only allele)
                if is_dosage:
                    # Should be the probability of genotype: [0.99, 0.01, 0.00] for [Hom_Ref, Het_Var, Hom_Var]
                    if "GP" in ind_format:
                        mat = float(col[c].split(":")[ind_format["DS"]])
                        fet = float(col[m].split(":")[ind_format["DS"]])

                    else:  # GP in ind_format
                        child_gp = list(map(float, col[c].split(":")[ind_format["GP"]].split(",")))
                        mother_gp = list(map(float, col[m].split(":")[ind_format["GP"]].split(",")))
                        mat = mother_gp[1] + 2 * mother_gp[2]
                        fet = child_gp[1] + 2 * child_gp[2]

                    h1 = child_gt[1] * fet / max(1, sum(child_gt))  # Het: divide 1; Hom: divide 2
                    h2 = (mother_gt[0] if (mother_gt[0] != child_gt[1]) else
                          mother_gt[1]) * mat / max(1, sum(mother_gt))
                    h3 = child_gt[0] * fet / max(1, sum(child_gt))  # paternal allele
                else:
                    mat = sum(mother_gt)
                    fet = sum(child_gt)
                    h1 = child_gt[1]  # The parent origin of variants of child has been distinguish in `TTC` process.
                    h2 = mother_gt[0] if mother_gt[0] != child_gt[1] else mother_gt[1]
                    h3 = child_gt[0]

                # `s_h1`: maternal transmitted haplotype genetic score
                # `s_h2`: maternal non-transmitted haplotype genetic score
                # `s_h3`: paternal (fetal only) transmitted haplotype genetic score
                # `s_mat`: maternal genotype score
                # `s_fet`: fetal genotype score
                # `g_mat`: Maternal genetic effect
                # `g_fet`: Fetal genetic effect
                s_h1 = h1 * beta
                s_h2 = h2 * beta
                s_h3 = h3 * beta
                s_mat = mat * beta  # s_h1 + s_h2
                s_fet = fet * beta  # s_h1 + s_h3
                g_mat = (s_h1 + s_h2 - s_h3) / 2.0
                g_fet = (s_h1 + s_h3 - s_h2) / 2.0

                k = index2sample[m] + "_" + index2sample[c]
                if k not in gs:
                    gs[k] = []

                gs[k].append([s_mat, s_fet, g_mat, g_fet, s_h1, s_h2, s_h3])

    elapse_time = datetime.now() - START_TIME
    sys.stderr.write("[INFO] All %d records loaded, %d seconds elapsed.\n" % (n, elapse_time.seconds))

    # Calculate the PRS for each type of allele
    print("#Mother\tChild\tmaternal_genotype_score\tchild_genotype_score\tmaternal_genetic_effect\t"
          "fetal_genetic_effect\th1\th2\th3\tsite_number")
    for m, c, _ in mother_child_idx:
        k = index2sample[m] + "_" + index2sample[c]
        genetic_score = np.mean(gs[k], axis=0)  # Average
        print("%s\t%s\t%s\t%d" % (index2sample[m],
                                  index2sample[c],
                                  "\t".join(map(str, genetic_score)),
                                  len(gs[k])))

    return


def calculate_genotype_score(in_vcf_fn, pos_beta_value, is_dosage=False):
    """Calculate the Genetic score (or call PRS) for individuals in VCF"""
    samples = []
    sample2index, index2sample = {}, {}
    gs, af_beta = {}, {}
    n = 0
    with gzip.open(in_vcf_fn, "rt") if in_vcf_fn.endswith(".gz") else open(in_vcf_fn, "rt") as IN:
        # VCF file
        for line in IN:
            if line.startswith("##"):
                continue

            col = line.strip().split()
            if line.startswith("#CHROM"):
                for i in range(9, len(col)):  # load sample ID and the index of sample
                    sample2index[col[i]] = i
                    index2sample[i] = col[i]
                    samples.append(col[i])
                continue

            n += 1
            if n % 100000 == 0:
                elapse_time = datetime.now() - START_TIME
                sys.stderr.write("[INFO] Processing %d records done, "
                                 "%d seconds elapsed\n" % (n, elapse_time.seconds))

            """
            #CHROM  POS     ID      REF     ALT
            chr7   44184122    rs730497        G       A
            """
            pos = col[0] + ":" + col[1]
            ref_allele = col[3].upper()
            alt_allele = col[4].upper()

            if "," in alt_allele:  # ignore multi-allelic
                continue

            if pos not in pos_beta_value:
                continue

            # a1 is the effective allele, a2 is the non-effective allele.
            a1, a2, beta = pos_beta_value[pos]
            if (ref_allele + alt_allele != a1 + a2) and (alt_allele + ref_allele != a1 + a2):
                raise ValueError("[ERROR] Alleles not matched: "
                                 "[%s, %s] != [%s, %s]" % (ref_allele, alt_allele, a1, a2))

            if ref_allele == a1:
                beta = -1.0 * beta

            # info = {c.split("=")[0]: c.split("=")[-1] for c in col[7].split(";") if "=" in c}
            # af = float(info["AF"])  # ALT allele frequency

            ind_format = {name: i for i, name in enumerate(col[8].split(":"))}
            if "GT" not in ind_format:
                raise ValueError("[ERROR] 'GT' filed is required in VCF for each individual.")

            if is_dosage and (("GP" not in ind_format) and ("DS" not in ind_format)):
                raise ValueError("[ERROR] 'GP' or 'DS' field is required for dosage "
                                 "for each individual.")

            for i in range(9, len(col)):
                # Genotype should be: [0, 0], [0, 1], [1, 0] or [1, 1]
                gt_str = col[i].split(":")[ind_format["GT"]]
                if "." in gt_str:
                    # non call genotype
                    continue

                gt = list(map(int, gt_str.replace("/", "|").split("|")))
                if is_dosage:
                    # Should be the probability of genotype: [0.99, 0.01, 0.00] for [Hom_Ref, Het_Var, Hom_Var]
                    if "DS" in ind_format:
                        g = float(col[i].split(":")[ind_format["DS"]])

                    else:  # GP in ind_format
                        gp = list(map(float, col[i].split(":")[ind_format["GP"]].split(",")))
                        g = gp[1] + 2 * gp[2]

                else:
                    g = sum(gt)

                g_score = g * beta
                k = index2sample[i]
                if k not in gs:
                    gs[k] = []

                gs[k].append(g_score)  # record score for each position

    elapse_time = datetime.now() - START_TIME
    sys.stderr.write("[INFO] All %d records loaded, %d seconds elapsed.\n" % (n, elapse_time.seconds))

    # Calculate the PRS for each type of allele
    print("#SampleID\tgenotype_score\tsite_number")
    for sample in samples:
        genetic_score = np.mean(gs[sample], axis=0)  # Average
        print("%s\t%f\t%d" % (sample, genetic_score, len(gs[k])))

    return


def phenotype_concat(in_gs_fn, in_pheno_file):
    gs_data = {}
    header = []
    with open(in_gs_fn, "rt") as I:
        """
        #Sample_G1	Sample_G2	maternal_genotype_score	child_genotype_score	M1(C1)	M2	C2
        00113051204M47BFF2	13120631BFF2	-18.16394	-18.3755	-10.6309	0.000978	-10.5710
        00113061190M24BFF2	13121060BFF2	-18.17056	-18.2710	-10.6823	-0.00194	-10.5112
        """
        n = 0
        for line in I:
            n += 1
            if n == 1:
                header = line.strip().split()
                continue

            col = line.strip().split()
            gs_data[col[0]] = col

    with open(in_pheno_file, "rt") as I:
        n = 0
        for line in I:
            """
            IID	AGE	SEX	height	weight	PREBMI	Weight_gain	......
            00113051204M47BFF2	29	2	159.0	44.0	17	......
            """
            n += 1
            if n == 1:
                header += line.strip().split()[1:]
                print("%s" % "\t".join(header))
                continue

            col = line.strip().split()
            # for i in range(len(col)):
            #     if col[i] == "-9" or col[i] == "-9.0":
            #         col[i] = "NA"

            sample_id = col[0]
            if sample_id in gs_data:
                print("%s" % "\t".join(gs_data[sample_id] + col[1:]))


def regression(data, y_name, x_names, covar_names):
    """
    :param data: A dataframe
    :param y_name: column name for Y
    :param x_names: column names for X
    :param covar_names: column name for covar
    :return:
    """
    d = x_names.strip().split(",")
    c = covar_names.strip().split(",")
    x = sm.add_constant(data[d + c].copy())
    y = data[y_name]
    regression = sm.OLS(y, x)
    model = regression.fit()

    fe = pd.concat([model.params, model.bse, model.pvalues, model.tvalues], axis=1).reset_index().rename(
        columns={"index": "Feature",
                 0: "Coef",
                 1: "Stderr",
                 2: "Pvalue",
                 3: "Zscore"}).sort_values(by=["Pvalue"], inplace=False)

    return model, fe


if __name__ == "__main__":
    cmd_parser = argparse.ArgumentParser(description="Usage: ")
    commands = cmd_parser.add_subparsers(dest="command", title="Commands")
    tc_cmd = commands.add_parser("TTC", help="Transform the phased genotype of child from "
                                             "'Haplotype-Block-A|Haplotype-Block-B' to "
                                             "'Paternal-Haplotype|Maternal-Haplotype'.")
    tc_cmd.add_argument("-I", "--target", dest="target", type=str, required=True,
                        help="Input a phased VCF. Required.")
    tc_cmd.add_argument("-w", "--window", dest="window", type=int, default=10000, required=False,
                        help="The phased block size. [10000]")
    tc_cmd.add_argument("--fam", dest="fam", type=str, required=True,
                        help="Input a .fam file with mother and children.")

    ss_cmd = commands.add_parser("Split", help="Distinguish the parental origin of individual genotype "
                                               "according to the parent_origin VCF (by ``TTC``).")
    ss_cmd.add_argument("-I", "--target", dest="target", type=str, required=True,
                        help="Input a phased VCF. Required.")
    ss_cmd.add_argument("--fam", dest="fam", type=str, required=True,
                        help="Input a .fam file with mother and children.")
    ss_cmd.add_argument("--dosage", dest="dosage", action="store_true", help="Use dosage.")

    gs_cmd = commands.add_parser("GeneticScore", help="Calculate Genetic score.")
    gs_cmd.add_argument("-I", "--target", dest="target", type=str, required=True,
                        help="Input VCF. Required.")
    gs_cmd.add_argument("-b", "--base", dest="base", type=str, required=True,
                        help="A POS file with beta value for each position")
    gs_cmd.add_argument("--fam", dest="fam", type=str, required=False,
                        help="Input a .fam file [option]. If provide .fam file, this module will only "
                             "calculate the genetic score for mother-child pairs according to the "
                             "parent_origin VCF, which create by 'TTC' module.")
    gs_cmd.add_argument("--dosage", dest="dosage", action="store_true", help="Use dosage.")

    # mr_cmd = commands.add_parser("MR", help="Mendelian Randomization")
    # mr_cmd.add_argument("-I", "--input", dest="input", type=str, required=True,
    #                     help="Input data file")
    # mr_cmd.add_argument("-x", dest="x_name", type=str, required=True,
    #                     help="Load the designated phenotype(s) as x from the '--input'.")
    # mr_cmd.add_argument("--covar", dest="covar_name", type=str, required=True,
    #                     help="Only load the designated covariate(s) from the '--input'.")
    # mr_cmd.add_argument("-y", dest="y_name", type=str, required=True,
    #                     help="Load the designated data as y from the '--input'.")

    add_cmd = commands.add_parser("ADD", help="Concat genetic score data together with phenotype data.")
    add_cmd.add_argument("-g", dest="genetic_score", type=str, required=True, help="input genetic score file.")
    add_cmd.add_argument("-p", dest="phenotype", type=str, required=True, help="input phenotype data file.")

    args = cmd_parser.parse_args()
    if args.command == "TTC":
        fam_data = load_fam_data(args.fam)
        determine_variant_parent_origin(args.target, fam_data, window=args.window)

    elif args.command == "Split":
        fam_data = load_fam_data(args.fam)
        distinguish_origin(args.target, fam_data, is_dosage=args.dosage)

    elif args.command == "GeneticScore":
        if args.fam:
            fam_data = load_fam_data(args.fam)
            beta_value = get_beta_value(args.base)
            calculate_genotype_and_haplotype_score(args.target, beta_value, fam_data,
                                                   is_dosage=args.dosage)
        else:
            beta_value = get_beta_value(args.base)
            calculate_genotype_score(args.target, beta_value, is_dosage=args.dosage)

    # elif args.command == "MR":
    #     data = pd.read_table(args.input, sep="\t")
    #     regression(data, args.y_name, args.x_name, args.covar_name)

    elif args.command == "ADD":
        phenotype_concat(args.genetic_score, args.phenotype)

    else:
        cmd_parser.print_help()

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write("\n** process done, %d seconds elapsed **\n" % elapsed_time.seconds)
