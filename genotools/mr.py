"""Mendelian randomization

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


def load_postion_beta_info(fname):
    data = {}
    with gzip.open(fname, "rt") if fname.endswith(".gz") else open(fname, "rt") as IN:
        """
        #CHROM	POS	BETA
        1	13273	12.8134	
        1	13649	23.4539	
        1	17375	63.2547
        """
        for line in IN:
            if line.startswith("#"):
                continue

            col = line.strip().split()
            col[0] = col[0] if col[0].startswith("chr") else "chr" + col[0]  # add 'chr'
            pos = col[0] + ":" + col[1]
            data[pos] = float(col[2])  # beta value.

    return data


def load_fam_file(fname):
    child_mother_pairs = {}
    with gzip.open(fname, "rt") if fname.endswith(".gz") else open(fname, "rt") as IN:
        """
        #Family IID FID MID Gender Phenotype
        6596178	16101233BFF2	0	00116011243M27BFF2	2	-9
        4459907	17200664BFF2	0	00116101038M15BFF2	1	-9
        1981164	00115091019M18BFF2	0	0	2	-9
        """
        for line in IN:
            if line.startswith("#"):
                continue

            col = line.strip().split()
            child, mother = col[1], col[3]
            if mother == "0":
                continue

            child_mother_pairs[child] = [mother, child]

    return child_mother_pairs


def calcute_variant_score(mother_gt_info, child_gt_info, format=None, is_dosage=True):
    """
    Format of ``mother_gt_info`` and ``child_gt_info`` should be "GT:xxx:GP"
    ``format`` is a dict, and must have 'GT' and 'GP' in it.
    """

    def haplotype_score(gt, alt_dosage_score, ref_score):
        if sum(gt) == 0:  # homo-ref
            score = [ref_score, ref_score]
        elif sum(gt) == 1:  # het
            score = [1.0 - alt_dosage_score, alt_dosage_score]
        else:  # homo variants
            score = [alt_dosage_score / 2, alt_dosage_score / 2]
        return score

    mother_data = mother_gt_info.split(":")
    if "." in mother_data[format["GT"]]:
        return None, None, None, None, None

    mother_gt = list(map(int, mother_data[format["GT"]].split("|")))  # Should be: [0,0], [0,1], [1,1]
    mother_gp = list(map(float, mother_data[format["GP"]].split(",")))  # probability for the genotype: [0.99, 0.01, 0]

    child_data = child_gt_info.split(":")
    if "." in child_data[format["GT"]]:
        return None, None, None, None, None

    child_gt = list(map(int, child_data[format["GT"]].split("|")))
    child_gp = list(map(float, child_data[format["GP"]].split(",")))

    # mother and child genotype score set to be the dosage
    m_ds = mother_gp[1] + 2 * mother_gp[2] if is_dosage else sum(mother_gt)
    c_ds = child_gp[1] + 2 * child_gp[2] if is_dosage else sum(child_gt)

    mother_hap_score = haplotype_score(mother_gt, m_ds, mother_gp[0] if is_dosage else 0)
    child_hap_score = haplotype_score(child_gt, c_ds, child_gp[0] if is_dosage else 0)

    if child_gt[0] == mother_gt[0]:
        m1 = mother_hap_score[mother_gt[0]]  # transmitted
        m2 = mother_hap_score[mother_gt[1]]  # un-transmitted allele
        c2 = child_hap_score[child_gt[1]]  # paternal allele
    elif child_gt[1] == mother_gt[1]:
        m1 = mother_hap_score[mother_gt[1]]  # transmitted
        m2 = mother_hap_score[mother_gt[0]]  # un-transmitted allele
        c2 = child_hap_score[child_gt[0]]  # paternal allele
    else:
        return None, None, None, None, None  # de novo mutation

    # mg: mother genotype-based genetic score (calculate by dosage)
    # cg: child genotype-based genetic score (calculate by dosage)
    # m1: haplotype-based genetic score of transmitted allele
    # m2: haplotype-based genetic score of un-transmitted allele
    # c2: haplotype-based genetic score of paternal transmitted allele
    return m_ds, c_ds, m1, m2, c2


def only_output(in_vcf_fn, pos_beta_value, child_mother_pairs):
    sample2index = {}
    mother_child_idx = []
    with gzip.open(in_vcf_fn, "rt") if in_vcf_fn.endswith(".gz") else open(in_vcf_fn, "rt") as IN:
        # VCF file
        for line in IN:
            if line.startswith("##"):
                continue

            col = line.strip().split()
            if line.startswith("#CHROM"):
                for i in range(9, len(col)):
                    sample2index[col[i]] = i

                header = col[:5] + ["AF", "BETA"]
                for i in range(9, len(col)):
                    sample_id = col[i]
                    if sample_id in child_mother_pairs:
                        mother, child = child_mother_pairs[sample_id]
                        if (mother not in sample2index) or (child not in sample2index):
                            raise ValueError("[ERROR] %s or %s not in VCF" % (mother, child))

                        mother_child_idx.append([sample2index[mother], sample2index[child]])
                        header.append(mother + "_" + child)
                print("%s" % "\t".join(header))
                continue

            pos = col[0] + ":" + col[1]
            if pos not in pos_beta_value:
                continue

            beta = pos_beta_value[pos]
            info = {c.split("=")[0]: c.split("=")[-1] for c in col[7].split(";") if "=" in c}
            out = col[:5] + [info["AF"], beta]  # [CHROM, POS, ID, REF, ALT, AF, Beta]

            ind_format = {name: i for i, name in enumerate(col[8].split(":"))}
            if ("GT" not in ind_format) or ("GP" not in ind_format):
                raise ValueError("[ERROR]VCF ERROR: GT or GP not in FORMAT.")

            for m, c in mother_child_idx:
                if ("|" not in col[m] and "." not in col[m]) or ("|" not in col[c] and "." not in col[c]):
                    raise ValueError("[ERORR] The VCF file must be phased. %s, %s" % (col[m], col[c]))

                mg, cg, m1, m2, c2 = calcute_variant_score(col[m], col[c], format=ind_format)
                if mg is None:
                    out.append(":".join(["."] * 5))
                else:
                    out.append(":".join(map(str, [mg, cg, m1, m2, c2])))

            print("%s" % "\t".join(map(str, out)))


def calculate_PRS(in_vcf_fn, pos_beta_value, child_mother_pairs, is_dosage=True):
    sample2index = {}
    index2sample = {}
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
                for i in range(9, len(col)):
                    sample2index[col[i]] = i
                    index2sample[i] = col[i]

                for i in range(9, len(col)):
                    sample_id = col[i]
                    if sample_id in child_mother_pairs:
                        mother, child = child_mother_pairs[sample_id]
                        if (mother not in sample2index) or (child not in sample2index):
                            raise ValueError("[ERROR] %s or %s not in VCF" % (mother, child))

                        k = mother + "_" + child
                        gs[k], af_beta[k] = [], []
                        mother_child_idx.append([sample2index[mother], sample2index[child]])
                continue

            n += 1
            if n % 100000 == 0:
                elapsed_time = datetime.now() - START_TIME
                sys.stderr.write("[INFO] Processing %d records done, "
                                 "%d seconds elapsed\n" % (n, elapsed_time.seconds))

            pos = col[0] + ":" + col[1]
            if pos not in pos_beta_value:
                continue

            beta = pos_beta_value[pos]
            info = {c.split("=")[0]: c.split("=")[-1] for c in col[7].split(";") if "=" in c}
            af = float(info["AF"])

            ind_format = {name: i for i, name in enumerate(col[8].split(":"))}
            if ("GT" not in ind_format) or ("GP" not in ind_format):
                raise ValueError("[ERROR]VCF ERROR: GT or GP not in FORMAT.")

            for m, c in mother_child_idx:
                if ("|" not in col[m] and "." not in col[m]) or ("|" not in col[c] and "." not in col[c]):
                    raise ValueError("[ERORR] The VCF file must be phased. %s, %s" % (col[m], col[c]))

                mg, cg, m1, m2, c2 = calcute_variant_score(col[m], col[c], format=ind_format,
                                                           is_dosage=is_dosage)
                k = index2sample[m] + "_" + index2sample[c]
                if mg is not None:
                    # mg: mother genotype-based genetic score (calculate by dosage)
                    # cg: child genotype-based genetic score (calculate by dosage)
                    # m1: haplotype-based genetic score of transmitted allele
                    # m2: haplotype-based genetic score of un-transmitted allele
                    # c2: haplotype-based genetic score of paternal transmitted allele
                    gs[k].append([mg, cg, m1, m2, c2])
                    af_beta[k].append([af, beta])

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write("[INFO] All %d records loaded, %d seconds elapsed\n" % (n, elapsed_time.seconds))

    # calculate the PRS for each type of allele
    print("#Sample_G1\tSample_G2\tmaternal_genotype_score\tchild_genotype_score\tM1(C1)\tM2\tC2")
    for m, c in mother_child_idx:
        k = index2sample[m] + "_" + index2sample[c]
        n = len(gs[k])
        gs[k] = np.array(gs[k])
        af_beta[k] = np.array(af_beta[k])  # [AF, Beta]

        prs = np.sum(af_beta[k][:, 1].reshape(n, 1) * gs[k], axis=0) / n
        print("%s\t%s\t%s" % (index2sample[m], index2sample[c], "\t".join(map(str, prs))))

    return


def add_phenotype(in_prs_fn, in_pheno_file):
    prs_data = {}
    header = []
    with open(in_prs_fn, "rt") as I:
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
            prs_data[col[0]] = col

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
            sample_id = col[0]
            if sample_id in prs_data:
                print("%s" % "\t".join(prs_data[sample_id] + col[1:]))


def mendelian_randomization(data, y_name, x_names, covar_names):
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

    print(model.summary())
    print(fe)
    return model, fe


if __name__ == "__main__":
    cmdparser = argparse.ArgumentParser(description="Usage: ")
    commands = cmdparser.add_subparsers(dest="command", title="Commands")
    prs_cmd = commands.add_parser("PRS", help="Calculate PRS")
    prs_cmd.add_argument("-I", "--input", dest="input", type=str, required=True,
                         help="Input VCF. Required.")
    prs_cmd.add_argument("--famfile", dest="famfile", type=str, required=True,
                         help="A .fam file")
    prs_cmd.add_argument("--betavaluefile", dest="betavaluefile", type=str, required=True,
                         help="A POS file with beta value for each position")

    mr_cmd = commands.add_parser("MR", help="Mendelian Randomization")
    mr_cmd.add_argument("-I", "--input", dest="input", type=str, required=True,
                        help="Input data file")
    mr_cmd.add_argument("--y-name", dest="y_name", type=str, required=True,
                        help="Load the designated data as y from the '--input'.")
    mr_cmd.add_argument("--x-name", dest="x_name", type=str, required=True,
                        help="Load the designated phenotype(s) as x from the '--input'.")
    mr_cmd.add_argument("--covar-name", dest="covar_name", type=str, required=True,
                        help="Only load the designated covariate(s) from the '--input'.")

    add_cmd = commands.add_parser("ADD", help="Concat PRS data together with phenotype data.")
    add_cmd.add_argument("--PRS", dest="prs", type=str, required=True, help="PRS result.")
    add_cmd.add_argument("--pheno", dest="pheno", type=str, required=True, help="phenotype data.")

    if len(sys.argv) < 2:
        cmdparser.print_help()
        sys.exit(1)

    args = cmdparser.parse_args()

    if args.command == "PRS":
        child_mother_pairs = load_fam_file(args.famfile)
        pos_beta_value = load_postion_beta_info(args.betavaluefile)
        sys.stderr.write("[INFO] Load beta value done.\n")
        calculate_PRS(args.input, pos_beta_value, child_mother_pairs, is_dosage=False)
    elif args.command == "MR":
        data = pd.read_table(args.input, sep="\t")
        mendelian_randomization(data, args.y_name, args.x_name, args.covar_name)
    elif args.command == "ADD":
        add_phenotype(args.prs, args.pheno)

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write("\n** process done, %d seconds elapsed **\n" % elapsed_time.seconds)
