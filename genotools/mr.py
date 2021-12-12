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


def get_child_mother_duos(fname):
    child_mother_pairs = {}
    with gzip.open(fname, "rt") if fname.endswith(".gz") else open(fname, "rt") as IN:
        """
        #Family IID FID MID Gender Phenotype
        6596178	16101233BFF2	0	00116011243M27BFF2	2	-9
        4459907	17200664BFF2	0	00116101038M15BFF2	1	-9
        """
        for line in IN:
            if line.startswith("#"):
                continue

            col = line.strip().split()
            child, mother = col[1], col[3]
            if mother == "0":  # Only get child-mother pairs.
                continue

            child_mother_pairs[child] = [mother, child]

    return child_mother_pairs


def calcute_variant_score(mother_gt_info, child_gt_info, format=None, is_dosage=True):
    """
    Format of ``mother_gt_info`` and ``child_gt_info`` should be "GT:xxx:GP"
    ``format`` is a dict, and must have 'GT' and 'GP' in it.
    """

    def haplotype_score(gt, alt_score):
        if gt == 0:  # homo-ref
            score = [0.0, 0.0]
        elif gt == 1:  # het
            score = [0.0, alt_score / 2]
        else:  # homo variants
            score = [alt_score / 2, alt_score / 2]
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

    mother_hap_score = haplotype_score(sum(mother_gt), m_ds)
    child_hap_score = haplotype_score(sum(child_gt), c_ds)

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


def distinguish_allele(maternal_GT, child_GT):
    if (sum(maternal_GT) == 2 and sum(child_GT) == 0) or (sum(maternal_GT) == 0 and sum(child_GT) == 2):
        # Mendelian error
        return None, None, None

    # 0,1 => 0,1
    if (sum(maternal_GT) == 1) and (sum(child_GT) == 1):
        # can not distinguish the maternal allele
        return None, None, None

    h1, h2, h3 = None, None, None
    if sum(maternal_GT) == 0:
        h1, h2 = 0, 0
        h3 = 0 if sum(child_GT) == 0 else 1

    elif sum(maternal_GT) == 1:  # 01 or 10
        if sum(child_GT) == 0:   # could only be 0 or 2
            h1, h2, h3 = 0, 1, 0
        elif sum(child_GT) == 2:
            h1, h2, h3 = 1, 0, 1
        else:
            raise ValueError("[ERROR] Child Genotype error.")

    elif sum(maternal_GT) == 2:
        h1, h2 = 1, 1
        h3 = sum(child_GT) - h2  # child_GT) could only be [0,1]/[1,0] or [1,1]

    else:
        raise ValueError("[ERROR] Maternal genotype error!")

    return h1, h2, h3


def calculate_genetic_score(in_vcf_fn, pos_beta_value, child_mother_pairs):
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
            af = float(info["AF"])

            ind_format = {name: i for i, name in enumerate(col[8].split(":"))}
            if ("GT" not in ind_format) or ("GP" not in ind_format):
                raise ValueError("[ERROR]VCF ERROR: GT or GP not in FORMAT.")

            for m, c in mother_child_idx:
                # Genotype should be: ["0", "0"], ["0", "1"], ["1", "0"] or ["1", "1"]
                mother_gt = col[m].split(":")[ind_format["GT"]].replace("/", "|").split("|")
                child_gt = col[c].split(":")[ind_format["GT"]].replace("/", "|").split("|")

                # Should be the probability of genotype: [0.99, 0.01, 0.00]
                # mother_gp = list(map(float, col[m].split(":")[ind_format["GP"]].split(",")))
                # child_gp = list(map(float, col[c].split(":")[ind_format["GP"]].split(",")))
                if ("." in mother_gt) or ("." in child_gt):
                    continue

                mother_gt = list(map(int, mother_gt))  # [0, 0], [0, 1], [1, 0] or [1, 1]
                child_gt = list(map(int, child_gt))  # [0, 0], [0, 1], [1, 0] or [1, 1]

                # `h1`: maternal transmitted allele
                # `h2`: maternal non-transmitted allele
                # `h3`: paternal transmitted allele (fetal only allele)
                h1, h2, h3 = distinguish_allele(mother_gt, child_gt)
                if h1 is None:
                    continue

                # `s_h1`: maternal transmitted haplotype genetic score
                # `s_h2`: maternal non-transmitted haplotype genetic score
                # `s_h3`: paternal (fetal only) transmitted haplotype genetic score
                # `s_mat`: maternal genotype score
                # `s_fet`: fetal genotype score
                s_h1 = h1 * beta
                s_h2 = h2 * beta
                s_h3 = h3 * beta
                s_mat = s_h1 + s_h2
                s_fet = s_h1 + s_h3

                k = index2sample[m] + "_" + index2sample[c]
                if k not in gs:
                    gs[k] = []

                gs[k].append([s_mat, s_fet, s_h1, s_h2, s_h3])

    elapse_time = datetime.now() - START_TIME
    sys.stderr.write("[INFO] All %d records loaded, %d seconds elapsed.\n" % (n, elapse_time.seconds))

    # calculate the PRS for each type of allele
    print("#Mother\tChild\tmaternal_genotype_score\tchild_genotype_score\th1\th2\th3\tsite_number")
    for m, c in mother_child_idx:
        k = index2sample[m] + "_" + index2sample[c]
        score = np.mean(gs[k], axis=0)
        print("%s\t%s\t%s\t%d" % (index2sample[m], index2sample[c], "\t".join(map(str, score)), len(gs[k])))

    return


def phenotype_concat(in_prs_fn, in_pheno_file):
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
            for i in range(len(col)):
                if col[i] == "-9" or col[i] == "-9.0":
                    col[i] = "NA"

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
    cmd_parser = argparse.ArgumentParser(description="Usage: ")
    commands = cmd_parser.add_subparsers(dest="command", title="Commands")

    gs_cmd = commands.add_parser("GeneticScore", help="Calculate Genetic score")
    gs_cmd.add_argument("-I", "--target", dest="target", type=str, required=True,
                        help="Input VCF. Required.")
    gs_cmd.add_argument("-b", "--base", dest="base", type=str, required=True,
                        help="A POS file with beta value for each position")
    gs_cmd.add_argument("--fam", dest="fam", type=str, required=True,
                        help="Input a .fam file with mother and children.")

    mr_cmd = commands.add_parser("MR", help="Mendelian Randomization")
    mr_cmd.add_argument("-I", "--input", dest="input", type=str, required=True,
                        help="Input data file")
    mr_cmd.add_argument("-y", dest="y_name", type=str, required=True,
                        help="Load the designated data as y from the '--input'.")
    mr_cmd.add_argument("-x", dest="x_name", type=str, required=True,
                        help="Load the designated phenotype(s) as x from the '--input'.")
    mr_cmd.add_argument("--covar", dest="covar_name", type=str, required=True,
                        help="Only load the designated covariate(s) from the '--input'.")

    add_cmd = commands.add_parser("ADD", help="Concat genetic score data together with phenotype data.")
    add_cmd.add_argument("-g", dest="genetic_score", type=str, required=True, help="input genetic score file.")
    add_cmd.add_argument("-p", dest="phenotype", type=str, required=True, help="input phenotype data file.")

    if len(sys.argv) < 2:
        cmd_parser.print_help()
        sys.exit(1)

    args = cmd_parser.parse_args()

    if args.command == "GeneticScore":
        child_mother_pairs = get_child_mother_duos(args.fam)
        beta_value = get_beta_value(args.base)
        calculate_genetic_score(args.target, beta_value, child_mother_pairs)

    elif args.command == "MR":
        data = pd.read_table(args.input, sep="\t")
        mendelian_randomization(data, args.y_name, args.x_name, args.covar_name)

    elif args.command == "ADD":
        phenotype_concat(args.genetic_score, args.phenotype)

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write("\n** process done, %d seconds elapsed **\n" % elapsed_time.seconds)
