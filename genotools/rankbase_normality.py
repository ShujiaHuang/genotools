"""Ranks the trait values and scales the ranks to (0,1), then transform the scaled ranks 
to a normal distribution using an inverse normal transformation.

Author: Shujia Huang
Date: 2021-07-09 16:31:52
"""
import argparse
import numpy as np
import pandas as pd
from scipy.stats import norm


def rankbase_normality(x, rmoutliner=True, nan_data_mark=-9):
    """Ranks the x value and scales the ranks to (0,1), then transform the scaled ranks 
    to a normal distribution using an inverse normal transformation.
    """
    if rmoutliner:
        dd = [d for d in x if d != nan_data_mark]
        mean, std = np.mean(dd), np.std(dd)
        x = np.array([nan_data_mark if ((abs(d-mean)>3*std) or (d==nan_data_mark)) else d for d in x])
    else:
        x = np.array(x)
    
    e_idx = [i for i, d in enumerate(x) if d !=nan_data_mark]
    e2x = {i:j for i, j in enumerate(e_idx)}  # index of e map to the index of x
    
    # ranks the values by increase order and return the index of `x`
    rank_idx = [e2x[k] for k in np.argsort([x[i] for i in e_idx])]
    rank_size = len(rank_idx)

    # scale the rank value to (0, 1) and transform the scaled ranks to a normal distribution 
    # by using an inverse normal transformation.
    rank_scale = np.array([(i+1.0)/(rank_size+1) for i in range(rank_size)])
    iid = norm.ppf(rank_scale)
    
    for i, k in enumerate(rank_idx):
        x[k] = round(iid[i], 4)
    
    return x


if __name__ == "__main__":
    """
    python rankbase_normality.py -I in.tsv --phenotype "pheno1,pheno2,pheno3" -O out.tsv
    """

    desc = "Transform the trait value to a normal distribution using rank-base inverse normal transformation."
    cmdparser = argparse.ArgumentParser(description=desc)
    cmdparser.add_argument("-I", "--input", dest="input", type=str, required=True,
                           help="A phenotypt file. Required.")
    cmdparser.add_argument("--phenotype", dest="phenotype", type=str, required=True, 
                           help="Phenotypes' labels seperate by comma which need to be normalization.")
    cmdparser.add_argument("-O", "--output", dest="output", type=str, required=True,
                           help="output file path. Required.")
    cmdparser.add_argument("--donot_remove_outliner", dest="donot_remove_outliner", action="store_true",
                           help="Do not remove outliner value.")

    args = cmdparser.parse_args()
    df = pd.read_table(args.input, sep="\t")
    phenotypes = args.phenotype.strip().split(",")

    for f in phenotypes:
        if args.donot_remove_outliner:
            df[f] = pd.Series(rankbase_normality(df[f], rmoutliner=False, nan_data_mark=-9))
        else:
            df[f] = pd.Series(rankbase_normality(df[f], rmoutliner=True, nan_data_mark=-9))
    
    df.to_csv(args.output, sep="\t")




