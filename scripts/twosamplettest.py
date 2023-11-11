"""Two sample t-test
"""
import argparse
import sys

import numpy as np
import scipy as sp


if __name__ == "__main__":
    # cmdparser = argparse.ArgumentParser(description="Perform a two sample t-test.")
    # cmdparser.add_argument("sample1", type=str, required=True,
                           # help="Statistic data of sample1. input format: mean,SE,samplesize")
    # cmdparser.add_argument("sample2", type=str, required=True,
                           # help="Statistic data of sample2. input format: mean,SE,samplesize")
    # cmdparser.add_argument("--two-tail-test", dest="is_two_tail", action="store_true",
                           # help="Set this option to perform two tail two-sample t-test.")
    # args = cmdparser.parse_args()
    # mean1, se1, n1 = list(map(float, args.sample1.split(",")))
    # mean2, se2, n2 = list(map(float, args.sample2.split(",")))

    mean1, se1, n1 = list(map(float, sys.argv[1].split(",")))
    mean2, se2, n2 = list(map(float, sys.argv[2].split(",")))
    is_two_tail = sys.argv[3]


    # Null hypothesis (H0): mean1 == mean2
    # Alternative hypothesis (H1): mean1 != mean2
    delta_mean = mean1 - mean2
    delta_std = np.sqrt(se1**2 + se2**2)

    # Computed T statistic
    t_value = delta_mean / delta_std

    # adjusted the degrees of freedom
    df = (se1**2 + se2**2)**2 / (se1**4/(n1 - 1) + se2**4/(n2 - 1))

    one_tail_p_value = sp.stats.t.cdf(x=t_value, df=df)
    if one_tail_p_value > 0.5:
        one_tail_p_value = 1 - one_tail_p_value


    p_value = 2 * one_tail_p_value if is_two_tail == '2' else one_tail_p_value
    print(f"DD: {delta_mean}\t{delta_std}\t{t_value}\t{df}\t{p_value}")

    print(f"delta_mean: {delta_mean}")
    print(f"standard deviation: {delta_std}")
    print(f"T value: {t_value}")
    print(f"Degree of freedom: {df}")  
    # print(f"CI: {delta_mean - 2 * delta_std}, {delta_mean + 2 * delta_std}")
    print(f"sample1: {mean1}, {se1}, {n1}")
    print(f"sample2: {mean2}, {se2}, {n2}")
    print(f"Two-sample t-test: {p_value}\n")




