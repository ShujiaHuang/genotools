## Detection of age-specific genetic effects on lipid levels among mothers and infants

We examined age-specific genetic influence on lipid levels between maternal and infant subjects with a two-sided two-sample t-test (see Code Availability) that tests the equivalence of genetic effects on the same traits between mothers and infants (Figure 5, Table S9) with the following hypotheses:

Null hypothesis H_0: β_m = β_i
Alternative hypothesis H_1: β_m ≠ β_i



```bash
#2.78e-9
python twosamplettest.py 0.4184,0.0677,1535 0.194675702,0.075442627,1457 2
python twosamplettest.py -0.4766,0.0594,1535 -0.102061516,0.06293341,1457 2
python twosamplettest.py 0.4151,0.0657,1532 0.3367,0.0652,1456 2
python twosamplettest.py -0.5984,0.0782,1532 -0.147362691,0.085472153,1455 2
python twosamplettest.py -0.7926,0.0606,1532 -0.708083103,0.061761567,1455 2
python twosamplettest.py -0.4559,0.0384,1531 0.088650799,0.042040878,1442 2

python twosamplettest.py 0.4214,0.0664,1456 0.420131273,0.06670078,1532 2
python twosamplettest.py 0.3372,0.0559,1455 -0.007205417,0.052879493,1532 2
python twosamplettest.py -0.7081,0.0618,1455 -0.792575853,0.06063332,1532 2
```

```bash
# 5e-8
python twosamplettest.py 0.4184,0.0677,1535 0.194675702,0.075442627,1457 2
python twosamplettest.py -0.4766,0.0594,1535 -0.102061516,0.06293341,1457 2
python twosamplettest.py 0.4151,0.0657,1532 0.3367,0.0652,1456 2
python twosamplettest.py -0.5984,0.0782,1532 -0.147362691,0.085472153,1455 2
python twosamplettest.py -0.7926,0.0606,1532 -0.708083103,0.061761567,1455 2
python twosamplettest.py -0.4559,0.0384,1531 0.088650799,0.042040878,1442 2

python twosamplettest.py 0.3371,0.0610,1456 0.0192,0.0610,1532 2
python twosamplettest.py 0.4214,0.0664,1456 0.4201,0.0667,1532 2
python twosamplettest.py 0.3372,0.0559,1455 -0.007205417,0.052879493,1532 2
python twosamplettest.py -0.7081,0.0618,1455 -0.792575853,0.06063332,1532 2
python twosamplettest.py -0.3597,0.0630,1442 -0.2297,0.0567,1531 2
python twosamplettest.py -0.2362,0.0422,1442 -0.0279,0.0389,1531 2
```


