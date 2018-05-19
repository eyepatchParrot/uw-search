#!/usr/bin/env python

from os import popen
from subprocess import call
from pandas import read_csv

call(['make', 'release'])
run_param = ['algorithm', 'n', 'param', 'distribution', 'record', 'n_thds']
df = read_csv(popen('./search i-seq.tsv'), sep='\t')
print(df.groupby(run_param)['ns'].describe())
