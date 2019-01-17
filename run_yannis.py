#!/usr/bin/env python

from os import popen 
from subprocess import call
from pandas import read_csv
import sys
import subprocess

call(['make', 'release'])
run_param = ['algorithm', 'n', 'param', 'distribution', 'record', 'n_thds']

with open("/mnt/tmp/log_file", "w") as log_file:
    subprocess.run(["./release","i-seq.tsv"], stdout=log_file)

df = read_csv("/mnt/tmp/log_file", sep='\t')

print(df.groupby(run_param)['ns'].describe(percentiles=[.25,0.75]))
