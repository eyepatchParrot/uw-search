#!/usr/bin/env python

from os import popen
from subprocess import call
from pandas import read_csv
import sys
import subprocess


lines = []
algorithms = ["b-lin","i-hyp-64", "b-eyt", "b-eyt-p"] #,"i-opt-8", "i-naive", "i-hyp-64"] #,"b-lin", "b-eyt", "b-eyt-p", "i-naive"]
recs = ["8"]

sizes = ["1000000"] #00000","1000000", "10000000"] #["1000","10000", "100000", "1000000", "10000000"] # "100000000", "1000000000"] 
#runs = {"1000":"10000","10000":"10000", "100000":"1000", "1000000":"1000", "10000000":"100"} #"100000000":"2", "1000000000":"1"}


dataset = "fal"
params = ["1.05", "1.25"] #["0.5","1.05","1.25","1.5"]
#params = ["/mnt/datasets/newman","/mnt/datasets/wiki"]
#param =  "1542497071"

for size in sizes:
#  times=runs[size]
  for alg in algorithms:
    for rec in recs:
      for param in params:
        s=size+"\t"+dataset+"\t"+param+"\t"+alg+"\t"+rec+"\t"+"1"
        for i in range(100): #range(int(times)):
          lines.append(s)




call(['rm', '-rf','i-seq.tsv'])
with open("i-seq.tsv", "w") as file:
  file.write("n\tdistribution\tparam\talgorithm\trecord\tthread\n")
  for line in lines:
    file.write(line+'\n')



