#!/usr/bin/python 

import sys

with open(sys.argv[1], "r") as f:
    firstline = f.readline()
    firstline = firstline.strip().split('\t')
    #next(f)
    for line in f:
        genes = line.strip().split('\t')
        for i in range(3):
            for j in range(3):
                print(genes[i] + "\t" + firstline[i] + "\t" + genes[j] + "\t" + firstline[j])
