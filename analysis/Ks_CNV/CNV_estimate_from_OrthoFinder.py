#!/usr/bin/python

"""This takes in as input the OrthoGroups.csv from OrthoFinder and calculates how many gene loss
and gene gain events have happened in each lineage, as compared to the lineage in the first tab-delimited column. 
The species in the first column must share the same number of genes as at least one other species, but they cannot 
all have the same number (indicating no gain and no loss). It then calculates how many gene loss or gains occured
in each lineage. 

usage: .py OrthoGroups.csv <outfile>"""

import sys

infile = open(sys.argv[1], "r")
outfile = open(sys.argv[2], "w")

sys.stdout = outfile

with infile as f:
    first_line = f.readline().strip()
    print(first_line)
    first_line = first_line.split('\t')
    gain_loss = {key:{"gain":[], "loss":[]} for key in first_line[1:]}
    for line in f:
        line = line.strip().split('\t')[1:]
        while len(line) < len(first_line):
            line.append(" ")
        c = [len(i.split(',')) for i in line]
        for i in range(len(first_line)):
            if line[i] == "" or line[i] == " ":
                c[i] = 0 #corrects for no genes reporting has having length of 1 (empty string)
        if len(set(c)) != 1:
            if c[0] == c[1] and c[0] != 0:
                if c[1] > c[2]: #loss in c[2]
                    dif = c[1] - c[2]
                    gain_loss[first_line[2]]["loss"].append(dif)
#                    print(*line, sep='\t')
                elif c[2] > c[1]: #gain in c[2]
                    dif = c[2] - c[1]
                    gain_loss[first_line[2]]["gain"].append(dif)
#                    print(*line, sep='\t')
            elif c[0] == c[2] and c[0] != 0:
                if c[1] > c[2]: #gain in c[1]
                    dif = c[1]-c[2]
                    gain_loss[first_line[1]]["gain"].append(dif)
#                    print(*line, sep='\t')
                elif c[2] > c[1]: #loss in c[1]
                    dif = c[2]-c[1]
                    gain_loss[first_line[1]]["loss"].append(dif)
#                    print(*line, sep='\t')

for x in gain_loss:
    print(x)
    for y in gain_loss[x]:
        print('\t', y)
        for z in sorted(set(gain_loss[x][y])):
            print('\t\t',z,":",gain_loss[x][y].count(z))


