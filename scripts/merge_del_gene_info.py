#!/usr/local/Cluster-Apps/python/2.7.10/bin/python

########################################################################
# Written by Olga Shamardina
#
# Auxiliary script to merge the genes annotation
#
########################################################################

import sys
import os.path

dels_file = sys.argv[1]
a = os.path.splitext(dels_file)
output_file = a[0] + "_merged" + a[1]

dels = {}
with open(dels_file, "r") as f, open(output_file, "w") as outf:
    line = f.readline().rstrip().split("\t")
    coord = line[0:3]
    genes = {line[3]}
    outf.write("\t".join(["CHROM", "START", "END", "GENES"]) + "\n")
    for l in f:
        line = l.rstrip().split("\t")
        coord_new = line[0:3]
        if coord == coord_new:
            genes.add(line[3])
        else:
            genes.discard(".")
            outf.write("\t".join(coord) + "\t" + ";".join(sorted(genes)) + "\n")
            genes = {line[3]}
        coord = coord_new
    genes.discard(".")
    outf.write("\t".join(coord) + "\t" + ";".join(sorted(genes)) + "\n")
