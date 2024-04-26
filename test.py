from ete3 import Tree
import os


fn = "data/results_bootstrapping/iqtree/felekesemitic_lexibank_cognate_full/part_bootstrap.splits.nex"

with open(fn, "r") as my_file:
    lines = my_file.readlines()

i = 0
while not lines[i].startswith("MATRIX"):
    i += 1
i += 1
supports = []
while i < len(lines):
    supports.append((int(lines[i].split("\t")[1])) / 100.0)
    i += 1
return supports
