import os
import json
import math
from lingdata import database
from ete3 import Tree
from tabulate import tabulate
import matplotlib.pyplot as plt
import numpy as np



def support_values(prefix):
    if not os.path.isfile(prefix + ".raxml.support"):
        return []
    support_tree = Tree(prefix + ".raxml.support")
    #print([node.dist for node in support_tree.traverse("postorder") if node.dist > 1.0])
    supports = [node.dist / 100 for node in support_tree.traverse("postorder")]
    return supports

def support_values_iqtree(prefix):
    fn = prefix + ".splits.nex"
    if not os.path.isfile(fn):
        return []
    with open(fn, "r") as my_file:
        lines = my_file.readlines()

    i = 0
    while not lines[i].startswith("MATRIX"):
        i += 1
    i += 1
    supports = []
    while i < len(lines) - 22:
        supports.append((int(lines[i].split("\t")[1])) / 100.0)
        i += 1
    return supports


def analyze(df, configs):
    r = []
    headers = ["ds_id"] + [model for _,model in configs] + ["BIN_iqtree", "GTR_iqtree"]
    for i, row in df.iterrows():
        res_row = [row["ds_id"]]
        for msa_type, model in configs:
            results_dir = os.path.join(results_base_dir, model)
            msa_path = row["msa_paths"][msa_type]
            if not os.path.isfile(msa_path):
                break
            msa_string = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
            prefix = os.path.join(results_dir, msa_string, "support")
            s = support_values(prefix)
            res_row.append(np.mean(s))
        iqtree_subdir = os.path.join(iqtree_dir, "BIN", "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]]))
        prefix = os.path.join(iqtree_subdir, "part_bootstrap")
        s = support_values_iqtree(prefix)
        if len(s) != 0:
            res_row.append(np.mean(s))
        iqtree_subdir = os.path.join(iqtree_dir, "GTR", "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]]))
        prefix = os.path.join(iqtree_subdir, "bootstrap")
        s = support_values_iqtree(prefix)
        if len(s) != 0:
            res_row.append(np.mean(s))
        if len(res_row) == 6:
            r.append(res_row)
    print(model)
    print(tabulate(r, tablefmt="pipe", headers = headers))




database.read_config("cognate_lingdata_config.json")
#database.compile()
df = database.data()
results_base_dir = "data/results_bootstrapping"
iqtree_dir = os.path.join(results_base_dir, "iqtree")
configs = [("bin_part_3", "BIN"), ("prototype_part_3", "GTR"), ("prototype_part_3", "COG")]
analyze(df, configs)
