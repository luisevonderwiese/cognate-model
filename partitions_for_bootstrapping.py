import os
import json
import math
from lingdata import database
from ete3 import Tree
from tabulate import tabulate
from Bio import AlignIO



def run_iqtree(msa_path, part_path, prefix):
    if os.path.isfile(prefix + ".splits.nex"):
        return
    command = "./bin/iqtree2"
    command += " -s " + msa_path
    command += " -p " + part_path
    command += " --prefix " + prefix
    command +=" -B 1000 --sampling GENE"
    os.system(command)

def run_iqtree_gtr(msa_path, prefix):
    #if os.path.isfile(prefix + ".splits.nex"):
    #            return
    command = "./bin/iqtree2"
    command += " -s " + msa_path
    command += " -st MORPH -m GTR"
    command += " --prefix " + prefix
    command +=" -B 1000 --redo"
    os.system(command)



def get_support_values(prefix):
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


def create_partitions(df, msa_type):
    kappa = int(msa_type.split("_")[-1])
    for i, row in df.iterrows():
        msa_path = row["msa_paths"][msa_type]
        if not os.path.isfile(msa_path):
            continue
        part_subdir = os.path.join(part_dir, "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]]))
        if not os.path.isdir(part_subdir):
            os.makedirs(part_subdir)
        part_path = os.path.join(part_subdir, msa_type + ".part")
        with open(msa_path, "r") as msa_file:
            num_sites = int(msa_file.readlines()[0].split(" ")[2])
        assert(num_sites % kappa == 0)
        num_parts = int(num_sites / kappa)
        with open(part_path, "w+") as part_file:
            for i in range(num_parts):
                start_idx = i * kappa + 1
                end_idx = (i + 1) * kappa
                part_idx = i + 1
                s = "BIN,  p" + str(part_idx) + "=" + str(start_idx) + "-" + str(end_idx) + "\n"
                part_file.write(s)


def bootstrap_iqtree(df, msa_type):
    kappa = int(msa_type.split("_")[-1])
    for i, row in df.iterrows():
        msa_path = row["msa_paths"][msa_type]
        if not os.path.isfile(msa_path):
            continue
        part_subdir = os.path.join(part_dir, "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]]))
        part_path = os.path.join(part_subdir, msa_type + ".part")
        iqtree_subdir = os.path.join(iqtree_dir, "BIN", "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]]))
        if not os.path.isdir(iqtree_subdir):
            os.makedirs(iqtree_subdir)
        print(iqtree_subdir)
        prefix = os.path.join(iqtree_subdir, "part_bootstrap")
        run_iqtree(msa_path, part_path, prefix)

def bootstrap_iqtree_gtr(df, msa_type):
    kappa = int(msa_type.split("_")[-1])
    for i, row in df.iterrows():
        msa_path = row["msa_paths"][msa_type]
        if not os.path.isfile(msa_path):
            continue
        iqtree_subdir = os.path.join(iqtree_dir, "GTR", "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]]))
        if not os.path.isdir(iqtree_subdir):
            os.makedirs(iqtree_subdir)
        print(iqtree_subdir)
        prefix = os.path.join(iqtree_subdir, "bootstrap")
        run_iqtree_gtr(msa_path, prefix)

def analyze(df, msa_type):
    kappa = int(msa_type.split("_")[-1])
    r = []
    headers = ["ds_id", "mean support"]
    for i, row in df.iterrows():
        msa_path = row["msa_paths"][msa_type]
        if not os.path.isfile(msa_path):
            continue
        iqtree_subdir = os.path.join(iqtree_dir, "BIN", "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]]))
        prefix = os.path.join(iqtree_subdir, "part_bootstrap")
        vals = get_support_values(prefix)
        if len(vals) == 0:
            continue
        r.append([row["ds_id"], sum(vals) / len(vals)])
    print(tabulate(r, tablefmt="pipe", headers = headers))

def analyze_gtr(df, msa_type):
    kappa = int(msa_type.split("_")[-1])
    r = []
    headers = ["ds_id", "mean support"]
    for i, row in df.iterrows():
        msa_path = row["msa_paths"][msa_type]
        if not os.path.isfile(msa_path):
            continue
        iqtree_subdir = os.path.join(iqtree_dir, "GTR", "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]]))
        prefix = os.path.join(iqtree_subdir, "bootstrap")
        vals = get_support_values(prefix)
        if len(vals) == 0:
            continue
        r.append([row["ds_id"], sum(vals) / len(vals)])
    print(tabulate(r, tablefmt="pipe", headers = headers))




database.read_config("cognate_lingdata_config.json")
#database.compile()
df = database.data()
part_dir = "data/lingdata_cognate/partitions"
iqtree_dir =  "data/results_bootstrapping/iqtree"
msa_type = "bin_part_3"
create_partitions(df, msa_type)
#bootstrap_iqtree(df, msa_type)
#analyze(df, msa_type)
msa_type = "prototype_part_3"
bootstrap_iqtree_gtr(df, msa_type)
analyze_gtr(df, msa_type)

