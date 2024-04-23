import os
import json
import math
from lingdata import database
from ete3 import Tree
from tabulate import tabulate
import matplotlib.pyplot as plt
import numpy as np


def run_inference(msa_path, model, prefix, args = ""):
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(prefix + ".raxml.bestTree"):
        args = args + " --redo"
    command = "./bin/raxml-ng-multiple-force"
    command += " --msa " + msa_path
    command += " --model " + model
    command += " --prefix " + prefix
    command += " --threads auto --seed 2 --force model_lh_impr"
    command += " " + args
    os.system(command)


def run_bootstrap(msa_path, model, prefix, args = ""):
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(prefix + ".raxml.bootstraps"):
        args = args + " --redo"
    command = "./bin/raxml-ng-multiple-force --bootstrap"
    command += " --msa " + msa_path
    command += " --model " + model
    command += " --prefix " + prefix
    command += " --threads auto --seed 2 --force model_lh_impr"
    command += " " + args
    os.system(command)


def run_support(inference_prefix, bootstrap_prefix, prefix, args = ""):
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(prefix + ".raxml.support"):
        args = args + " --redo"
    command = "./bin/raxml-ng-multiple-force  --support"
    command += " --tree " + inference_prefix + ".raxml.bestTree"
    command += " --bs-trees " + bootstrap_prefix + ".raxml.bootstraps"
    command += " --prefix " + prefix
    command += " --threads auto  --seed 2 --force model_lh_impr"
    command += " " + args
    os.system(command)


def support_values(prefix):
    if not os.path.isfile(prefix + ".raxml.support"):
        return []
    support_tree = Tree(prefix + ".raxml.support")
    supports = [node.dist / 100 for node in support_tree.traverse("postorder")]
    return supports




def final_llh(prefix):
    if not os.path.isfile(prefix + ".raxml.log"):
        return float("nan")
    with open(prefix + ".raxml.log", "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("Final LogLikelihood: "):
            return float(line.split(": ")[1])
    return float('nan')


def inferences(df, msa_type, model):
    results_dir = os.path.join(results_base_dir, model)
    kappa = int(msa_type.split("_")[-1])
    x = int(math.pow(2, kappa))
    if model == "COG":
        model = "COG" + str(x)
    if model == "GTR":
        model = "MULTI" + str(x) + "_GTR"
    for i, row in df.iterrows():
        msa_path = row["msa_paths"][msa_type]
        if not os.path.isfile(msa_path):
            continue
        msa_string = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
        prefix = os.path.join(results_dir, msa_string, "inference")
        run_inference(msa_path, model, prefix)


def bootstraps(df, msa_type, model):
    results_dir = os.path.join(results_base_dir, model)
    kappa = int(msa_type.split("_")[-1])
    x = int(math.pow(2, kappa))
    if model == "COG":
        model = "COG" + str(x)
    if model == "GTR":
        model = "MULTI" + str(x) + "_GTR"
    for i, row in df.iterrows():
        msa_path = row["msa_paths"][msa_type]
        if not os.path.isfile(msa_path):
            continue
        msa_string = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
        prefix = os.path.join(results_dir, msa_string, "bootstrap")
        run_bootstrap(msa_path, model, prefix)


def supports(df, msa_type, model):
    results_dir = os.path.join(results_base_dir, model)
    for i, row in df.iterrows():
        msa_path = row["msa_paths"][msa_type]
        if not os.path.isfile(msa_path):
            continue
        msa_string = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
        bootstrap_prefix = os.path.join(results_dir, msa_string, "bootstrap")
        inference_prefix = os.path.join(results_dir, msa_string, "inference")
        prefix = os.path.join(results_dir, msa_string, "support")
        run_support(inference_prefix, bootstrap_prefix, prefix)

def analyze(df, msa_type, model):
    results_dir = os.path.join(results_base_dir, model)
    r = []
    headers = ["ds_id", "mean", "median", "var"]
    for i, row in df.iterrows():
        msa_path = row["msa_paths"][msa_type]
        if not os.path.isfile(msa_path):
            continue
        msa_string = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
        prefix = os.path.join(results_dir, msa_string, "support")
        s = support_values(prefix)
        r.append(row["ds_id"], s.mean(), s.median(), s.var())
    print(model)
    print(tabulate(r, tablefmt="pipe", headers = headers))


database.read_config("cognate_lingdata_config.json")
df = database.data()
results_base_dir = "data/results_bootstrapping"
for (msa_type, model) in [("bin_part_3", "BIN"), ("prototype_part_3", "GTR"), ("prototype_part_3", "COG")]:
    inferences(df, msa_type, model)
    bootstraps(df, msa_type, model)
    supports(df, msa_type, model)
    analyze(df, msa_type, model)
