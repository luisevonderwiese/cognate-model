from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter
import random
import os
import math
import numpy as np
from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import matplotlib
import seaborn
import pandas as pd
import statistics

def run_inference(msa_path, model, prefix, args = ""):
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(prefix + ".raxml.bestTree"):
        args = args + " --redo"
    command = "./bin/raxml-ng-COG"
    command += " --msa " + msa_path
    command += " --model " + model
    command += " --prefix " + prefix
    command += " --threads auto --seed 2 --force model_lh_impr -blopt nr_safe"
    command += " " + args
    os.system(command)


def run_evaluate(msa_path, prefix, ref_prefix, args = ""):
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(ref_prefix + ".raxml.bestModel"):
        return
    with open(ref_prefix + ".raxml.bestModel", "r") as model_file:
        model =  model_file.readlines()[0].split(",")[0]
    command = "./bin/raxml-ng-COG --evaluate "
    command += " --msa " + msa_path
    command += " --tree " + ref_prefix + ".raxml.bestTree"
    command += " --model " + model
    command += " --prefix " + prefix
    command += " --threads auto --seed 2 --opt-model off --opt-branches off"
    command += " " + args
    os.system(command)


def final_llh(prefix):
    if not os.path.isfile(prefix + ".raxml.log"):
        return float("nan")
    with open(prefix + ".raxml.log", "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("Final LogLikelihood: "):
            return float(line.split(": ")[1])
    return float('nan')

def relative_llh(msa_path, prefix):
    with open(msa_path, "r") as msa_file:
        num_sites = int(msa_file.readlines()[0].split(" ")[2])
    return final_llh(prefix) / num_sites
    #return final_llh(prefix)


def train_raxml_ng(msa_dir, target_dir):
    for t in range(10):
        bin_msa_path = os.path.join(msa_dir, "bin_cv_train_" + str(t) + ".phy")
        bin_prefix = os.path.join(target_dir, "bin_cv_train_" + str(t) + "_BIN")
        run_inference(bin_msa_path, "BIN", bin_prefix)



def test_raxml_ng(msa_dir, target_dir):
    for t in range(10):
        bin_msa_path = os.path.join(msa_dir, "bin_cv_test_" + str(t) + ".phy")
        bin_prefix = os.path.join(target_dir, "bin_cv_train_" + str(t) + "_BIN")
        bin_test_prefix = os.path.join(target_dir, "bin_cv_test_" + str(t) + "_BIN")
        run_evaluate(bin_msa_path, bin_test_prefix, bin_prefix)



def analysis(msa_dir, target_dir):
    results = [[] for _ in range(2)]
    for t in range(10):
        for m, (model, msa_type) in enumerate([("BIN", "bin")]):
            train_msa_path = os.path.join(msa_dir, msa_type + "_cv_train_" + str(t) + ".phy")
            train_prefix = os.path.join(target_dir, msa_type + "_cv_train_" + str(t) + "_" + model)
            results[m * 2].append(relative_llh(train_msa_path, train_prefix))
            test_msa_path = os.path.join(msa_dir, msa_type + "_cv_test_" + str(t) + ".phy")
            test_prefix = os.path.join(target_dir, msa_type + "_cv_test_" + str(t) + "_" + model)
            results[m * 2 + 1].append(relative_llh(test_msa_path, test_prefix))
    return [sum(el) / len(el) for el in results]


def differences_analysis(msa_dir, target_dir):
    results = [[] for _ in range(1)]
    for t in range(10):
        for m, (model, msa_type) in enumerate([("BIN", "bin")]):
            train_msa_path = os.path.join(msa_dir, msa_type + "_cv_train_" + str(t) + ".phy")
            train_prefix = os.path.join(target_dir, msa_type + "_cv_train_" + str(t) + "_" + model)
            rel_train_llh = relative_llh(train_msa_path, train_prefix)
            test_msa_path = os.path.join(msa_dir, msa_type + "_cv_test_" + str(t) + ".phy")
            test_prefix = os.path.join(target_dir, msa_type + "_cv_test_" + str(t) + "_" + model)
            rel_test_llh = relative_llh(test_msa_path, test_prefix)
            results[m].append((rel_test_llh - rel_train_llh) / rel_train_llh)
            #results[m].append(rel_train_llh - rel_test_llh)
    return [sum(el) / len(el) for el in results]





def violin_plots(results, path):
    models = ["BIN"]
    results_transformed = [[] for _ in range(1)]
    for row in results:
        if row[1] != row[1]:
            continue
        for i in range(1, 2):
            results_transformed[i-1].append(row[i])
    print(models[0], str(statistics.median(results_transformed[0])))
    plt.rcParams["figure.figsize"] = (3.5,6)
    ax = seaborn.violinplot(data = results_transformed, palette = [cm.to_hex(plt.cm.Set2(num)) for num in range(5)])
    ax.set_xticklabels(models)
    #plt.ylabel(r"$e$ (average)")
    plt.savefig(path + "/violin.png")
    plt.clf()
    plt.close()


msa_super_dir = "data/lexibench/bin_cross_validation"
raxmlng_super_dir = "data/bin_cross_validation"
plots_super_dir = "data/bin_cross_validation_plots"
all_diff_res = []
diff_headers = ("dataset", "diff_BIN")
for ds_name in os.listdir(msa_super_dir):
    msa_dir = os.path.join(msa_super_dir, ds_name)
    target_dir = os.path.join(raxmlng_super_dir, ds_name)
    #train_raxml_ng(msa_dir, target_dir)
    #test_raxml_ng(msa_dir, target_dir)
    all_diff_res.append([ds_name] + differences_analysis(msa_dir, target_dir))
violin_plots(all_diff_res, plots_super_dir)
print(tabulate(all_diff_res, tablefmt="pipe", headers = diff_headers))
