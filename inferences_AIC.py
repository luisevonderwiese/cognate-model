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
import matplotlib




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
    command += " --threads auto --seed 2 --force model_lh_impr -blopt nr_safe"
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

def relative_llh(msa_path, prefix, kappa, model):
    with open(msa_path, "r") as msa_file:
        num_sites = int(msa_file.readlines()[0].split(" ")[2])
    if model == "BIN":
        num_sites = int(num_sites / kappa)
    return final_llh(prefix) / num_sites
    #return final_llh(prefix)


def AIC(prefix):
    logpath = prefix + ".raxml.log"
    if not os.path.isfile(logpath):
        return float('nan')
    with open(logpath, "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("AIC"):
            parts = line.split(" / ")
            scores = []
            for part in parts:
                scores.append(float(part.split(" ")[2]))
            return scores[0]
    return float('nan')



def run_raxml_ng(msa_dir, target_dir, kappa):
    bin_msa_type = "bin_part_" + str(kappa)
    prototype_msa_type = "prototype_part_" + str(kappa)
    bin_msa_path = os.path.join(msa_dir, bin_msa_type + ".phy")
    prototype_msa_path = os.path.join(msa_dir, prototype_msa_type + ".phy")
    x = int(math.pow(2, kappa))
    bin_prefix = os.path.join(target_dir, bin_msa_type + "_BIN")
    run_inference(bin_msa_path, "BIN", bin_prefix)
    prototype_prefix = os.path.join(target_dir, prototype_msa_type + "_COG")
    run_inference(prototype_msa_path, "COG" + str(x), prototype_prefix)
    gtr_prefix = os.path.join(target_dir, prototype_msa_type + "_GTR")
    run_inference(prototype_msa_path, "MULTI" + str(x - 1) + "_GTR", gtr_prefix)
    mk_prefix = os.path.join(target_dir, prototype_msa_type + "_MK")
    run_inference(prototype_msa_path, "MULTI" + str(x - 1) + "_MK", mk_prefix)




def AIC_analysis(target_dir, kappa, cv = False):
    results = []
    bin_msa_type = "bin_part_" + str(kappa)
    prototype_msa_type = "prototype_part_" + str(kappa)
    for m, (model, msa_type) in enumerate([("BIN", bin_msa_type), ("COG", prototype_msa_type), ("GTR", prototype_msa_type), ("MK", prototype_msa_type)]):
        if cv:
            prefix = os.path.join(target_dir, msa_type + "_cv_train_0_" + model)
        else:
            prefix = os.path.join(target_dir, msa_type + "_" + model)
        results.append(AIC(prefix))
    return results

def llh_analysis(target_dir, kappa, cv = False):
    results = []
    bin_msa_type = "bin_part_" + str(kappa)
    prototype_msa_type = "prototype_part_" + str(kappa)
    for m, (model, msa_type) in enumerate([("BIN", bin_msa_type), ("COG", prototype_msa_type), ("GTR", prototype_msa_type), ("MK", prototype_msa_type)]):
        if cv:
            prefix = os.path.join(target_dir, msa_type + "_cv_train_0_" +  model)
        else:
            prefix = os.path.join(target_dir, msa_type + "_" +  model)
        results.append(final_llh(prefix))
    return results

def violin_plots(AIC_results, path):
    results_transformed = [[] for _ in range(4)]
    for row in AIC_results:
        for i in range(1, 5):
            results_transformed[i-1].append(row[i])
    fig, ax = plt.figure()
    ax.violinplot(results_transformed)
    plt.savefig(path)
    plt.clf()
    plt.close()



msa_super_dir = "data/lingdata_cognate/msa"
raxmlng_super_dir = "data/inferences"
cv_super_dir = "data/cross_validation"
plots_super_dir = "data/AIC_plots"
if not os.path.isdir(plots_super_dir):
    os.makedirs(plots_super_dir)
kappa = 3
random.seed(2)
AIC_res = []
AIC_cv_res = []
llh_res = []
AIC_headers = ("dataset", "AIC BIN", "AIC COG", "AIC GTR", "AIC MK")
llh_headers = ("dataset", "llh BIN", "llh COG", "llh GTR", "llh MK")
for ds_name in os.listdir(msa_super_dir):
    msa_dir = os.path.join(msa_super_dir, ds_name)
    target_dir = os.path.join(raxmlng_super_dir, ds_name)
    cv_target_dir = os.path.join(cv_super_dir, ds_name)
    bin_msa_type = "bin_part_" + str(kappa)
    prototype_msa_type = "prototype_part_" + str(kappa)
    bin_msa_path = os.path.join(msa_dir, bin_msa_type + ".phy")
    prototype_msa_path = os.path.join(msa_dir, prototype_msa_type + ".phy")
    if not os.path.isfile(bin_msa_path) or not os.path.isfile(prototype_msa_path):
        continue
    #run_raxml_ng(msa_dir, target_dir, kappa)
    AIC_res.append([ds_name] + AIC_analysis(target_dir, kappa))
    AIC_cv_res.append([ds_name] + AIC_analysis(cv_target_dir, kappa, true))
    #llh_res.append([ds_name] + llh_analysis(target_dir, kappa))
violin_plots(AIC_res, os.path.join(plots_super_dir, str(kappa) + "_all.png"))
violin_plots(AIC_cv_res, os.path.join(plots_super_dir, str(kappa) + "_cv.png"))
print(tabulate(AIC_res, tablefmt="pipe", headers = AIC_headers))
print(tabulate(AIC_cv_res, tablefmt="pipe", headers = AIC_headers))
