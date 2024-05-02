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
import seaborn



def run_inference(msa_path, model, prefix, s = False):
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    args = ""
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(prefix + ".raxml.bestTree"):
        args = args + " --redo"
    if s:
        command = "./bin/raxml-ng-single-noforce"
    else:
        command = "./bin/raxml-ng-multiple-force"
    command += " --msa " + msa_path
    command += " --model " + model
    command += " --prefix " + prefix
    command += " --threads auto --seed 2 --force model_lh_impr -blopt nr_safe"
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


def substitution_rates(prefix, x):
    file_path = prefix + ".raxml.log"
    if not os.path.isfile(file_path):
        return []
    with open(file_path, "r") as logfile:
        lines = logfile.readlines()
    rates = []
    for line in lines:
        if line.startswith("   Substitution rates"):
            parts = line.split(": ")[1].split(" ")[:-1]
            rates = [float(part) for part in parts]
            break
    if rates == []:
        print("Empty rates")
        return rates
    if x == -1: #code for single lamda rate
        return [rates[0], rates[1]]
    elif x == 4:
        return [rates[0], rates[1]]
    elif x == 8:
        return [rates[0], rates[1], rates[14]]
    elif x == 16:
        return [rates[0], rates[1], rates[30], rates[76]]
    elif x == 32:
        return [rates[0], rates[1], rates[62], rates[172], rates[344]]
    elif x == 64:
        return [rates[0], rates[1], rates[126], rates[364], rates[792], rates[1456]]
    else:
        print("Illegal x")
        return []

def base_frequencies(prefix, x):
    file_path = prefix + ".raxml.log"
    if not os.path.isfile(file_path):
        return []
    with open(file_path, "r") as logfile:
        lines = logfile.readlines()
    frequencies = []
    for line in lines:
        if line.startswith("   Base frequencies"):
            parts = line.split(": ")[1].split(" ")[:-1]
            frequencies = [float(part) for part in parts]
            break
    if frequencies == []:
        print("Empty frequences")
        return []
    elif x == 4:
        return [frequencies[0], frequencies[2]]
    elif x == 8:
        return [frequencies[0], frequencies[2], frequencies[6]]
    elif x == 16:
        return [frequencies[0], frequencies[2], frequencies[6], frequencies[14]]
    elif x == 32:
        return [frequencies[0], frequencies[2], frequencies[6], frequencies[14], frequencies[30]]
    elif x == 64:
        return [frequencies[0], frequencies[2], frequencies[6], frequencies[14], frequencies[30], frequencies[62]]
    else:
        print("Illegal x")
        return []


def sites_taxa_ratio(msa_path):
    with open(msa_path, "r") as msa_file:
        parts = msa_file.readlines()[0].split(" ")
    return int(parts[2]) / int(parts[1])



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
    prototype_s_prefix = os.path.join(target_dir, prototype_msa_type + "_COGs")
    run_inference(prototype_msa_path, "COG" + str(x), prototype_s_prefix, s = True)
    gtr_prefix = os.path.join(target_dir, prototype_msa_type + "_GTR")
    run_inference(prototype_msa_path, "MULTI" + str(x - 1) + "_GTR", gtr_prefix)
    mk_prefix = os.path.join(target_dir, prototype_msa_type + "_MK")
    run_inference(prototype_msa_path, "MULTI" + str(x - 1) + "_MK", mk_prefix)




def AIC_scores(target_dir, kappa):
    results = []
    bin_msa_type = "bin_part_" + str(kappa)
    prototype_msa_type = "prototype_part_" + str(kappa)
    for m, (model, msa_type) in enumerate([("BIN", bin_msa_type), ("COG", prototype_msa_type), ("COGs", prototype_msa_type), ("GTR", prototype_msa_type), ("MK", prototype_msa_type)]):
        prefix = os.path.join(target_dir, msa_type + "_" + model)
        results.append(AIC(prefix))
    return results



def get_all_substitution_rates(raxmlng_super_dir, kappa, s = False):
    r = []
    for ds_name in os.listdir(raxmlng_super_dir):
        target_dir = os.path.join(raxmlng_super_dir, ds_name)
        prototype_msa_type = "prototype_part_" + str(kappa)
        prefix = os.path.join(target_dir, prototype_msa_type + "_COG")
        if s:
            prefix += "s"
        if s:
            x = -1
        else:
            x = int(math.pow(2, kappa))
        f = substitution_rates(prefix, x) 
        #if len(f) == 0:
        #    continue
        r.append(f)
    return r

def get_all_base_frequencies(raxmlng_super_dir, kappa, s = False):
    r = []
    for ds_name in os.listdir(raxmlng_super_dir):
        target_dir = os.path.join(raxmlng_super_dir, ds_name)
        prototype_msa_type = "prototype_part_" + str(kappa)
        prefix = os.path.join(target_dir, prototype_msa_type + "_COG")
        if s:
            prefix += "s"
        x = int(math.pow(2, kappa))
        f = base_frequencies(prefix, x)
        #if len(f) == 0:
        #    continue
        r.append(f)
    return r


def violin_plots(results, path):
    models = ["BIN", "COG", "COGs", "GTR", "MK"]
    results_transformed = [[] for _ in range(5)]
    for row in results:
        if row[1] != row[1]:
            continue
        for i in range(1, 6):
            results_transformed[i-1].append(row[i])
    ax = seaborn.violinplot(data = results_transformed)
    ax.set_xticklabels(models)
    plt.savefig(path)
    plt.clf()
    plt.close()



def rates_stacked_plot(all_rates, path, plot_type):
    all_rates = [rates for rates in all_rates if rates != []]
    if plot_type == "bf":
        all_rates = [[0] + rates for rates in all_rates] #because of colors
    max_num = max([len(rates) for rates in all_rates])
    fig,ax = plt.subplots()
    x = range(len(all_rates))
    y_old = [0 for el in x]
    if plot_type == "bf":
        label_list = ['dummy', r'$\pi_1$', r'$\pi_2$', r'$\pi_3$', r'$\pi_4$', r'$\pi_5$', r'$\pi_6$']
    elif plot_type == "sr":
        label_list = [r'$\lambda_0$', r'$\lambda_1$', r'$\lambda_2$', r'$\lambda_3$', \
                r'$\lambda_4$', r'$\lambda_5$', r'$\lambda_6$']
    else:
        print("Illegal plot type")
        return
    for num in range(max_num):
        y_new = []
        for rates in all_rates:
            if len(rates) > num:
                y_new.append(rates[num] / sum(rates))
            else:
                y_new.append(0)
        if plot_type == "bf" and num == 0:
            ax.bar(x, y_new, bottom=y_old)
        else:
            ax.bar(x, y_new, bottom=y_old, label = label_list[num])
        for i in x:
            y_old[i] = y_old[i] + y_new[i]
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=max_num)
    plt.tick_params(labelbottom = False, bottom = False)
    plt.xlabel("Datasets")
    if plot_type == "bf":
        plt.ylabel("Base Frequencies (relative)")
    if plot_type == "sr":
        plt.ylabel("Substitution Rates (relative)")
    plt.savefig(path)
    plt.clf()
    plt.close()


def raxml_ng():
    for ds_name in os.listdir(msa_super_dir):
        msa_dir = os.path.join(msa_super_dir, ds_name)
        target_dir = os.path.join(raxmlng_super_dir, ds_name)
        bin_msa_type = "bin_part_" + str(kappa)
        prototype_msa_type = "prototype_part_" + str(kappa)
        bin_msa_path = os.path.join(msa_dir, bin_msa_type + ".phy")
        prototype_msa_path = os.path.join(msa_dir, prototype_msa_type + ".phy")
        if not os.path.isfile(bin_msa_path) or not os.path.isfile(prototype_msa_path):
            continue
        run_raxml_ng(msa_dir, target_dir, kappa)

def AIC_analysis():
    AIC_res = []
    for ds_name in os.listdir(msa_super_dir):
        msa_dir = os.path.join(msa_super_dir, ds_name)
        target_dir = os.path.join(raxmlng_super_dir, ds_name)
        bin_msa_type = "bin_part_" + str(kappa)
        prototype_msa_type = "prototype_part_" + str(kappa)
        bin_msa_path = os.path.join(msa_dir, bin_msa_type + ".phy")
        prototype_msa_path = os.path.join(msa_dir, prototype_msa_type + ".phy")
        if not os.path.isfile(bin_msa_path) or not os.path.isfile(prototype_msa_path):
            continue
        AIC_res.append([ds_name] + AIC_scores(target_dir, kappa))
    violin_plots(AIC_res, os.path.join(plots_super_dir, "AIC_" + str(kappa) + ".png"))
    print(tabulate(AIC_res, tablefmt="pipe", headers = ["dataset", "BIN", "COG", "COGs", "GTR", "MK"]))


version = "multiple-force"
msa_super_dir = "data/lingdata_cognate/msa"
plots_super_dir = os.path.join("data", "plots")
raxmlng_super_dir = os.path.join("data","inferences")
if not os.path.isdir(plots_super_dir):
    os.makedirs(plots_super_dir)
kappa = 4 

raxml_ng()
AIC_analysis()

rates_stacked_plot(get_all_substitution_rates(raxmlng_super_dir, kappa), os.path.join(plots_super_dir, "substitution_rates_" + str(kappa) + ".png"), "sr")
rates_stacked_plot(get_all_substitution_rates(raxmlng_super_dir, kappa, True), os.path.join(plots_super_dir, "substitution_rates_" + str(kappa) + "_s.png"), "sr")
rates_stacked_plot(get_all_base_frequencies(raxmlng_super_dir, kappa), os.path.join(plots_super_dir, "base_frequencies_" + str(kappa) + ".png"), "bf")
