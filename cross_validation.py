from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter
import random
import os
import math
import numpy as np
import pandas as pd
from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import matplotlib
import seaborn
import statistics

def split_indices(num_sites, ratio, num_samples = 10):
    num_sites_train = math.ceil(num_sites * ratio)
    res = []
    for _ in range(num_samples):
        li = [__ for __ in range(num_sites)]
        random.shuffle(li)
        train_indices = li[:num_sites_train]
        res.append(train_indices)
    return res


def empty_align(ref_align):
    new_records = [SeqRecord("", id=ref_align[i].id) for i in range(len(ref_align))]
    return MultipleSeqAlignment(new_records, annotations={}, column_annotations={})

def concat_align(a1, a2):
    new_sequences = []
    assert(len(a1) == len(a2))
    for i in range(len(a1)):
        seq1 = a1[i].seq
        seq2 = a2[i].seq
        new_sequences.append(seq1 + seq2)
    new_records = [SeqRecord(new_sequences[i], id=a1[i].id) for i in range(len(a1))]
    return MultipleSeqAlignment(new_records, annotations={}, column_annotations={})



def run_inference(msa_path, model, prefix, s = False):
    args = ""
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(prefix + ".raxml.bestTree"):
        args = args + " --redo"
    if s:
        command = "./bin/raxml-ng-COGs"
    else:
        command = "./bin/raxml-ng-COG"
    command += " --msa " + msa_path
    command += " --model " + model
    command += " --prefix " + prefix
    command += " --threads auto --seed 2 --force model_lh_impr -blopt nr_safe"
    command += " " + args
    os.system(command)


def run_evaluate(msa_path, prefix, ref_prefix, s = False):
    args = ""
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(ref_prefix + ".raxml.bestModel"):
        return
    with open(ref_prefix + ".raxml.bestModel", "r", encoding = "utf-8") as model_file:
        model =  model_file.readlines()[0].split(",")[0]
    if s:
        command = "./bin/raxml-ng-COGs --evaluate "
    else:
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
    with open(prefix + ".raxml.log", "r", encoding = "utf-8") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("Final LogLikelihood: "):
            return float(line.split(": ")[1])
    return float('nan')

def relative_llh(msa_path, prefix, kappa, model):
    with open(msa_path, "r", encoding = "utf-8") as msa_file:
        num_sites = int(msa_file.readlines()[0].split(" ")[2])
    if model == "BIN":
        num_sites = int(num_sites / kappa)
    return final_llh(prefix) / num_sites
    #return final_llh(prefix)

def create_samples(kappa, ratio, msa_dir, cv_msa_dir):
    bin_msa_type = "bin_part_" + str(kappa)
    bv_msa_type = "bv_part_" + str(kappa)
    bin_msa_path = os.path.join(msa_dir, bin_msa_type + ".phy")
    bv_msa_path = os.path.join(msa_dir, bv_msa_type + ".phy")
    with open(bv_msa_path, "r", encoding = "utf-8") as msa_file:
        num_sites = int(msa_file.readlines()[0].split(" ")[2])
    smaller_part_size = num_sites * (1 - ratio)
    if smaller_part_size < 10:
        return False
    with open(bin_msa_path, "r", encoding = "utf-8") as msa_file:
        num_sites_bin = int(msa_file.readlines()[0].split(" ")[2])
    assert(num_sites_bin == kappa * num_sites)
    try:
        bin_align = AlignIO.read(bin_msa_path, "phylip-relaxed")
    except:
        print(msa_dir, "Failed")
        return False
    try:
        bv_align = AlignIO.read(bv_msa_path, "phylip-relaxed")
    except:
        print(msa_dir, "Failed")
        return False
    if not os.path.isdir(cv_msa_dir):
        os.makedirs(os.path.join(cv_msa_dir, "train"))
        os.makedirs(os.path.join(cv_msa_dir, "test"))
    indices_list = split_indices(num_sites, ratio)
    for (t, train_indices) in enumerate(indices_list):
        bin_train_align = empty_align(bin_align)
        bin_test_align = empty_align(bin_align)
        bv_train_align = empty_align(bv_align)
        bv_test_align = empty_align(bv_align)
        for s in range(num_sites):
            if s in train_indices :
                bv_train_align = concat_align(bv_train_align, bv_align[:, s:s+1])
                bin_train_align = concat_align(bin_train_align, bin_align[:, s*kappa : (s+1) * kappa])
            else:
                bv_test_align = concat_align(bv_test_align, bv_align[:, s:s+1])
                bin_test_align = concat_align(bin_test_align, bin_align[:, s*kappa : (s+1) * kappa])
        with open(os.path.join(cv_msa_dir, "train", bin_msa_type +  str(t) + ".phy"),
                "w+", encoding = "utf-8") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(bin_train_align)
        with open(os.path.join(cv_msa_dir, "test", bin_msa_type + str(t) + ".phy"),
                "w+", encoding = "utf-8") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(bin_test_align)
        with open(os.path.join(cv_msa_dir, "train", bv_msa_type + str(t) + ".phy"),
                "w+", encoding = "utf-8") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(bv_train_align)
        with open(os.path.join(cv_msa_dir, "test", bv_msa_type + str(t) + ".phy"),
                "w+", encoding = "utf-8") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(bv_test_align)
    print(msa_dir, "done")
    return True


def train_raxml_ng(msa_dir, target_dir, kappa):
    bin_msa_type = "bin_part_" + str(kappa)
    bv_msa_type = "bv_part_" + str(kappa)
    x = int(math.pow(2, kappa))
    for t in range(10):
        bin_msa_path = os.path.join(msa_dir, "train", bin_msa_type + str(t) + ".phy")
        bin_prefix = os.path.join(target_dir,  "train", bin_msa_type + str(t),  "BIN")
        run_inference(bin_msa_path, "BIN", bin_prefix)
        bv_msa_path = os.path.join(msa_dir, "train", bv_msa_type + str(t) + ".phy")
        cog_prefix = os.path.join(target_dir, "train", bv_msa_type + str(t),  "COG")
        run_inference(bv_msa_path, "COG" + str(x), cog_prefix)
        cogs_prefix = os.path.join(target_dir, "train", bv_msa_type + str(t),  "COGs")
        run_inference(bv_msa_path, "COG" + str(x), cogs_prefix, s = True)
        gtr_prefix = os.path.join(target_dir, "train", bv_msa_type + str(t),  "GTR")
        run_inference(bv_msa_path, "MULTI" + str(x - 1) + "_GTR", gtr_prefix)
        mk_prefix = os.path.join(target_dir, "train", bv_msa_type + str(t),  "MK")
        run_inference(bv_msa_path, "MULTI" + str(x - 1) + "_MK", mk_prefix)



def test_raxml_ng(msa_dir, target_dir, kappa):
    bin_msa_type = "bin_part_" + str(kappa)
    bv_msa_type = "bv_part_" + str(kappa)
    for t in range(10):
        bin_msa_path = os.path.join(msa_dir, "test", bin_msa_type + str(t) + ".phy")
        bin_prefix = os.path.join(target_dir,  "train", bin_msa_type + str(t),  "BIN")
        bin_test_prefix = os.path.join(target_dir, "test", bin_msa_type + str(t),  "BIN")
        run_evaluate(bin_msa_path, bin_test_prefix, bin_prefix)
        bv_msa_path = os.path.join(msa_dir, "test", bv_msa_type + str(t) + ".phy")
        cog_prefix = os.path.join(target_dir, "train", bv_msa_type + str(t),  "COG")
        cog_test_prefix = os.path.join(target_dir, "test", bv_msa_type + str(t),  "COG")
        run_evaluate(bv_msa_path, cog_test_prefix, cog_prefix)
        cogs_prefix = os.path.join(target_dir, "train", bv_msa_type + str(t),  "COGs")
        cogs_test_prefix = os.path.join(target_dir, "test", bv_msa_type + str(t),  "COGs")
        run_evaluate(bv_msa_path, cogs_test_prefix, cogs_prefix, s = True)
        gtr_prefix = os.path.join(target_dir, "train", bv_msa_type + str(t),  "GTR")
        gtr_test_prefix = os.path.join(target_dir, "test", bv_msa_type + str(t),  "GTR")
        run_evaluate(bv_msa_path, gtr_test_prefix, gtr_prefix)
        mk_prefix = os.path.join(target_dir, "train", bv_msa_type + str(t),  "MK")
        mk_test_prefix = os.path.join(target_dir, "test", bv_msa_type + str(t),  "MK")
        run_evaluate(bv_msa_path, mk_test_prefix, mk_prefix)



def analysis(msa_dir, target_dir, kappa):
    bin_msa_type = "bin_part_" + str(kappa)
    bv_msa_type = "bv_part_" + str(kappa)
    results = [[] for _ in range(8)]
    for t in range(10):
        for m, (model, msa_type) in enumerate([("BIN", bin_msa_type), ("MK", bv_msa_type), ("GTR", bv_msa_type), ("COGs", bv_msa_type), ("COG", bv_msa_type)]):
            train_msa_path = os.path.join(msa_dir, "train", msa_type + str(t) + ".phy")
            train_prefix = os.path.join(target_dir, "train", msa_type +  str(t), model)
            results[m * 2].append(relative_llh(train_msa_path, train_prefix, kappa, model))
            test_msa_path = os.path.join(msa_dir, "test", msa_type + str(t) + ".phy")
            test_prefix = os.path.join(target_dir, "test", msa_type +  str(t), model)
            results[m * 2 + 1].append(relative_llh(test_msa_path, test_prefix, kappa, model))
    return [sum(el) / len(el) for el in results]


def differences_analysis(msa_dir, target_dir, kappa):
    bin_msa_type = "bin_part_" + str(kappa)
    bv_msa_type = "bv_part_" + str(kappa)
    results = [[] for _ in range(5)]
    for t in range(10):
        for m, (model, msa_type) in enumerate([("BIN", bin_msa_type), ("MK", bv_msa_type), ("GTR", bv_msa_type), ("COGs", bv_msa_type), ("COG", bv_msa_type)]):
            train_msa_path = os.path.join(msa_dir, "train", msa_type + str(t) + ".phy")
            train_prefix = os.path.join(target_dir, "train", msa_type +  str(t), model)
            rel_train_llh = relative_llh(train_msa_path, train_prefix, kappa, model)
            test_msa_path = os.path.join(msa_dir, "test", msa_type + str(t) + ".phy")
            test_prefix = os.path.join(target_dir, "test", msa_type +  str(t), model)
            rel_test_llh = relative_llh(test_msa_path, test_prefix, kappa, model)
            results[m].append((rel_test_llh - rel_train_llh) / rel_train_llh)
            #results[m].append(rel_train_llh - rel_test_llh)
    rv = [sum(el) / len(el) for el in results]
    return rv

def plots(msa_dir, target_dir, kappa, plots_super_dir, ds_name):
    bin_msa_type = "bin_part_" + str(kappa)
    bv_msa_type = "bv_part_" + str(kappa)
    ind = np.arange(10)
    width = 0.08
    offsets = [-0.36, -0.28, -0.20, -0.12, -0.04, 0.04, 0.12, 0.20, 0.28, 0.36]
    cmap_train = matplotlib.cm.get_cmap('Set1')
    cmap_test = matplotlib.cm.get_cmap('Pastel1')
    bin_msa_type = "bin_part_" + str(kappa)
    bv_msa_type = "bv_part_" + str(kappa)
    results = [[] for _ in range(10)]
    plots_dir = os.path.join(plots_super_dir, str(kappa))
    if not os.path.isdir(plots_dir):
        os.makedirs(plots_dir)
    for t in range(10):
        for m, (model, msa_type) in enumerate([("BIN", bin_msa_type), ("COG", bv_msa_type), ("COGs", bv_msa_type), ("GTR", bv_msa_type), ("MK", bv_msa_type)]):
            train_msa_path = os.path.join(msa_dir, "train", msa_type + str(t) + ".phy")
            train_prefix = os.path.join(target_dir, "train", msa_type +  str(t), model)
            results[m * 2].append(relative_llh(train_msa_path, train_prefix, kappa, model))
            test_msa_path = os.path.join(msa_dir, "test", msa_type + str(t) + ".phy")
            test_prefix = os.path.join(target_dir, "test", msa_type +  str(t), model)
            results[m * 2 + 1].append(relative_llh(test_msa_path, test_prefix, kappa, model))
    _, ax = plt.subplots()
    ax.bar(ind + offsets[0], results[0], width, label='train BIN', color = cmap_train(0))
    ax.bar(ind + offsets[1], results[1], width, label='test BIN', color = cmap_test(0))
    ax.bar(ind + offsets[2], results[2], width, label='train COG', color = cmap_train(1))
    ax.bar(ind + offsets[3], results[3], width, label='test COG', color = cmap_test(1))
    ax.bar(ind + offsets[4], results[4], width, label='train GTR', color = cmap_train(2))
    ax.bar(ind + offsets[5], results[5], width, label='test GTR', color = cmap_test(2))
    ax.bar(ind + offsets[6], results[6], width, label='train MK', color = cmap_train(3))
    ax.bar(ind + offsets[7], results[7], width, label='test MK', color = cmap_test(3))
    ax.bar(ind + offsets[8], results[8], width, label='train MK', color = cmap_train(4))
    ax.bar(ind + offsets[9], results[9], width, label='test MK', color = cmap_test(4))

    ax.set_ylabel('relative llh')
    ax.set_xticks(ind)
    ax.set_xticklabels(range(10))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
      box.width, box.height * 0.9])

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=10)
    ax.legend()
    plt.savefig(os.path.join(plots_dir, ds_name  + ".png"))
    plt.clf()
    plt.close()



def box_plots(results, path):
    models = ["BIN", "MK", "GTR", "COGs", "COG"]
    results_transformed = [[] for _ in range(5)]
    for row in results:
        if row[1] != row[1]:
            continue
        for i in range(1, 6):
            results_transformed[i-1].append(row[i])
    for i in range(5):
        print(models[i], str(statistics.median(results_transformed[i])))
    ax = seaborn.boxplot(data = results_transformed, palette = [cm.to_hex(plt.cm.Set2(num)) for num in range(5)])
    ax.set_xticklabels(models)
    plt.ylabel(r"$e$ (average)")
    plt.savefig(path + "_box.png")
    plt.clf()
    plt.close()



plt.rcParams.update({'font.size': 14})

msa_super_dir = "data/lexibench/character_matrices"
cv_msa_super_dir = "data/bv_cross_validation_data"
raxmlng_super_dir = "data/cross_validation"
plots_super_dir = "data/cross_validation_plots"
tabels_dir = "data/tabels"
if not os.path.isdir(tabels_dir):
    os.makedirs(tabels_dir)
for kappa in range(3, 4): #possible to include other kappa subset sizes here
    random.seed(2)
    diff_headers = ("dataset", "diff_BIN", "diff_MK", "diff_GTR", "diff_COGs", "diff_COG")
    for train_ratio in [60]: # possible to add different split ratios here
        all_diff_res = []
        plots_dir = os.path.join(plots_super_dir, "cv_" + str(train_ratio))
        if not os.path.isdir(plots_dir):
            os.makedirs(plots_dir)
        for ds_name in os.listdir(msa_super_dir):
            msa_dir = os.path.join(msa_super_dir, ds_name)
            target_dir = os.path.join(raxmlng_super_dir, ds_name, "cv_" + str(train_ratio))
            bin_msa_type = "bin_part_" + str(kappa)
            bv_msa_type = "bv_part_" + str(kappa)
            bin_msa_path = os.path.join(msa_dir, bin_msa_type + ".phy")
            bv_msa_path = os.path.join(msa_dir, bv_msa_type + ".phy")
            if not os.path.isfile(bin_msa_path) or not os.path.isfile(bv_msa_path):
                continue
            cv_msa_dir = os.path.join(cv_msa_super_dir, "cv_" + str(train_ratio))
            success = create_samples(kappa, train_ratio / 100, msa_dir, cv_msa_dir)
            if not success:
                continue
            train_raxml_ng(cv_msa_dir, target_dir, kappa)
            test_raxml_ng(cv_msa_dir, target_dir, kappa)
            diff_res = differences_analysis(cv_msa_dir, target_dir, kappa)
            if diff_res[0] == diff_res[0]:
                all_diff_res.append([ds_name] + diff_res)
            else:
                print(diff_res)
            plots(cv_msa_dir, target_dir, kappa, plots_dir, ds_name)
        if len(all_diff_res) > 0:
            box_plots(all_diff_res, os.path.join(plots_dir, str(kappa)))
        print(tabulate(all_diff_res, tablefmt="latex", headers = diff_headers))
        res_df = pd.DataFrame(all_diff_res, columns = diff_headers)
        res_df.to_csv(os.path.join(tabels_dir, "cross_validation_errors_kappa=" + str(kappa) + ".csv"))
