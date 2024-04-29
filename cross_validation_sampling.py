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

def split_indices(num_sites, num_samples = 10, ratio = 0.6):
    num_sites_train = math.ceil(num_sites * ratio)
    res = []
    for i in range(num_samples):
        l = [_ for _ in range(num_sites)]
        random.shuffle(l)
        train_indices = l[:num_sites_train]
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
    command = "./bin/raxml-ng-multiple-force --evaluate "
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


def create_samples(kappa, msa_dir):
    bin_msa_type = "bin_part_" + str(kappa)
    prototype_msa_type = "prototype_part_" + str(kappa)
    bin_msa_path = os.path.join(msa_dir, bin_msa_type + ".phy")
    prototype_msa_path = os.path.join(msa_dir, prototype_msa_type + ".phy")
    with open(prototype_msa_path, "r") as msa_file:
        num_sites = int(msa_file.readlines()[0].split(" ")[2])
    with open(bin_msa_path, "r") as msa_file:
        num_sites_bin = int(msa_file.readlines()[0].split(" ")[2])
    assert(num_sites_bin == kappa * num_sites)
    try:
        bin_align = AlignIO.read(bin_msa_path, "phylip-relaxed")
    except:
        print(msa_dir, "Failed")
        return False
    try:
        prototype_align = AlignIO.read(bin_msa_path, "phylip-relaxed")
    except:
        print(msa_dir, "Failed")
        return False
    indices_list = split_indices(num_sites)
    for (t, train_indices) in enumerate(indices_list):
        bin_train_align = empty_align(bin_align)
        bin_test_align = empty_align(bin_align)
        prototype_train_align = empty_align(prototype_align)
        prototype_test_align = empty_align(prototype_align)
        for s in range(num_sites):
            if s in train_indices :
                prototype_train_align = concat_align(prototype_train_align, prototype_align[:, s:s+1])
                bin_train_align = concat_align(bin_train_align, bin_align[:, s*kappa : (s+1) * kappa])
            else:
                prototype_test_align = concat_align(prototype_test_align, prototype_align[:, s:s+1])
                bin_test_align = concat_align(bin_test_align, bin_align[:, s*kappa : (s+1) * kappa])
        with open(os.path.join(msa_dir, bin_msa_type + "_cv_train_" + str(t) + ".phy"),"w+") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(bin_train_align)
        with open(os.path.join(msa_dir, bin_msa_type + "_cv_test_" + str(t) + ".phy"),"w+") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(bin_test_align)
        with open(os.path.join(msa_dir, prototype_msa_type + "_cv_train_" + str(t) + ".phy"),"w+") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(prototype_train_align)
        with open(os.path.join(msa_dir, prototype_msa_type + "_cv_test_" + str(t) + ".phy"),"w+") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(prototype_test_align)
    print(msa_dir, "done")
    return True


def train_raxml_ng(msa_dir, target_dir, kappa):
    bin_msa_type = "bin_part_" + str(kappa)
    prototype_msa_type = "prototype_part_" + str(kappa)
    x = int(math.pow(2, kappa))
    for t in range(10):
        bin_msa_path = os.path.join(msa_dir, bin_msa_type + "_cv_train_" + str(t) + ".phy")
        bin_prefix = os.path.join(target_dir, bin_msa_type + "_cv_train_" + str(t) + "_BIN")
        run_inference(bin_msa_path, "BIN", bin_prefix)
        prototype_msa_path = os.path.join(msa_dir, bin_msa_type + "_cv_train_" + str(t) + ".phy")
        prototype_prefix = os.path.join(target_dir, prototype_msa_type + "_cv_train_" + str(t) + "_COG")
        run_inference(prototype_msa_path, "COG" + str(x), prototype_prefix)
        gtr_prefix = os.path.join(target_dir, prototype_msa_type + "_cv_train_" + str(t) + "_GTR")
        run_inference(prototype_msa_path, "MULTI" + str(x - 1) + "_GTR", gtr_prefix)


def test_raxml_ng(msa_dir, target_dir, kappa):
    bin_msa_type = "bin_part_" + str(kappa)
    prototype_msa_type = "prototype_part_" + str(kappa)
    x = int(math.pow(2, kappa))
    for t in range(10):
        bin_msa_path = os.path.join(msa_dir, bin_msa_type + "_cv_test_" + str(t) + ".phy")
        bin_prefix = os.path.join(target_dir, bin_msa_type + "_cv_train_" + str(t) + "_BIN")
        bin_test_prefix = os.path.join(target_dir, bin_msa_type + "_cv_test_" + str(t) + "_BIN")
        run_evaluate(bin_msa_path, bin_test_prefix, bin_prefix)
        prototype_msa_path = os.path.join(msa_dir, bin_msa_type + "_cv_test_" + str(t) + ".phy")
        prototype_prefix = os.path.join(target_dir, prototype_msa_type + "_cv_train_" + str(t) + "_COG")
        prototype_test_prefix = os.path.join(target_dir, prototype_msa_type + "_cv_test_" + str(t) + "_COG")
        run_evaluate(bin_msa_path, prototype_test_prefix, prototype_prefix)
        gtr_prefix = os.path.join(target_dir, prototype_msa_type + "_cv_train_" + str(t) + "_GTR")
        gtr_test_prefix = os.path.join(target_dir, prototype_msa_type + "_cv_test_" + str(t) + "_GTR")
        run_evaluate(bin_msa_path, gtr_test_prefix, gtr_prefix)


def analysis(target_dir, kappa):
    bin_msa_type = "bin_part_" + str(kappa)
    prototype_msa_type = "prototype_part_" + str(kappa)
    results = [[] for _ in range(6)]
    for t in range(10):
        for m, (model, msa_type) in enumerate([("BIN", bin_msa_type), ("COG", prototype_msa_type), ("GTR", prototype_msa_type)]):
            results[m * 2].append(final_llh(os.path.join(target_dir, msa_type + "_cv_train_" + str(t) + "_" + model)))
            results[m * 2 + 1].append(final_llh(os.path.join(target_dir, msa_type + "_cv_test_" + str(t) + "_" +  model )))
    return [sum(el) / len(el) for el in results]


def differences_analysis(target_dir, kappa):
    bin_msa_type = "bin_part_" + str(kappa)
    prototype_msa_type = "prototype_part_" + str(kappa)
    results = [[] for _ in range(3)]
    for t in range(10):
        for m, (model, msa_type) in enumerate([("BIN", bin_msa_type), ("COG", prototype_msa_type), ("GTR", prototype_msa_type)]):
            results[m].append(final_llh(os.path.join(target_dir, msa_type + "_cv_train_" + str(t) + "_" + model)) - \
            final_llh(os.path.join(target_dir, msa_type + "_cv_test_" + str(t) + "_" + model)))
    return [sum(el) / len(el) for el in results]


def plots(target_dir, kappa, plots_super_dir, ds_name):
    bin_msa_type = "bin_part_" + str(kappa)
    prototype_msa_type = "prototype_part_" + str(kappa)
    ind = np.arange(10)
    width = 0.35
    for (model, msa_type) in [("BIN", bin_msa_type), ("COG", prototype_msa_type), ("GTR", prototype_msa_type)]:
        plots_dir = os.path.join(plots_super_dir, model)
        if not os.path.isdir(plots_dir):
            os.makedirs(plots_dir)
        train_llhs = [final_llh(os.path.join(target_dir, msa_type + "_cv_train_" + str(t) + "_" + model)) for t in range(10)]
        test_llhs = [final_llh(os.path.join(target_dir, msa_type + "_cv_test_" + str(t) + "_" + model)) for t in range(10)]
        fig, ax = plt.subplots()
        rects1 = ax.bar(ind - width/2, train_llhs, width, label='train')
        rects2 = ax.bar(ind + width/2, test_llhs, width, label='test')
        ax.set_ylabel('final llh')
        ax.set_xticks(ind)
        ax.set_xticklabels(range(10))
        ax.legend()
        plt.savefig(os.path.join(plots_dir, ds_name  + ".png"))
        plt.clf()
        plt.close()






msa_super_dir = "data/lingdata_cognate/msa"
raxmlng_super_dir = "data/cross_validation"
plots_super_dir = "data/cross_validation_plots"
kappa = 3
random.seed(2)
all_res = []
all_diff_res = []
headers = ("dataset", "train_BIN", "test_BIN", "train_COG", "test_COG", "train_GTR", "test_GTR")
diff_headers = ("dataset", "diff_BIN", "diff_COG", "diff_GTR")
for ds_name in os.listdir(msa_super_dir):
    msa_dir = os.path.join(msa_super_dir, ds_name)
    target_dir = os.path.join(raxmlng_super_dir, ds_name)
    bin_msa_type = "bin_part_" + str(kappa)
    prototype_msa_type = "prototype_part_" + str(kappa)
    bin_msa_path = os.path.join(msa_dir, bin_msa_type + ".phy")
    prototype_msa_path = os.path.join(msa_dir, prototype_msa_type + ".phy")
    if not os.path.isfile(bin_msa_path) or not os.path.isfile(prototype_msa_path):
        continue
    success = create_samples(kappa, msa_dir)
    if not success:
        continue
    #train_raxml_ng(msa_dir, target_dir, kappa)
    #test_raxml_ng(msa_dir, target_dir, kappa)
    all_res.append([ds_name] + analysis(target_dir, kappa))
    all_diff_res.append([ds_name] + differences_analysis(target_dir, kappa))
    plots(target_dir, kappa, plots_super_dir, ds_name)
print(tabulate(all_res, tablefmt="pipe", headers = headers))
print(tabulate(all_diff_res, tablefmt="pipe", headers = diff_headers))
