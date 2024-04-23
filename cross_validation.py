import os
import json
import math
from lingdata import database
from ete3 import Tree
from tabulate import tabulate
import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter


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

def run_check(msa_path, ref_prefix, args = ""):
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    if not os.path.isfile(ref_prefix + ".raxml.bestModel"):
        return
    with open(ref_prefix + ".raxml.bestModel", "r") as model_file:
        model =  model_file.readlines()[0].split(",")[0]
    command = "./bin/raxml-ng-multiple-force --check "
    command += " --msa " + msa_path
    command += " --model " + ref_prefix + ".raxml.bestModel"
    command += " --prefix check --redo --threads auto --seed 2"
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



def create_cv_data(df, msa_type):
    kappa = int(msa_type.split("_")[-1])
    for i, row in df.iterrows():
        msa_path = row["msa_paths"][msa_type]
        if not os.path.isfile(msa_path):
            continue
        msa_prefix = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
        align = AlignIO.read(msa_path, "phylip-relaxed")
        num_sites = align.get_alignment_length()
        half = int(num_sites/2)
        train_align = align[:, :half]
        train_dir = os.path.join(cv_data_dir, msa_prefix, "train")
        if not os.path.isdir(train_dir):
            os.makedirs(train_dir)
        with open(os.path.join(train_dir, msa_type + ".phy"),"w+") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(train_align)
        test_align = align[:, half:]
        test_dir = os.path.join(cv_data_dir, msa_prefix, "test")
        if not os.path.isdir(test_dir):
            os.makedirs(test_dir)
        with open(os.path.join(test_dir, msa_type + ".phy"),"w+") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(test_align)



def train_raxml_ng(df, msa_type):
    kappa = int(msa_type.split("_")[-1])
    x = int(math.pow(2, kappa))
    model = "COG" + str(x)
    model = "MULTI" + str(x) + "_GTR"
    for i, row in df.iterrows():
        msa_prefix = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
        train_dir = os.path.join(cv_data_dir, msa_prefix, "train")
        train_msa_path = os.path.join(train_dir, msa_type + ".phy")
        if not os.path.isfile(train_msa_path):
            continue
        prefix = os.path.join(cv_results_dir, msa_prefix, "train", msa_type)
        run_inference(train_msa_path, model, prefix)


def test_raxml_ng(df, msa_type):
    for i, row in df.iterrows():
        msa_prefix = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
        train_dir = os.path.join(cv_data_dir, msa_prefix, "train")
        train_msa_path = os.path.join(train_dir, msa_type + ".phy")
        if not os.path.isfile(train_msa_path):
            continue
        ref_prefix = os.path.join(cv_results_dir, msa_prefix, "train", msa_type)
        prefix = os.path.join(cv_results_dir, msa_prefix, "test", msa_type)
        test_dir = os.path.join(cv_data_dir, msa_prefix, "test")
        msa_path = os.path.join(test_dir, msa_type + ".phy")
        #run_check(msa_path, ref_prefix)
        run_evaluate(msa_path, prefix, ref_prefix)


def validate(df, msa_type):
    r = []
    headers = ["msa", "train_llh", "test_llh"]
    for i, row in df.iterrows():
        msa_prefix = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
        train_dir = os.path.join(cv_data_dir, msa_prefix, "train")
        train_msa_path = os.path.join(train_dir, msa_type + ".phy")
        if not os.path.isfile(train_msa_path):
            continue
        ref_prefix = os.path.join(cv_results_dir, msa_prefix, "train", msa_type)
        prefix = os.path.join(cv_results_dir, msa_prefix, "test", msa_type)
        train_llh = final_llh(ref_prefix)
        test_llh = final_llh(prefix)
        r.append([msa_prefix, train_llh, test_llh])
    print(tabulate(r, tablefmt="pipe", headers = headers))







database.read_config("cognate_lingdata_config.json")
df = database.data()
cv_data_dir = "data/lingdata_cognate/msa_cv"
cv_results_dir = "data/results_cv_gtr"
msa_type = "prototype_part_3"
#create_cv_data(df, msa_type)
train_raxml_ng(df, msa_type)
test_raxml_ng(df, msa_type)
validate(df, msa_type)
