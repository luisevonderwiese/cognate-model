import os
import sys
import pandas as pd
from lingdata import database


configfile: "config.yaml"

raxmlng_command = config["software"]["raxml-ng"]["command"]
pythia_command = config["software"]["pythia"]["command"]
pythia_predictor = config["software"]["pythia"]["predictor"]

database.read_config(config["lingdata_config_path"])
df = database.data()

prefixes = []
msa_path_dict = {}
model_dict = {}
prob_msa_dict = {}

pythia_prefixes = []
pythia_msa_path_dict = {}

multi_models = ["MK"]


def add_for_raxml(row, msa_prefix, msa_type, model, gamma = False, partition_mode = ""):
    model_name = model
    if msa_type == "ambig":
        model_name += "+M"
    if gamma:
        model_name += "+G"
    if paritition_mode != "":
        assert(partition_mode in ["x", "2"])
        model_name += "_" + paritition_mode
    run_prefix = os.path.join(msa_prefix, msa_type, model_name)
    prefixes.append(run_prefix)
    msa_path_dict[run_prefix] = row["msa_paths"][msa_type]
    if partition_mode != "":
        model_dict[run_prefix] = row["partition_paths"][msa_type + "_" + model_name]
    else:
        if model in ["MK", "GTR"]:
            if msa_type == "ambig":
                model_dict[run_prefix] = row["MULTIx_" + model + "+M"]
            else:
                assert(msa_type in ["multi", "catg_multi"])
                model_dict[run_prefix] = row["MULTIx_" + model]
        else:
            assert(model == "BIN")
            assert(msa_type in ["bin", "catg_bin"])
            model_dict[run_prefix] = model
        if gamma:
            model_dict[run_prefix] += "+G"
    if msa_type.startswith("catg_"):
        prob_msa_dict[run_prefix] = "on"
    else:
        prob_msa_dict[run_prefix] = "off"

def add_for_pythia(row, msa_prefix, msa_type):
    pythia_prefix = os.path.join(msa_prefix, msa_type)
    pythia_prefixes.append(pythia_prefix)
    pythia_msa_path_dict[pythia_prefix] = row["msa_paths"][msa_type]


for i, row in df.iterrows():
    msa_prefix = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])

    add_for_pythia(row, msa_prefix, "bin")

    add_for_raxml(row, msa_prefix, "bin", "BIN")
    add_for_raxml(row, msa_prefix, "bin", "BIN", True)
    add_for_raxml(row, msa_prefix, "catg_bin", "BIN")
    add_for_raxml(row, msa_prefix, "catg_bin", "BIN", True)

    if row["msa_paths"]["multi"] != "":
        add_for_pythia(row, msa_prefix, "multi")
        for model in multi_models:
           add_for_raxml(row, msa_prefix, "multi", model)
           add_for_raxml(row, msa_prefix, "multi", model, True)
           add_for_raxml(row, msa_prefix, "multi", model, False, x)
           add_for_raxml(row, msa_prefix, "multi", model, True, x)
           add_for_raxml(row, msa_prefix, "multi", model, False, 2)
           add_for_raxml(row, msa_prefix, "multi", model, True, 2)

    if row["msa_paths"]["catg_multi"] != "":
        for model in multi_models:
           add_for_raxml(row, msa_prefix, "catg_multi", model)
           add_for_raxml(row, msa_prefix, "catg_multi", model, True)

    if row["msa_paths"]["ambig"] != "":
        add_for_pythia(row, msa_prefix, "ambig")
        for model in multi_models:
           add_for_raxml(row, msa_prefix, "ambig", model)
           add_for_raxml(row, msa_prefix, "ambig", model, True)


seed = 2

outdir = config["outdir"]
results_dir = outdir + "{run_prefix}/"

rule all:
    input:
        expand(f"{results_dir}results.json", run_prefix=prefixes),
        expand(f"{results_dir}pythia.difficulty", run_prefix=pythia_prefixes)


include: "rules.smk"
