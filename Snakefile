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
multi_part_models = ["part2MK", "partxMK"]


def add_for_raxml(row, msa_prefix, msa_type, model, partitioned = False):
    run_prefix = os.path.join(msa_prefix, msa_type, model)
    prefixes.append(run_prefix)
    msa_path_dict[run_prefix] = row["msa_path"][msa_type]
    if partitioned:
        model_dict[run_prefix] = row["partition_paths"][model]
    elif model.startswith("MULTI_"):
        model_dict[run_prefix] = row["MULTIx_" + model]
    else:
        model_dict[run_prefix] = model
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
    add_for_raxml(row, msa_prefix, "bin", "BIN+G")
    add_for_raxml(row, msa_prefix, "catg_bin", "BIN")
    add_for_raxml(row, msa_prefix, "catg_bin", "BIN+G")

    if row["msa_paths"]["multi"] != "":
        add_for_pythia(row, msa_prefix, "multi")
        for model in multi_models:
           add_for_raxml(row, msa_prefix, "multi", "MULTIx_" + model)
           add_for_raxml(row, msa_prefix, "multi", "MULTIx_" + model + "+G")
        for model in multi_part_models:
            add_for_raxml(row, msa_prefix, "multi", model, partitioned = True)

    if row["msa_paths"]["catg_multi"] != "":
        for model in multi_models:
           add_for_raxml(row, msa_prefix, "catg_multi", "MULTIx_" + model)
           add_for_raxml(row, msa_prefix, "catg_multi", "MULTIx_" + model + "+G")

    if row["msa_paths"]["ambig"] != "":
        add_for_pythia(row, msa_prefix, "ambig")
        for model in multi_models:
           add_for_raxml(row, msa_prefix, "ambig", "MULTIx_" + model)
           add_for_raxml(row, msa_prefix, "ambig", "MULTIx_" + model + "+G")
        for model in multi_part_models:
            add_for_raxml(row, msa_prefix, "ambig", model, partitioned = True)


seed = 2

outdir = config["outdir"]
results_dir = outdir + "{run_prefix}/"

rule all:
    input:
        expand(f"{results_dir}results.json", run_prefix=prefixes),
        expand(f"{results_dir}pythia.difficulty", run_prefix=pythia_prefixes)


include: "rules.smk"
