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



for i, row in df.iterrows():
    msa_prefix = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
    run_prefix = os.path.join(msa_prefix, "prototype", "GTR")
    prefixes.append(run_prefix)
    msa_path_dict[run_prefix] = row["msa_paths"]["prototype"]
    model_dict[run_prefix] = row["MULTIx_GTR_prototype"]
    prob_msa_dict[run_prefix] = "off"

    pythia_prefix = os.path.join(msa_prefix, "prototype")
    pythia_prefixes.append(pythia_prefix)
    pythia_msa_path_dict[pythia_prefix] = row["msa_paths"]["prototype"]




seed = 2

outdir = config["outdir"]
results_dir = outdir + "{run_prefix}/"

rule all:
    input:
        expand(f"{results_dir}results.json", run_prefix=prefixes),
        expand(f"{results_dir}pythia.difficulty", run_prefix=pythia_prefixes)


include: "rules.smk"
