import os
import math
from lingdata import database
from tabulate import tabulate
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error


def substitution_rates(prefix):
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
    return rates




database.read_config("cognate_lingdata_config.json")
df = database.data()
#configurations = [("single", "noforce"), ("multiple", "force")]
configurations = [("multiple", "force")]
msa_types = ["prototype_part_3"]

for rate_mode, force_mode in configurations:
    out_dir = "data/results_cognate_" + rate_mode + "_" + force_mode + "/"
    for msa_type in msa_types:
        gtr_dir = os.path.join("data", "results_prototype", msa_type) #needs to be moved here from repo data-properties
        for i, row in df.iterrows():
            msa_prefix = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
            prefix = os.path.join(out_dir, msa_prefix, msa_type, "COG", "inference")
            rates = substitution_rates(prefix)
            gtr_prefix = os.path.join(gtr_dir, "raxmlng", "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]]), msa_type)
            gtr_rates = substitution_rates(gtr_prefix)
            if len(rates) == 0 or len(gtr_rates) == 0:
                continue
            print(rates)
            print(gtr_rates)
            print(mean_squared_error(rates, gtr_rates))
            
    


