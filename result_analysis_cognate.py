import os
import json
from lingdata import database
from ete3 import Tree


def rf_distance(t1, t2):
    if t1 is None or t2 is None:
        return float('nan')
    if t1 != t1 or t2 != t2:
        return float("nan")
    rf, max_rf, common_leaves, parts_t1, parts_t2,discard_t1, discart_t2 = t1.robinson_foulds(t2, unrooted_trees = True)
    if max_rf == 0:
        return float('nan')
    return rf/max_rf

def substitution_rates(prefix, x):
    with open(prefix + ".raxml.log", "r") as logfile:
        lines = logfile.readlines()
    rates = []
    for line in lines:
        if line.startswith("   Substitution rates"):
            parts = line.split(": ")[1].split(" ")[:-1]
            rates = [float(part) for part in parts]
            break
    if rates == []:
        return rates
    if x == 2:
        return rates
    elif x == 4:
        return [rates[2], rates[0], rates[4]]
    elif x == 8:
        return [rates[2], rates[0], rates[8], rates[21]]
    elif x == 16:
        return [rates[2], rates[0], rates[16], rates[45], rates[91]]
    elif x == 32:
        return [rates[2], rates[0], rates[32], rates[93], rates[203], rates[375]]
    elif x == 64:
        return [rates[2], rates[0], rates[64], rates[189], rates[427], rates[855], rates[1519]]
    else:
        print("Illegal x")
        return []

out_dir = "data/results_cognate/"
database.read_config("cognate_lingdata_config.json")
df = database.data()
big_results_dict = {}
for i, row in df.iterrows():
    msa_prefix = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
    glottolog_tree_path = row["glottolog_tree_path"]
    small_results_dict = {}
    ds_path = os.path.join(out_dir, msa_prefix)
    for msa_type in os.listdir(ds_path):
        msa_path = os.path.join(ds_path, msa_type)
        for run_name in os.listdir(msa_path):
            best_tree_path = os.path.join(msa_path, run_name, "inference.raxml.bestTree")
            if os.path.isfile(glottolog_tree_path)  and os.path.isfile(best_tree_path):
                small_results_dict["RF_" + msa_type + "_" + run_name] = rf_distance(Tree(best_tree_path), Tree(glottolog_tree_path))
            if run_name == "COG":
                prefix = os.path.join(msa_path, run_name, "inference")
                small_results_dict["substitution_rates"] = substitution_rates(prefix, row["max_values_prototype"])

    big_results_dict[msa_prefix] = small_results_dict

print(big_results_dict)
