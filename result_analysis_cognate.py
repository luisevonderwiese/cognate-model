import os
import json
from lingdata import database
from ete3 import Tree
from tabulate import tabulate
import matplotlib.pyplot as plt


def rf_distance(t1, t2):
    if t1 is None or t2 is None:
        return float('nan')
    if t1 != t1 or t2 != t2:
        return float("nan")
    rf, max_rf, common_leaves, parts_t1, parts_t2,discard_t1, discart_t2 = t1.robinson_foulds(t2, unrooted_trees = True)
    if max_rf == 0:
        return float('nan')
    return rf/max_rf

def gq_distance(tree_name1, tree_name2):
    if tree_name1 is None or tree_name2 is None:
        return float('nan')
    if tree_name1 != tree_name1 or tree_name2 != tree_name2:
        return float("nan")
    os.system("./bin/qdist " + tree_name1 + " " + tree_name2 + " >out.txt")
    lines = open("out.txt").readlines()
    if len(lines) < 2: #error occurred
        return float('nan')
    res_q = float(lines[1].split("\t")[-3])
    qdist = 1 - res_q
    os.remove("out.txt")
    return qdist

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
        return rates
    if x == -1: #code for single lamda rate
        if len(rates) <= 2:
            return rates
        return [rates[2], rates[0]]
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
        return []
    final_f = []
    cursor = 1
    while cursor <= x:
        final_f.append(frequencies[cursor-1])
        cursor *= 2
    return final_f




out_dir = "data/results_cognate/"
plots_dir = "data/plots_cognate"
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)
database.read_config("cognate_lingdata_config.json")
df = database.data()
results = []
gq_distances = {"BIN" : [], "COG": []}
for i, row in df.iterrows():
    msa_prefix = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
    glottolog_tree_path = row["glottolog_tree_path"]
    r = []
    r.append(row["ds_id"])
    ds_path = os.path.join(out_dir, msa_prefix)
    for msa_type in ["bin", "prototype"]:
        msa_path = os.path.join(ds_path, msa_type)
        for run_name in os.listdir(msa_path):
            if run_name.startswith("pythia"):
                    continue
            best_tree_path = os.path.join(msa_path, run_name, "inference.raxml.bestTree")
            if glottolog_tree_path == glottolog_tree_path and os.path.isfile(glottolog_tree_path)  and os.path.isfile(best_tree_path):
                d = gq_distance(best_tree_path, glottolog_tree_path)
                r.append(d)
                gq_distances[run_name].append(d)
            else:
                r.append(float("nan"))
                gq_distances[run_name].append(float("nan"))
            if run_name == "COG":
                prefix = os.path.join(msa_path, run_name, "inference")
                r.append(substitution_rates(prefix, row["max_values_prototype"])) #-1 for single
                r.append(base_frequencies(prefix, row["max_values_prototype"]))

    results.append(r)

print(tabulate(results, tablefmt="pipe", floatfmt=".3f", headers = ["ds_id", "GQ BIN", "GQ COG", "subtitution rates", "base_frequencies"]))


plt.axline([0, 0], slope=1, color = 'lightgray', linewidth = 1, linestyle = "--")
plt.scatter(gq_distances["BIN"], gq_distances["COG"], s=10)
plt.xlabel('BIN')
plt.ylabel('COG')
plt.savefig(os.path.join(plots_dir, "scatter_cognate.png")) 
plt.clf()

all_rates = [r[3] for r in results]
max_num = max([len(rates) for rates in all_rates])
fig,ax = plt.subplots(figsize=(40, 30))
x = range(len(all_rates))
y_old = [0 for el in x]
for num in range(max_num):
    y_new = []
    for rates in all_rates:
        if len(rates) > num:
            y_new.append(rates[num] /sum(rates))
        else:
            y_new.append(0)
    ax.bar(x, y_new, bottom=y_old, label = str(num))
    for i in x:
        y_old[i] = y_old[i] + y_new[i]
ax.legend()
plt.savefig(os.path.join(plots_dir, "stacked_substitution_rates.png"))
plt.clf()

all_rates = [r[4] for r in results]
max_num = max([len(rates) for rates in all_rates])
fig,ax = plt.subplots(figsize=(40, 30))
x = range(len(all_rates))
y_old = [0 for el in x]
for num in range(max_num):
    y_new = []
    for rates in all_rates:
        if len(rates) > num:
            y_new.append(rates[num] / sum(rates))
        else:
            y_new.append(0)
    ax.bar(x, y_new, bottom=y_old, label = str(num))
    for i in x:
        y_old[i] = y_old[i] + y_new[i]
ax.legend()
plt.savefig(os.path.join(plots_dir, "stacked_base_frequencies.png"))
plt.clf()

