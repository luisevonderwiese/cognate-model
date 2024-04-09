import os
import json
import math
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

def gq_distance_analysis(df):
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
            if not os.path.isdir(msa_path):
                continue
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
        results.append(r)

    print(tabulate(results, tablefmt="pipe", floatfmt=".3f", headers = ["ds_id", "GQ BIN", "GQ COG"]))
    plt.axline([0, 0], slope=1, color = 'lightgray', linewidth = 1, linestyle = "--")
    plt.scatter(gq_distances["BIN"], gq_distances["COG"], s=10)
    plt.xlabel('BIN')
    plt.ylabel('COG')
    plt.savefig(os.path.join(plots_dir, "scatter_cognate.png"))
    plt.clf()
    plt.close()

def get_all_substitution_rates(df, msa_type):
    r = []
    for i, row in df.iterrows():
        msa_prefix = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
        prefix = os.path.join(out_dir, msa_prefix, msa_type, "COG", "inference")
        if rate_mode == "single":
            r.append(substitution_rates(prefix, -1))
        elif msa_type.startswith("prototype_part"):
            x = int(math.pow(2, int(msa_type.split("_")[-1])))
            r.append(substitution_rates(prefix, x))
        else:
            r.append(substitution_rates(prefix, row["max_values_prototype"])) 
    return r

def get_all_base_frequencies(df, msa_type):
    r = []
    for i, row in df.iterrows():
        msa_prefix = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
        prefix = os.path.join(out_dir, msa_prefix, msa_type, "COG", "inference")
        if msa_type.startswith("prototype_part"):
            x = int(math.pow(2, int(msa_type.split("_")[-1])))
            r.append(base_frequencies(prefix, x))
        else:
            r.append(base_frequencies(prefix, row["max_values_prototype"]))
    return r

def rates_stacked_plot(all_rates, file_name, plot_type):
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
    plt.savefig(os.path.join(plots_dir, file_name))
    plt.clf()
    plt.close()

database.read_config("cognate_lingdata_config.json")
df = database.data()

configurations = [("single", "noforce"), ("multiple", "force")]
for rate_mode, force_mode in configurations:
    out_dir = "data/results_cognate_" + rate_mode + "_" + force_mode + "/"
    plots_dir = "data/plots_cognate_" + rate_mode + "_" + force_mode + "/"
    if not os.path.isdir(plots_dir):
        os.makedirs(plots_dir)

    #gq_distance_analysis(df)

    #msa_types = ["prototype"] + ["prototype_part_" + str(i) for i in range(2, 7)]
    msa_types = ["prototype_part_" + str(i) for i in range(3, 6)] 
    for msa_type in msa_types:
        rates_stacked_plot(get_all_substitution_rates(df, msa_type), "substitution_rates_" + msa_type + ".png", "sr")
        rates_stacked_plot(get_all_base_frequencies(df, msa_type), "base_frequencies_" + msa_type + ".png", "bf")
