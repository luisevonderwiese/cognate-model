import os
import math
from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import seaborn



def run_inference(msa_path, model, prefix, s = False):
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    args = ""
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
    os.system(command)


def AIC(prefix):
    logpath = prefix + ".raxml.log"
    if not os.path.isfile(logpath):
        return float('nan')
    with open(logpath, "r", encoding = "utf-8") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("AIC"):
            parts = line.split(" / ")
            scores = []
            for part in parts:
                scores.append(float(part.split(" ")[2]))
            return scores[0]
    return float('nan')


def substitution_rates(prefix, x):
    file_path = prefix + ".raxml.log"
    if not os.path.isfile(file_path):
        return []
    with open(file_path, "r", encoding = "utf-8") as logfile:
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
    if x == 4:
        return [rates[0], rates[1]]
    if x == 8:
        return [rates[0], rates[1], rates[14]]
    if x == 16:
        return [rates[0], rates[1], rates[30], rates[76]]
    if x == 32:
        return [rates[0], rates[1], rates[62], rates[172], rates[344]]
    if x == 64:
        return [rates[0], rates[1], rates[126], rates[364], rates[792], rates[1456]]
    print("Illegal x")
    return []

def base_frequencies(prefix, x):
    file_path = prefix + ".raxml.log"
    if not os.path.isfile(file_path):
        return []
    with open(file_path, "r", encoding = "utf-8") as logfile:
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
    if x == 4:
        return [frequencies[0], frequencies[2]]
    if x == 8:
        return [frequencies[0], frequencies[2], frequencies[6]]
    if x == 16:
        return [frequencies[0], frequencies[2], frequencies[6], frequencies[14]]
    if x == 32:
        return [frequencies[0], frequencies[2], frequencies[6], frequencies[14], frequencies[30]]
    if x == 64:
        return [frequencies[0], frequencies[2], frequencies[6], frequencies[14], frequencies[30], frequencies[62]]
    print("Illegal x")
    return []



def run_raxml_ng(msa_dir, target_dir):
    bv_msa_path = os.path.join(msa_dir, "bv_six.phy")
    x = 64
    cog_prefix = os.path.join(target_dir, "bv_six_COG")
    run_inference(bv_msa_path, "COG" + str(x), cog_prefix)
    cogs_prefix = os.path.join(target_dir, "bv_six_COGs")
    run_inference(bv_msa_path, "COG" + str(x), cogs_prefix, s = True)



def AIC_scores(target_dir, kappa):
    results = []
    #bin_msa_type = "bin_part_" + str(kappa)
    bv_msa_type = "bv_part_" + str(kappa)    
    for _, (model, msa_type) in enumerate([("MK", bv_msa_type), ("GTR", bv_msa_type), ("COGs", bv_msa_type), ("COG", bv_msa_type)]):
        prefix = os.path.join(target_dir, msa_type + "_" + model)
        results.append(AIC(prefix))
    return results



def get_all_substitution_rates(raxmlng_super_dir, s = False, kappa = 6):
    r = {}
    for ds_name in os.listdir(raxmlng_super_dir):
        target_dir = os.path.join(raxmlng_super_dir, ds_name)
        bv_msa_type = "bv_part_" + str(kappa)
        msa_dir = os.path.join(msa_super_dir, ds_name)
        bv_msa_path = os.path.join(msa_dir, bv_msa_type + ".phy")
        try:
            with open(bv_msa_path, "r") as msa_file:
                parts = msa_file.readlines()[0].split(" ")
                num_chars = int(parts[2])
                num_langs = int(parts[1])
        except FileNotFoundError:
            continue
        if num_chars < num_langs:
            continue
        prefix = os.path.join(target_dir, bv_msa_type + "_COG")
        if s:
            prefix += "s"
        if s:
            x = -1
        else:
            x = int(math.pow(2, kappa))
        f = substitution_rates(prefix, x)
        #if len(f) == 0:
        #    continue
        r[ds_name] = f
    return r

def get_all_base_frequencies(raxmlng_super_dir, s = False, kappa = 6):
    r = {}
    for ds_name in os.listdir(raxmlng_super_dir):
        target_dir = os.path.join(raxmlng_super_dir, ds_name)
        bv_msa_type = "bv_part_" + str(kappa)
        msa_dir = os.path.join(msa_super_dir, ds_name)        
        bv_msa_path = os.path.join(msa_dir, bv_msa_type + ".phy")
        try:
            with open(bv_msa_path, "r") as msa_file:
                parts = msa_file.readlines()[0].split(" ")
                num_chars = int(parts[2])
                num_langs = int(parts[1])
        except FileNotFoundError:
            continue
        if num_chars < num_langs:
            continue
        prefix = os.path.join(target_dir, bv_msa_type + "_COG")
        if s:
            prefix += "s"
        x = int(math.pow(2, kappa))
        f = base_frequencies(prefix, x)
        #if len(f) == 0:
        #    continue
        r[ds_name]= f
    return r


def violin_plots(results, path):
    models = ["MK", "GTR", "COGs", "COG"]
    results_transformed = [[] for _ in range(4)]
    for row in results:
        if row[1] != row[1]:
            continue
        for i in range(1, 5):
            results_transformed[i-1].append(row[i])
    ax = seaborn.violinplot(data = results_transformed, palette = [cm.to_hex(plt.cm.Set2(num)) for num in range(1, 5)])
    ax.set_xticklabels(models)
    plt.xlabel("model")
    plt.ylabel("AIC")
    plt.savefig(path)
    plt.clf()
    plt.close()



def rates_stacked_plot(all_rates, path, plot_type):
    all_rates = dict((ds, rates) for ds, rates in all_rates.items() if rates != [])
    if plot_type == "bf":
        for ds, rates in all_rates.items():
            all_rates[ds] = [0] + rates  #because of colors
    max_num = max([len(rates) for _, rates in all_rates.items()])
    _, ax = plt.subplots(figsize=(15, 10))
    y_old = [0 for _ in all_rates.keys()]
    all_rates = {k: v for k, v in sorted(list(all_rates.items()))}
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
        for __, rates in all_rates.items():
            if len(rates) > num:
                y_new.append(rates[num] / sum(rates))
            else:
                y_new.append(0)
        if plot_type == "bf" and num == 0:
            ax.bar(all_rates.keys(), y_new, bottom=y_old)
        else:
            ax.bar(all_rates.keys(), y_new, bottom=y_old, label = label_list[num])
        for i in range(len(all_rates)):
            y_old[i] = y_old[i] + y_new[i]
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
          fancybox=True, shadow=True, ncol=max_num)
    plt.xticks(rotation=30, ha='right')
    plt.xlabel("Datasets")
    if plot_type == "bf":
        plt.ylabel("Base Frequencies (relative)")
    if plot_type == "sr":
        plt.ylabel("Substitution Rates (relative)")
    plt.savefig(path)
    plt.clf()
    plt.close()


def raxml_ng():
    for ds_name in os.listdir(msa_super_dir):
        msa_dir = os.path.join(msa_super_dir, ds_name)
        target_dir = os.path.join(raxmlng_super_dir, ds_name)
        bv_msa_path = os.path.join(msa_dir, "bv_six.phy")
        if not os.path.isfile(bv_msa_path):
            continue
        run_raxml_ng(msa_dir, target_dir)

def AIC_analysis(kappa):
    AIC_res = []
    for ds_name in os.listdir(msa_super_dir):
        msa_dir = os.path.join(msa_super_dir, ds_name)
        target_dir = os.path.join(raxmlng_super_dir, ds_name)
        bin_msa_type = "bin_part_" + str(kappa)
        bv_msa_type = "bv_part_" + str(kappa)
        bin_msa_path = os.path.join(msa_dir, bin_msa_type + ".phy")
        bv_msa_path = os.path.join(msa_dir, bv_msa_type + ".phy")
        if not os.path.isfile(bin_msa_path) or not os.path.isfile(bv_msa_path):
            continue
        with open(bv_msa_path, "r") as msa_file:
            parts = msa_file.readlines()[0].split(" ")
            num_chars = int(parts[2])
            num_langs = int(parts[1])
        if num_chars < num_langs:
            continue
        AIC_res.append([ds_name] + AIC_scores(target_dir, kappa))
    violin_plots(AIC_res, os.path.join(plots_super_dir, "AIC_" + str(kappa) + ".png"))
    print(tabulate(AIC_res, tablefmt="pipe", headers = ["dataset", "MK", "GTR", "COGs", "COG"]))



msa_super_dir = "data/lexibench_six/character_matrices"
plots_super_dir = os.path.join("data", "plots_six")
raxmlng_super_dir = os.path.join("data","inferences_six")
if not os.path.isdir(plots_super_dir):
    os.makedirs(plots_super_dir)

raxml_ng()
#AIC_analysis()

rates_stacked_plot(get_all_substitution_rates(raxmlng_super_dir), os.path.join(plots_super_dir, "substitution_rates.png"), "sr")
rates_stacked_plot(get_all_substitution_rates(raxmlng_super_dir, True), os.path.join(plots_super_dir, "substitution_rates_s.png"), "sr")
rates_stacked_plot(get_all_base_frequencies(raxmlng_super_dir), os.path.join(plots_super_dir, "base_frequencies.png"), "bf")
