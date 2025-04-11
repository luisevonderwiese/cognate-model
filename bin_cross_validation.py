import os
from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import seaborn
import statistics

def run_inference(msa_path, model, prefix, args = ""):
    if not os.path.isfile(msa_path):
        print("MSA " + msa_path + " does not exist")
        return
    prefix_dir = "/".join(prefix.split("/")[:-1])
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)
    if not os.path.isfile(prefix + ".raxml.bestTree"):
        args = args + " --redo"
    command = "./bin/raxml-ng-COG"
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
    with open(ref_prefix + ".raxml.bestModel", "r", encoding = "utf-8") as model_file:
        model =  model_file.readlines()[0].split(",")[0]
    command = "./bin/raxml-ng-COG --evaluate "
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
    with open(prefix + ".raxml.log", "r", encoding = "utf-8") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("Final LogLikelihood: "):
            return float(line.split(": ")[1])
    return float('nan')

def relative_llh(msa_path, prefix):
    with open(msa_path, "r", encoding = "utf-8") as msa_file:
        num_sites = int(msa_file.readlines()[0].split(" ")[2])
    return final_llh(prefix) / num_sites
    #return final_llh(prefix)


def train_raxml_ng(msa_dir, target_dir):
    for t in range(10):
        bin_msa_path = os.path.join(msa_dir, "train", "bin" + str(t) + ".phy")
        bin_prefix = os.path.join(target_dir, "train", "bin" + str(t), "BIN")
        run_inference(bin_msa_path, "BIN", bin_prefix)



def test_raxml_ng(msa_dir, target_dir):
    for t in range(10):
        bin_msa_path = os.path.join(msa_dir, "test", "bin" + str(t) + ".phy")
        bin_prefix = os.path.join(target_dir, "train", "bin" + str(t), "BIN")
        bin_test_prefix = os.path.join(target_dir, "test", "bin" + str(t), "BIN")
        run_evaluate(bin_msa_path, bin_test_prefix, bin_prefix)



def analysis(msa_dir, target_dir):
    results = [[] for _ in range(2)]
    for t in range(10):
        for m, (model, msa_type) in enumerate([("BIN", "bin")]):
            train_msa_path = os.path.join(msa_dir, "train", msa_type + str(t) + ".phy")
            train_prefix = os.path.join(target_dir, "train", msa_type +  str(t), model)
            results[m * 2].append(relative_llh(train_msa_path, train_prefix))
            test_msa_path = os.path.join(msa_dir, "test", msa_type + str(t) + ".phy")
            test_prefix = os.path.join(target_dir, "test", msa_type +  str(t), model)
            results[m * 2 + 1].append(relative_llh(test_msa_path, test_prefix))
    return [sum(el) / len(el) for el in results]


def differences_analysis(msa_dir, target_dir):
    results = [[] for _ in range(1)]
    for t in range(10):
        for m, (model, msa_type) in enumerate([("BIN", "bin")]):
            train_msa_path = os.path.join(msa_dir, "train", msa_type + str(t) + ".phy")
            train_prefix = os.path.join(target_dir, "train", msa_type +  str(t), model)
            rel_train_llh = relative_llh(train_msa_path, train_prefix)
            test_msa_path = os.path.join(msa_dir, "test", msa_type + str(t) + ".phy")
            test_prefix = os.path.join(target_dir, "test", msa_type +  str(t), model)
            rel_test_llh = relative_llh(test_msa_path, test_prefix)
            results[m].append((rel_test_llh - rel_train_llh) / rel_train_llh)
            #results[m].append(rel_train_llh - rel_test_llh)
    return [sum(el) / len(el) for el in results]





def box_plots(results, path):
    if not os.path.isdir(path):
        os.makedirs(path)
    models = ["BIN"]
    results_transformed = [[] for _ in range(1)]
    for row in results:
        if row[1] != row[1]:
            continue
        for i in range(1, 2):
            results_transformed[i-1].append(row[i])
    print(models[0], str(statistics.median(results_transformed[0])))
    plt.rcParams["figure.figsize"] = (3.5,6)
    ax = seaborn.boxplot(data = results_transformed, palette = [cm.to_hex(plt.cm.Set2(num)) for num in range(5)])
    ax.set_xticklabels(models)
    #plt.ylabel(r"$e$ (average)")
    plt.savefig(os.path.join(path, "box.png"), bbox_inches='tight')
    plt.clf()
    plt.close()


msa_super_dir = "data/lexibench/bin_cross_validation"
raxmlng_super_dir = "data/bin_cross_validation"
plots_super_dir = "data/bin_cross_validation_plots"
all_diff_res = []
diff_headers = ("dataset", "diff_BIN")

# to use same datasets as in study with kappa subsets
studied_datasets = ["lieberherrkhobwa-sinotibetan", "bodtkhobwa-sinotibetan", "chaconcolumbian-chibchan", "felekesemitic-afroasiatic", "kesslersignificance-indoeuropean", "chaconcolumbian-huitotoan", "walworthpolynesian-austronesian", "oskolskayatungusic-tungusic", "robbeetstriangulation-mongolic", "listsamplesize-indoeuropean", "robbeetstriangulation-koreanic", "chaconbaniwa-arawakan", "dhakalsouthwesttibetic-sinotibetan", "savelyevturkic-turkic", "leeainu-ainu", "chaconcolumbian-chocoan", "chaconcolumbian-arawakan", "cals-turkic", "constenlachibchan-misumalpan", "nagarajakhasian-austroasiatic", "mixtecansubgrouping-otomanguean", "robbeetstriangulation-tungusic", "mannburmish-sinotibetan", "liusinitic-sinotibetan", "gerarditupi-tupian", "robbeetstriangulation-japonic"]

for train_ratio in [60]:
    for ds_name in os.listdir(msa_super_dir):
        if ds_name not in studied_datasets:
            continue
        msa_dir = os.path.join(msa_super_dir, ds_name, "cv_" + str(train_ratio))
        target_dir = os.path.join(raxmlng_super_dir, ds_name, "cv_" + str(train_ratio))
        #train_raxml_ng(msa_dir, target_dir)
        #test_raxml_ng(msa_dir, target_dir)
        all_diff_res.append([ds_name] + differences_analysis(msa_dir, target_dir))
    box_plots(all_diff_res, os.path.join(plots_super_dir, "cv_" + str(train_ratio)))
    print(tabulate(all_diff_res, tablefmt="pipe", headers = diff_headers))
