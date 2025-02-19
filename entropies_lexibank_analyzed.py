import os
import traceback

from cognate import CognateData
from lingpy import *

import matplotlib.pyplot as plt


def detect_cognates(wordlist_path, wordlist_cognate_path):
    assert(os.path.isfile(wordlist_path))
    print("Detecting cognates...")
    lex = LexStat(wordlist_path, check = False)
    lex.get_scorer()
    lex.cluster(method="lexstat", threshold=0.6, ref="cognates")
    lex.output('tsv', filename=wordlist_cognate_path.split(".")[0], ignore = "all", prettify = False)
    

base_dir = os.path.join("results", "lexibank-analyzed_families")
wordlist_dir = "data/lexibank-analyzed/wordlists/families/"
wordlist_cognate_dir = "data/lexibank-analyzed/wordlists_cognate/families/"
plots_dir = "data/properties_plots"
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)
families = set()
for file_name in os.listdir(wordlist_dir):
    families.add(file_name.split("_")[0])

all_entropies = []
for file_name in os.listdir(wordlist_dir):
    full_name = file_name.split("_")[0]
    wordlist_path = os.path.join(wordlist_dir, file_name)
    if not os.path.isfile(wordlist_path):
        continue
    print(full_name)
    wordlist_cognate_path = os.path.join(wordlist_cognate_dir, full_name + "_wordlist_cognate.tsv")
    if not os.path.isfile(wordlist_cognate_path):
        try:
            detect_cognates(wordlist_path, wordlist_cognate_path)
        except Exception as e:
            traceback.print_exc()
            print(e)
            continue
    cd = CognateData.from_edictor_tsv(wordlist_cognate_path)
    if cd.num_languages() < 4:
        print(full_name, "too small")
        continue
    all_entropies.append(cd.bin_entropy())

plt.hist(all_entropies, bins = 20)
plt.xlabel("entropy")
plt.ylabel("num datasets")
plt.savefig(os.path.join(plots_dir, "hist_entropies_lexibank-analyzed.png"))
plt.clf()
plt.close()

