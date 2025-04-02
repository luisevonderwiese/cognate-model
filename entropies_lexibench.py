import os
import pandas as pd
import matplotlib.pyplot as plt

import util

plots_dir = "data/properties_plots"
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)

entropies = []
metadata_df = pd.read_csv("data/lexibench/character_matrices/stats.tsv", sep = "\t")
for i, row in metadata_df.iterrows():
    try:
        align = util.safe_msa_read(row["bin.phy"])
    except:
        continue
    entropies.append(util.bin_entropy(align))


plt.hist(entropies, bins = 20)
plt.xlabel("entropy")
plt.ylabel("num datasets")
plt.savefig(os.path.join(plots_dir, "hist_entropies_lexibench.png"))
plt.clf()
plt.close()
