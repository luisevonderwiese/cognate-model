import os
import pandas as pd
import matplotlib.pyplot as plt
plots_dir = "data/properties_plots"
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)
df = pd.read_csv("data/lexibench/metadata.tsv", sep = "\t")
plt.hist(df["bin_entropy"], bins = 20)
plt.xlabel("entropy")
plt.ylabel("num datasets")
plt.savefig(os.path.join(plots_dir, "hist_entropies_lexibench.png"))
plt.clf()
plt.close()



