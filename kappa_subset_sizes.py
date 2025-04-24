import matplotlib.pyplot as plt
import matplotlib.colors as cm
import os
import seaborn
import util


plots_super_dir = "data/properties_plots/"
msa_super_dir = "data/lexibench/character_matrices/"


if not os.path.isdir(plots_super_dir):
    os.makedirs(plots_super_dir)
subset_sizes = [[] for kappa in range(7)]
for dataset in os.listdir(msa_super_dir):
    for kappa in range(2, 7):
        msa_path = os.path.join(msa_super_dir, dataset, "bv_part_" + str(kappa) + ".phy")
        if os.path.isfile(msa_path):
            alignment = util.safe_msa_read(msa_path)
            subset_sizes[kappa].append(alignment.get_alignment_length())
        else:
            subset_sizes[kappa].append(0)



ax = seaborn.boxplot(data = subset_sizes[2:], palette = [cm.to_hex(plt.cm.Set1(num)) for num in range(5)])
ax.set_xticklabels(range(2, 7))
plt.ylabel(r"sizes of the $\kappa$-subsets")
plt.xlabel(r"$\kappa$")
plt.savefig(os.path.join(plots_super_dir, "kappa_subset_sizes.png"))
plt.clf()
plt.close()

