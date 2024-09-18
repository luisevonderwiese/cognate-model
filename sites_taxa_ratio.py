import lingdata.database as database
from Bio import AlignIO
import matplotlib.pyplot as plt
import os
plots_super_dir = "data/properties_plots/"
msa_super_dir = "data/lexibench/msa/"
all_counts = []
datasets = [dataset for dataset in os.listdir(msa_super_dir)]
for dataset in os.listdir(msa_super_dir):
    counts = [0, 0]
    for kappa in range(2, 7):
        msa_path = os.path.join(msa_super_dir, dataset, "prototype_part_" + str(kappa) + ".phy")
        if os.path.isfile(msa_path):
            alignment = AlignIO.read(msa_path, "phylip-relaxed")
            counts.append(alignment.get_alignment_length())
        else:
            counts.append(0)
    all_counts.append(counts)





fig,ax = plt.subplots(figsize=(15, 10))
y_old = [0 for el in datasets]
for num in range(2, 7):
    y_new = []
    for counts in all_counts:
        y_new.append(counts[num])
    ax.bar(datasets, y_new, bottom=y_old, label = str(num))
    for i in range(len(datasets)):
        y_old[i] = y_old[i] + y_new[i]
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
      box.width, box.height * 0.9])
plt.xticks(rotation=30, ha='right')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
          fancybox=True, shadow=True, ncol=11)
plt.xlabel("datasets")
plt.ylabel("#concepts")

plt.savefig(os.path.join(plots_super_dir, "partition_sizes.png"))
plt.clf()

plots_dir = os.path.join(plots_super_dir, "sites_taxa_ratio")
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)
fig = plt.figure(figsize=(13.5, 9))
axes = []
axes.append(plt.subplot2grid(shape=(2,6), loc=(0,1), colspan=2, fig=fig))
axes.append(plt.subplot2grid((2,6), (0,3), colspan=2, fig=fig))
axes.append(plt.subplot2grid((2,6), (1,0), colspan=2, fig=fig))
axes.append(plt.subplot2grid((2,6), (1,2), colspan=2, fig=fig))
axes.append(plt.subplot2grid((2,6), (1,4), colspan=2, fig=fig))
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
for ax in axes:
    ax.set_xlim([-1, 80])
    ax.set_ylim([-1, 80])
    ax.axline([0, 0], slope=1, color = 'lightgray', linewidth = 1, linestyle = "--")
axes[3].set_xlabel("#languages")
axes[0].set_ylabel("#colums")
axes[2].set_ylabel("#colums")
all_handles = []
all_labels = []
for kappa in range(2, 7):
    x_values = []
    y_values = []
    for dataset in os.listdir(msa_super_dir):
        msa_path = os.path.join(msa_super_dir, dataset, "prototype_part_" + str(kappa) + ".phy")
        if os.path.isfile(msa_path):
            with open(msa_path, "r") as msa_file:
                parts = msa_file.readlines()[0].split(" ")
            x_values.append(int(parts[1]))
            y_values.append(int(parts[2]))
    axes[kappa-2].scatter(x_values, y_values, s = 8, label=r'$\kappa=' + str(kappa) + '$', color = cycle[kappa-2])
    handles, labels = axes[kappa-2].get_legend_handles_labels()
    all_handles += handles
    all_labels += labels
fig.legend(all_handles, all_labels, loc='lower center', ncol = 5)
#box = ax.get_position()
#ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), fancybox=True, shadow=True, ncol=11)
plt.savefig(os.path.join(plots_dir, "scatter_all.png"))
plt.clf()
