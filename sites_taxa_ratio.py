import lingdata.database as database
from Bio import AlignIO
import matplotlib.pyplot as plt
import os
plots_super_dir = "data/properties_plots/"
config_path = "cognate_lingdata_config.json"
database.read_config(config_path)
#database.download()
#database.compile()
df = database.data()
all_counts = []
for i, row in df.iterrows():
    counts = [0, 0]
    for x in range(2, 7):
        msa_path = row["msa_paths"]["prototype_part_" + str(x)]
        if msa_path == msa_path and os.path.isfile(msa_path):
            alignment = AlignIO.read(msa_path, "phylip-relaxed")
            counts.append(alignment.get_alignment_length())
        else:
            counts.append(0)
    all_counts.append(counts)





fig,ax = plt.subplots(figsize=(40, 30))
x = range(len(all_counts))
y_old = [0 for el in x]
for num in range(2, 7):
    y_new = []
    for counts in all_counts:
        y_new.append(counts[num])
    ax.bar(x, y_new, bottom=y_old, label = str(num))
    for i in x:
        y_old[i] = y_old[i] + y_new[i]
ax.legend()
plt.savefig(os.path.join(plots_super_dir, "partition_sizes.png"))
plt.clf()

plots_dir = os.path.join(plots_super_dir, "sites_taxa_ratio")
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)
fig,ax = plt.subplots(figsize=(10, 7.5))
for x in range(2, 7):
    x_values = []
    y_values = []
    for i, row in df.iterrows():
        x_values.append(row["num_taxa"])
        msa_path = row["msa_paths"]["prototype_part_" + str(x)]
        if msa_path == msa_path and os.path.isfile(msa_path):
            alignment = AlignIO.read(msa_path, "phylip-relaxed")
            y_values.append(alignment.get_alignment_length())
        else:
            y_values.append(0)
    plt.scatter(x_values, y_values)
    plt.axline([0, 0], slope=1, color = 'lightgray', linewidth = 1, linestyle = "--")
    plt.xlabel("num taxa")
    plt.ylabel("num sites")
    plt.savefig(os.path.join(plots_dir, "scatter_" + str(x) + ".png"))
    plt.clf()
