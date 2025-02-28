import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.colors as cm


plots_dir = "data/properties_plots"
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)
metadata_path = "data/lexibench/metadata.tsv"
metadata_df = pd.read_csv(metadata_path, sep = "\t")

all_counts = {}
for i, row in metadata_df.iterrows():
    dataset = row["dataset"]
    counts_string = row["polymorphism_counts"][:-1].strip("[]")
    if counts_string == "[]":
        all_counts[dataset] = []
    else:
        all_counts[dataset] = [int(el) for el in counts_string.split(", ")]


all_counts = dict(sorted(all_counts.items()))
max_num = max([len(counts) for dataset, counts in all_counts.items()])
_, ax = plt.subplots(figsize=(15, 10))
x = [dataset for dataset,counts in all_counts.items()]
y_old = [0 for el in x]
for num in range(max_num):
    y_new = []
    for dataset, counts in all_counts.items():
        if len(counts) > num:
            y_new.append(counts[num])
        else:
            y_new.append(0)
    if num < 10:
        color = cm.to_hex(plt.cm.tab10(num))
    else:
        color = "black"
    ax.bar(x, y_new, bottom=y_old, label=r'$\nu=' + str(num) + '$', color = color)
    for i in range(len(x)):
        y_old[i] = y_old[i] + y_new[i]
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
      box.width, box.height * 0.9])

ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
          fancybox=True, shadow=True, ncol=11)
plt.xticks(rotation=30, ha='right')
plt.xlabel("datasets")
plt.ylabel("#language-concept-pairs")

plt.savefig(os.path.join(plots_dir, "polymorphism_counts.png"))
plt.clf()
