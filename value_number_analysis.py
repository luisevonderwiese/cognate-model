import pandas as pd
import os
import matplotlib.pyplot as plt
from tabulate import tabulate
import math
from scipy import stats
import numpy as np

from lingdata import database





plots_dir = "data/properties_plots"
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)
msa_super_dir = "data/lexibench/msa"
vn_path = "data/lexibench/value_number_counts.csv"
with open(vn_path, "r") as vn_file:
    lines = vn_file.readlines()



all_counts = {}
for line in lines:
    parts = line.split(";")
    dataset = parts[0]
    counts_string = parts[1][:-1].strip("[]")
    if counts_string == "[]":
        all_counts[dataset] = []
    else:
        all_counts[dataset] = [int(el) for el in counts_string.split(", ")]


max_num = max([len(counts) for dataset, counts in all_counts.items()])
fig,ax = plt.subplots(figsize=(15, 10))
x = range(len(all_counts))
y_old = [0 for el in x]
for num in range(max_num):
    y_new = []
    for dataset, counts in all_counts.items():
        if len(counts) > num:
            y_new.append(counts[num])
        else:
            y_new.append(0)
    ax.bar(x, y_new, bottom=y_old, label = str(num))
    for i in x:
        y_old[i] = y_old[i] + y_new[i]
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
      box.width, box.height * 0.9])

ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=11)

plt.savefig(os.path.join(plots_dir, "value_number_analysis.png"))
plt.clf()
