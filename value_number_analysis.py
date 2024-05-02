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

config_path = "cognate_lingdata_config.json"



database.read_config(config_path)
df = database.data()


big_matrix = []
for i,row in df.iterrows():
    matrix_string = row["value_number_matrix"]
    matrix = []
    for el in matrix_string[2:-2].split("], ["):
        if el == "[]" or el == "":
            matrix.append([])
        else:
            matrix.append([int(inner_el) for inner_el in el.split(", ")])

    #matrix = [[] if el == "[]" else [int(inner_el) for inner_el in el.split(", ")] for el in matrix_string.strip("[]").split("], [")]
    while len(big_matrix) < len(matrix):
        big_matrix.append([])
    for i, counts in enumerate(matrix):
        while len(counts) >= len(big_matrix[i]):
            big_matrix[i].append(0)
        for j, count in enumerate(counts):
            big_matrix[i][j] += count


all_counts = [row["value_number_counts"] for i,row in df.iterrows()]
max_num = max([len(counts) for counts in all_counts])
fig,ax = plt.subplots(figsize=(15, 10))
x = range(len(all_counts))
y_old = [0 for el in x]
for num in range(max_num):
    y_new = []
    for counts in all_counts:
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
