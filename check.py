import os
import numpy as np
import matplotlib.pyplot as plt

plots_dir = "data/plots"
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)

def to_matrix(rates):
    if len(rates) == 1:
        x = 2
    elif len(rates) == 6:
        x = 4
    elif len(rates) == 28:
        x = 8
    elif len(rates) == 120:
        x = 16
    elif len(rates) == 496:
        x = 32
    elif len(rates) == 2016:
        x = 64
    else:
        print("Illegal num rates " + str(len(rates)))
    matrix = [[0 for i in range(x)] for j in range(x)]
    for col_idx in range(x):
        for row_idx in range(col_idx):
            idx = sum(range(x-row_idx, x)) + col_idx - row_idx - 1
            rate = rates[idx]
            if rate > 30:
                print(rate)
                print(row_idx)
                print(col_idx)
                print(" ")
            matrix[col_idx][row_idx] = rate
            matrix[row_idx][col_idx] = rate
    return matrix

for ds in os.listdir("data/results"):
    if os.path.isfile(os.path.join("data/results", ds, "prototype", "GTR", "inference.raxml.bestTree")):
            print(ds)
            with open(os.path.join("data/results", ds, "prototype", "GTR", "inference.raxml.log"), "r") as logfile:
                lines = logfile.readlines()
            for line in lines:
                if line.startswith("   Substitution rates"):
                    parts = line.split(": ")[1].split(" ")[:-1]
                    rates = [float(part) for part in parts]
                    break
            rate_matrix = to_matrix(rates)
            rate_matrix = np.array(rate_matrix)
            plt.imshow(rate_matrix)
            plt.colorbar()
            plt.savefig(os.path.join(plots_dir, ds))
            plt.clf()
print("COG")
for ds in os.listdir("data/results"):
    if os.path.isfile(os.path.join("data/results", ds, "cognate", "COG", "inference.raxml.bestTree")):
        print(ds)

