import os
from tabulate import tabulate

import util

msa_super_dir = "data/lexibench/character_matrices/"

result1 = []
result2 = []
for dataset in os.listdir(msa_super_dir):
    row1 = [dataset]
    row2 = [dataset]
    msa_path = os.path.join(msa_super_dir, dataset, "bin.phy")
    if not os.path.isfile(msa_path):
        continue
    for kappa in range(2, 7):
        msa_path = os.path.join(msa_super_dir, dataset, "bv_part_" + str(kappa) + ".phy")
        if os.path.isfile(msa_path):
            alignment = util.safe_msa_read(msa_path)
            num_concepts = alignment.get_alignment_length()
            num_languages = len(alignment)
            row1.append(num_concepts)
            row2.append(num_concepts / num_languages)
        else:
            row1.append(0)
            row2.append(0)

    result1.append(row1)
    result2.append(row2)


headers = ["dataset"] + ["kappa=" + str(s) for s in range(2, 7)]

table = tabulate(result1, headers, tablefmt = "latex")
print(table)

table = tabulate(result2, headers, tablefmt = "latex")
print(table)

