import os
import shutil

results_path = "data/results_cognate_single_noforce/"
msa_dir = "data/lingdata_cognate/msa/"
for ds in os.listdir(results_path):
    print(ds)
    for x in range(2, 7):
        msa_path = os.path.join(msa_dir, ds , "prototype_part_" + str(x) + ".phy")
        if not os.path.isfile(msa_path):
            cur_dir = os.path.join(results_path, ds, "prototype_part" + str(x))
            print(cur_dir)
            if os.path.isdir(cur_dir):
                shutil.rmtree(cur_dir)
