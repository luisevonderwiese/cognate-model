import os
import shutil

results_path = "data/results_cognate/"

for ds in os.listdir(results_path):
    shutil.rmtree(os.path.join(results_path, ds, "prototype"))
