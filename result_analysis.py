import os
import json
from lingdata import database



out_dir = "data/results/"
database.read_config("setup_comparison_lingdata_config.json")
df = database.data()
big_results_dict = {}
for i, row in df.iterrows():
    msa_prefix = "_".join([row["ds_id"], row["source"], row["ling_type"], row["family"]])
    small_results_dict = {}
    ds_path = os.path.join(out_dir, msa_prefix)
    for msa_type in os.listdir(ds_path):
        msa_path = os.path.join(ds_path, msa_type)
        for run_name in os.listdir(msa_path):
            results_path = os.path.join(msa_path, run_name, "results.json")
            if not os.path.isfile(results_path):
                continue
            with open(results_path) as f:
                result = json.load(f)
            small_results_dict[msa_type + "_" + run_name] = result
    big_results_dict[msa_prefix] = small_results_dict

print(big_results_dict)
