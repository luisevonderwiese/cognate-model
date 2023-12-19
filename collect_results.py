
from ete3 import Tree
import os
import json
#from database import *

def get_alpha(inference_log_path):
    with open(inference_log_path, "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("   Rate heterogeneity:"):
            parts = line.split(",  ")
            if len(parts) < 2:
                return float("nan")
            return float(parts[1].split(" ")[1])
    return float('nan')

def get_num_topos(rf_distances_log_path):
    with open(rf_distances_log_path, "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("Number of unique topologies in this tree set:"):
            return int(line.split(": ")[1])
    return float('nan')

def get_average_rf_distance(rf_distances_log_path):
    with open(rf_distances_log_path, "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("Average relative RF distance in this tree set:"):
            return float(line.split(": ")[1])
    return float('nan')

def get_scores(inference_log_path):
    with open(inference_log_path, "r") as logfile:
        lines = logfile.readlines()
    for line in lines:
        if line.startswith("AIC"):
            parts = line.split(" / ")
            scores = []
            for part in parts:
                scores.append(float(part.split(" ")[2]))
            return scores
    return [float('nan'), float('nan'), float('nan')]


#db.init(snakemake.output.database)
#db.connect()
inference_log_path = snakemake.input.inference_log
rf_distances_log_path = snakemake.input.rf_distances_log
best_tree_path = snakemake.input.best_tree

outpath = snakemake.output.results

results = {}
with open(best_tree_path, "r") as best_tree_file:
    results["best_tree_newick"] = best_tree_file.read()
results["alpha"] = get_alpha(inference_log_path)
scores = get_scores(inference_log_path)
results["AIC"] = scores[0]
results["AICc"] = scores[1]
results["BIC"] = scores[2]
results["num_topos"] = get_num_topos(rf_distances_log_path)
results["average_rf_distance"] = get_average_rf_distance(rf_distances_log_path)
json.dump(results, open(outpath,'w+'))


