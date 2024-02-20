
rule raxml:
    output:
        best_tree     = f"{results_dir}inference.raxml.bestTree",
        ml_trees     = f"{results_dir}inference.raxml.mlTrees",
        inference_log           = f"{results_dir}inference.raxml.log",
    params:
        prefix  = results_dir,
        msa     = lambda wildcards: msa_path_dict[wildcards.run_prefix],
        model   = lambda wildcards: model_dict[wildcards.run_prefix],
        prob_msa = lambda wildcards: prob_msa_dict[wildcards.run_prefix],
        threads = config["software"]["raxml-ng"]["threads"],
        seed    = seed,
    log:
        f"{results_dir}inference.snakelog",
    shell:
        "{raxmlng_command} "
        "--msa {params.msa} "
        "--model {params.model} "
        "--prefix {params.prefix}inference "
        "--seed {params.seed} "
        "--threads {params.threads} "
        "--prob-msa {params.prob_msa} --redo "
        "> {output.inference_log} "


rule pythia:
    output:
        pythia_difficulty     = f"{results_dir}pythia.difficulty",
    params:
        msa     = lambda wildcards: pythia_msa_path_dict[wildcards.run_prefix],
        raxmlng_path = raxmlng_command,
        predictor_path = pythia_predictor,

    log:
        f"{results_dir}pythia.snakelog",
    script:
        "run_pythia.py"

rule rf_distances:
    input:
        ml_trees = rules.raxml.output.ml_trees
    output:
        rf_distances      = f"{results_dir}rf_distances.raxml.rfDistances",
        rf_distances_log  = f"{results_dir}rf_distances.raxml.log",
    params:
        prefix = results_dir
    log:
        f"{results_dir}rf_distances.snakelog",
    shell:
        "{raxmlng_command} "
        "--rfdist "
        "--tree {input.ml_trees} "
        "--prefix {params.prefix}rf_distances "
        "> {output.rf_distances_log} "


rule result_collection:
    input:
        inference_log =  rules.raxml.output.inference_log,
        best_tree = rules.raxml.output.best_tree,
        rf_distances_log = rules.rf_distances.output.rf_distances_log,
    output:
        results      = f"{results_dir}results.json",
    log:
        f"{results_dir}result_collection.snakelog",
    script:
        "collect_results.py"




