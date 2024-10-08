from Bio import AlignIO
import matplotlib.pyplot as plt
import os
import math

symbols = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "G",
         "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X",
         "Y", "Z", "!", "\"", "#", "$", "%", "&", "\'", "(", ")", "*", "+", ",", "/", ":",
         ";", "<", "=", ">", "@", "[", "\\", "]", "^", "_", "{", "|", "}", "~"]

symbols_dict  = dict((symbol, (i+1).bit_count()) for (i, symbol) in enumerate(symbols))

plots_super_dir = "data/properties_plots/"
msa_super_dir = "data/lexibench/msa/"


all_counts = []
plots_dir = os.path.join(plots_super_dir, "symbol_counts")
if not os.path.isdir(plots_dir):
    os.makedirs(plots_dir)
for dataset in os.listdir(msa_super_dir):
    for kappa in range(2, 7):
        counts_plots_dir = os.path.join(plots_super_dir, str(kappa), "counts")
        symbol_counts_plots_dir = os.path.join(plots_super_dir, str(kappa), "symbol_counts")
        for d in [counts_plots_dir, symbol_counts_plots_dir]:
            if not os.path.isdir(d):
                os.makedirs(d)
        msa_path = os.path.join(msa_super_dir, dataset, "bv_part_" + str(kappa) + ".phy")
        if os.path.isfile(msa_path):
            counts = [0 for _ in range(kappa + 1)]
            limit = int(math.pow(2, kappa) - 1)
            symbol_counts = dict((symbol, 0) for symbol in symbols[:limit])
            alignment = AlignIO.read(msa_path, "phylip-relaxed")
            for record in alignment:
                for el in record.seq:
                    if el in ['-', '?']:
                        continue
                    counts[symbols_dict[el]] += 1
                    symbol_counts[el] += 1
            fig, ax = plt.subplots(figsize = (15, 6))
            #ax1.bar(range(kappa + 1), counts)
            #ax1.set_xlabel("#cognate classes")
            #ax1.set_ylabel("#language-concept-pairs")
            #ax1.set_title("(a)")
            labels = []
            for i, symbol in enumerate(symbol_counts.keys()):
                bitvector = bin(i+1)[2:].zfill(kappa)
                labels.append(symbol + " [" + str(bitvector) + "]")
            ax.bar(labels, symbol_counts.values())
            ax.set_xlabel("symbol")
            ax.set_ylabel("#occurrences")
            plt.xticks(rotation=30, ha='right')
            plt.subplots_adjust(bottom=0.19)
            plt.savefig(os.path.join(plots_dir, str(kappa) + "_" + dataset + ".png"))
            plt.clf()
            plt.close()
