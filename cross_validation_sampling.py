from Bio import AlignIO
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter
import random
import os
import math

def split_indices(num_sites, num_samples = 10, ratio = 0.6):
    num_sites_train = math.ceil(num_sites * ratio)
    res = []
    for i in range(num_samples):
        shuffle = random.shuffle([_ for _ in range(num_sites)])
        train_indices = shuffle[:num_chars_train]
        res.append(train_indices)
    return res


def create_samples(kappa, msa_dir):
    bin_msa_type = "bin_part_" + str(kappa)
    prototype_msa_type = "prototype_part_" + str(kappa)
    bin_msa_path = os.path.join(msa_dir, bin_msa_type + ".phy")
    prototype_msa_type = os.path.join(msa_dir, prototype_msa_type + ".phy")
    with open(prototype_msa_type, "r") as msa_file:
        num_sites = int(msa_file.readlines()[0].split(" ")[2])
    with open(bin_msa_path, "r") as msa_file:
        num_sites_bin = int(msa_file.readlines()[0].split(" ")[2])
    assert(num_sites_bin == kappa * num_sites)
    try:
        bin_align = AlignIO.read(bin_msa_path, "phylip-relaxed")
    except:
        print(msa_dir)
        continue
    try:
        prototype_align = AlignIO.read(bin_msa_path, "phylip-relaxed")
    except:
        print(msa_dir)
        continue

    indices_list = split_indices(num_sites)
    for (t, train_indices) in enumerate(indices_list):
        bin_train_align = Bio.Align.MultipleSeqAlignment([])
        bin_test_align = Bio.Align.MultipleSeqAlignment([])
        prototype_train_align = Bio.Align.MultipleSeqAlignment([])
        prototype_test_align = Bio.Align.MultipleSeqAlignment([])
        for s in range(num_sites):
            if s in train_indices :
                prototype_train_align += prototype_align[:, s]
                bin_train_align += bin_align[:, s*kappa : (s+1) * kappa]
            else:
                prototype_train_align += prototype_align[:, s]
                bin_train_align += bin_align[:, s*kappa : (s+1) * kappa]
        with open(os.path.join(msa_dir, bin_msa_type + "_cv_train_" + str(t) + ".phy"),"w+") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(bin_train_align)
        with open(os.path.join(msa_dir, bin_msa_type + "_cv_test_" + str(t) + ".phy"),"w+") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(bin_test_align)
        with open(os.path.join(msa_dir, prototype_msa_type + "_cv_train_" + str(t) + ".phy"),"w+") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(prototype_train_align)
        with open(os.path.join(msa_dir, prototype_msa_type + "_cv_test_" + str(t) + ".phy"),"w+") as f:
            writer = RelaxedPhylipWriter(f)
            writer.write_alignment(prototype_test_align)



msa_super_dir = "data/lingdata_cognate/msa"
kappa = 3
for ds_name in os.listdir(msa_super_dir):
    msa_dir = os.path.join(msa_super_dir, ds_name)
    create_samples(msa_dir, kappa)
