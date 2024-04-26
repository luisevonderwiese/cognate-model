from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter
import random
import os
import math

def split_indices(num_sites, num_samples = 10, ratio = 0.6):
    num_sites_train = math.ceil(num_sites * ratio)
    res = []
    for i in range(num_samples):
        l = [_ for _ in range(num_sites)]
        random.shuffle(l)
        train_indices = l[:num_sites_train]
        res.append(train_indices)
    return res

def empty_align(ref_align):
    new_records = [SeqRecord([]), id=ref_align[i].id) for i in range(len(ref_align))]
    return MultipleSeqAlignment(new_records, annotations={}, column_annotations={})

def concat_align(a1, a2):
    new_sequences = []
    assert(len(a1) == len(a2))
    for i in range(len(a1)):
        assert(a1[i].id == a2[i].id)
        seq1 = a1[i].seq
        seq2 = a2[i].seq
        new_sequences.append(seq1 + seq2)
    new_records = [SeqRecord([new_sequences[i]]), id=a1[i].id) for i in range(len(a1))]
    return MultipleSeqAlignment(new_records, annotations={}, column_annotations={})




def create_samples(kappa, msa_dir):
    bin_msa_type = "bin_part_" + str(kappa)
    prototype_msa_type = "prototype_part_" + str(kappa)
    bin_msa_path = os.path.join(msa_dir, bin_msa_type + ".phy")
    prototype_msa_path = os.path.join(msa_dir, prototype_msa_type + ".phy")
    if not os.path.isfile(bin_msa_path) or not os.path.isfile(prototype_msa_path):
        return
    with open(prototype_msa_path, "r") as msa_file:
        num_sites = int(msa_file.readlines()[0].split(" ")[2])
    with open(bin_msa_path, "r") as msa_file:
        num_sites_bin = int(msa_file.readlines()[0].split(" ")[2])
    assert(num_sites_bin == kappa * num_sites)
    try:
        bin_align = AlignIO.read(bin_msa_path, "phylip-relaxed")
    except:
        print(msa_dir)
        return
    try:
        prototype_align = AlignIO.read(bin_msa_path, "phylip-relaxed")
    except:
        print(msa_dir)
        return
    indices_list = split_indices(num_sites)
    for (t, train_indices) in enumerate(indices_list):
        bin_train_align = empty_align(bin_align)
        bin_test_align = empty_align(bin_align)
        prototype_train_align = empty_align(prototype_align)
        prototype_test_align = empty_align(prototype_align)
        for s in range(num_sites):
            if s in train_indices :
                prototype_train_align = concat_align(prototype_train_align, prototype_align[:, s])
                bin_train_align += concat_align(bin_train_align, bin_align[:, s*kappa : (s+1) * kappa])
            else:
                prototype_test_align = bin_train_align(prototype_test_align, prototype_align[:, s])
                bin_test_align += bin_train_align(bin_test_align, bin_align[:, s*kappa : (s+1) * kappa])
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
    create_samples(kappa, msa_dir)
