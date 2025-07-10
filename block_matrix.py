import os

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.AlignIO.PhylipIO import RelaxedPhylipWriter

kappa = 3

# determine total_num_sites

alignments = {}
total_num_sites = 0
base_dir = "data/lexibench/character_matrices"
for ds in os.listdir(base_dir):
    msa_path = os.path.join(base_dir, ds, "bv_part_" + str(kappa) + ".phy")
    if not os.path.isfile(msa_path):
        print("Skipping", ds)
        continue
    try:
        alignment = AlignIO.read(msa_path, "phylip-relaxed")
    except Exception as e:
        print(e)
        print("Error in", ds)
        continue
    total_num_sites += alignment.get_alignment_length()
    alignments[ds] = alignment

new_records = []
passed_sites = 0
for ds, alignment in alignments.items():
    num_sites = alignment.get_alignment_length()
    for record in alignment:
        seq = ("-" * passed_sites) + record.seq + ("-" * (total_num_sites - passed_sites - num_sites))
        new_records.append(SeqRecord(seq, id=ds + "_" + record.id))
    passed_sites += num_sites

block_alignment = MultipleSeqAlignment(new_records, annotations={}, column_annotations={})
with open("data/block_matrix_" + str(kappa) + ".phy", "w+", encoding = "utf-8") as f:
    writer = RelaxedPhylipWriter(f)
    writer.write_alignment(block_alignment)
