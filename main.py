# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    gap_open = -10.0
    gap_extend = -1.0
    sub_mat_file = "./substitution_matrices/BLOSUM62.mat"
    NW = NeedlemanWunsch(sub_mat_file, gap_open, gap_extend)

    species_seqs = {"Gallus gallus": gg_seq,
                    "Mus musculus": mm_seq,
                    "Balaeniceps rex": br_seq,
                    "Tursiops truncatus": tt_seq}
    
    species_align_scores = {}

    for species, seq in species_seqs.items():
        score, _, _ = NW.align(hs_seq, seq)
        species_align_scores[species] = score

    species_ranked = sorted(species_align_scores, key=lambda k: species_align_scores[k], reverse=True)
    print("Species in order of decreasing similarity of BRD2 to human BRD2")
    print(species_ranked)

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print("\nAlignment scores between human BRD2 and species' BRD2")
    for species in species_ranked:
        print(f"{species} alignment score to human: {species_align_scores[species]}")
    

if __name__ == "__main__":
    main()