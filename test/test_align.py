# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    sub_mat_file = "./substitution_matrices/BLOSUM62.mat"
    gap_open = -10.0
    gap_extend = -1.0

    NW = NeedlemanWunsch(sub_mat_file, gap_open, gap_extend)
    _, _, _ = NW.align(seq1, seq2)

    # check that all -inf in the main body of matrices were replaced
    assert np.isfinite(NW._align_matrix[0, 0]), "Initial align_matrix value still -inf"
    assert np.all(np.isfinite(NW._align_matrix[1:, 1:])), "Some interior align_matrix values are still -inf!"

    assert np.all(np.isfinite(NW._gapA_matrix[0, 1:])), "Initial gapA_matrix values are still -inf!"
    assert np.all(np.isfinite(NW._gapA_matrix[2:, 2:])), "Some interior gapA_matrix values are still -inf!"

    assert np.all(np.isfinite(NW._gapB_matrix[1:, 0])), "Initial gapB_matrix values are still -inf!"
    assert np.all(np.isfinite(NW._gapB_matrix[2:, 2:])), "Some interior gapB_matrix values are still -inf!"

    # assert that initial gap boundary conditions were handled correctly
    assert NW._gapA_matrix[0, 1] == gap_open + gap_extend
    assert NW._gapA_matrix[0, 2] == gap_open + 2*gap_extend

    assert NW._gapB_matrix[1, 0] == gap_open + gap_extend
    assert NW._gapB_matrix[2, 0] == gap_open + 2*gap_extend
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    sub_mat_file = "./substitution_matrices/BLOSUM62.mat"
    gap_open = -10.0
    gap_extend = -1.0

    NW = NeedlemanWunsch(sub_mat_file, gap_open, gap_extend)
    align_score, seq3_align, seq4_align = NW.align(seq3, seq4)

    assert align_score == 17
    assert seq3_align == "MAVHQLIRRP"
    assert seq4_align == "M---QLIRHP"



