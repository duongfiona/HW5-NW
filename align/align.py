# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing
        n = len(self._seqA)
        m = len(self._seqB)

        # Initialize align matrix with -inf, except for 0,0 = 0 --> to work with later max()
        self._align_matrix = np.full((n+1, m+1), -np.inf, dtype=float)
        self._align_matrix[0,0] = 0

        # Initialize gap matrices, filling in i=0 and j=0 scenarios appropriately
        self._gapA_matrix = np.full((n+1, m+1), -np.inf, dtype=float)
        for j in range(1, m+1): # worst case scenario: we start a gap and just extend it for all seq B (no seq A to align)
            self._gapA_matrix[0, j] = self.gap_open + (j-1)*self.gap_extend

        self._gapB_matrix = np.full((n+1, m+1), -np.inf, dtype=float)
        for i in range(1, n+1): # converse worst case scenario all gaps along seq A
            self._gapB_matrix[i, 0] = self.gap_open + (i-1)*self.gap_extend

        # Initialize backtracing matrices with zeroes 
        self._back = np.zeros((n+1, m+1))
        self._back_A = np.zeros((n+1, m+1))
        self._back_B = np.zeros((n+1, m+1))


        # TODO: Implement global alignment here
        for i in range(1, n+1):
            for j in range(1, m+1):
                ij_pair = (self._seqA[i-1], self._seqB[j-1]) # funky bc of zero indexing
                ij_align_score = self.sub_dict[ij_pair]

                ## (1) Update align matrix
                prev_opts = ( # prior to this ij alignment, could have come from:
                    self._align_matrix[i-1, j-1], # an alignment
                    self._gapA_matrix[i-1, j-1], # a gap in seq A
                    self._gapB_matrix[i-1, j-1] # a gap in seq B
                )

                # new ij alignment score will be best of those previous options, plus the score for ij alignment
                self._align_matrix[i, j] = max(prev_opts) + ij_align_score

                # store backtrace for which previous option new alignment is based on 
                self._back[i, j] = np.argmax(prev_opts) # 0 = prev align, 1 = gap in A, 2 = gap in B


                ## (2) Update gap A matrix
                openA   = self._align_matrix[i, j-1] + self.gap_open + self.gap_extend
                extendA = self._gapA_matrix[i, j-1] + self.gap_extend

                if openA >= extendA: # if opening a new gap after a prev alignment is better than coming from a prev gap 
                    self._gapA_matrix[i, j] = openA
                    self._back_A[i, j] = 0  # 0 = started new gap from prev being an alignment
                else:
                    self._gapA_matrix[i, j] = extendA
                    self._back_A[i, j] = 1  # 1 = extend a prev gap from gapA
                

                ## (3) Update gap B matrix
                openB   = self._align_matrix[i-1, j] + self.gap_open + self.gap_extend
                extendB = self._gapB_matrix[i-1, j] + self.gap_extend

                if openB >= extendB:
                    self._gapB_matrix[i, j] = openB
                    self._back_B[i, j] = 0
                else:
                    self._gapB_matrix[i, j] = extendB
                    self._back_B[i, j] = 1       		
        		    
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # start from the end of alignment
        n = len(self._seqA)
        m = len(self._seqB)

        i = n
        j = m

        final_scores = (
            self._align_matrix[n, m],
            self._gapA_matrix[n, m],
            self._gapB_matrix[n, m]
        )
        self.alignment_score = max(final_scores)
        current_move = np.argmax(final_scores) # 0 = prev align, 1 = gap in A, 2 = gap in B

        # move through backtracing until reaching beginning of each sequence
        while i > 0 or j > 0:
            if current_move == 0:  # came from alignment
                if i > 0 and j > 0:
                    self.seqA_align = self._seqA[i-1] + self.seqA_align # funky seqA/B indexing bc "1st" aa is index 0
                    self.seqB_align = self._seqB[j-1] + self.seqB_align
                    current_move = self._back[i, j]
                elif i > 0:  # only seqA left
                    self.seqA_align = self._seqA[i-1] + self.seqA_align
                    self.seqB_align = "-" + self.seqB_align
                    current_move = 2 # only gaps in seqB remaining
                elif j > 0:  # only seqB left
                    self.seqA_align = "-" + self.seqA_align
                    self.seqB_align = self._seqB[j-1] + self.seqB_align
                    current_move = 1 # only gaps in seqA remaining

                # update position in seq A and B
                i -= 1
                j -= 1

            elif current_move == 1:  # came from gap in A
                if j > 0:
                    self.seqA_align = "-" + self.seqA_align
                    self.seqB_align = self._seqB[j-1] + self.seqB_align

                    # check stored gap behavior (open or extend), update next backtrace move accordingly
                    gap_behavior = self._back_A[i, j] 
                        # 0 = started new gap from prev being an alignment
                        # 1 = extend a prev gap from gapA
                    if gap_behavior == 0: 
                        current_move = 0 # set current_move as from an alignment
                    elif gap_behavior == 1: 
                        current_move = 1 # set current_move as from gap in A
                    j -= 1
    
                else:
                    # at left edge, only seqA left
                    self.seqA_align = self._seqA[i-1] + self.seqA_align
                    self.seqB_align = "-" + self.seqB_align
                    current_move = 2
                    i -= 1

            elif current_move == 2:  # came from gap in B
                if i > 0:
                    self.seqA_align = self._seqA[i-1] + self.seqA_align
                    self.seqB_align = "-" + self.seqB_align

                    # check stored gap behavior, update next backtrace move accordingly
                    gap_behavior = self._back_B[i, j]
                    if gap_behavior == 0: 
                        current_move = 0 # set current_move as from an alignment
                    elif gap_behavior == 1: 
                        current_move = 2 # set current_move as from gap in B
                    i -= 1
                
                else:
                    # at top edge, only seqB left
                    self.seqA_align = "-" + self.seqA_align
                    self.seqB_align = self._seqB[j-1] + self.seqB_align
                    current_move = 1
                    j -= 1

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
