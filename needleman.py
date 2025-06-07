import numpy as np

def generate_array(strand1, strand2, match, mismatch, gap):
    """
    Generates the Needleman-Wunsch scoring matrix for global sequence alignment.

    Parameters:
        strand1 (str): First nucleotide sequence.
        strand2 (str): Second nucleotide sequence.
        match (int): Score for a match.
        mismatch (int): Penalty for a mismatch.
        gap (int): Penalty for a gap.

    Returns:
        np.ndarray: Filled scoring matrix.
    """
    strand1 = strand1.upper()
    strand2 = strand2.upper()
    rows = len(strand1) + 1
    cols = len(strand2) + 1
    matrix = np.zeros((rows, cols), dtype=int)

    #initialize first row and column with cumulative gap penalties
    for i in range(1, rows):
        matrix[i][0] = matrix[i - 1][0] + gap
    for j in range(1, cols):
        matrix[0][j] = matrix[0][j - 1] + gap

    #fill in the scoring matrix based on the best path to each cell
    for i in range(1, rows):
        for j in range(1, cols):
            if strand1[i - 1] == strand2[j - 1]:
                diag = matrix[i - 1][j - 1] + match
            else:
                diag = matrix[i - 1][j - 1] + mismatch
            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap
            matrix[i][j] = max(diag, up, left)

    return matrix

def analyze_alignment(matrix, strand1, strand2, match, mismatch, gap):
    """
    Performs traceback to retrieve aligned sequences from a Needleman-Wunsch matrix.

    Parameters:
        matrix (np.ndarray): Scoring matrix.
        strand1 (str): First nucleotide sequence.
        strand2 (str): Second nucleotide sequence.
        match (int): Score for a match.
        mismatch (int): Penalty for a mismatch.
        gap (int): Penalty for a gap.

    Returns:
        tuple[list[str], list[str]]: The aligned sequences with gaps.
    """
    i, j = len(strand1), len(strand2)
    aligned1, aligned2 = [], []

    #trace back from bottom-right to top-left of matrix
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            #determine if current cell is from a diagonal match/mismatch
            match_score = match if strand1[i - 1] == strand2[j - 1] else mismatch
            if matrix[i][j] == matrix[i - 1][j - 1] + match_score:
                aligned1.append(strand1[i - 1])
                aligned2.append(strand2[j - 1])
                i -= 1
                j -= 1
                continue

        #checck if movement was from above (gap in strand2)
        if i > 0 and matrix[i][j] == matrix[i - 1][j] + gap:
            aligned1.append(strand1[i - 1])
            aligned2.append('-')
            i -= 1
        else:
            #movement must be from left (gap in strand1)
            aligned1.append('-')
            aligned2.append(strand2[j - 1])
            j -= 1

    return aligned1[::-1], aligned2[::-1]

def alignment_score(aligned1, aligned2, match, mismatch, gap):
    """
    Calculates the total alignment score for two aligned sequences.

    Parameters:
        aligned1 (list[str]): First aligned sequence.
        aligned2 (list[str]): Second aligned sequence.
        match (int): Score for a match.
        mismatch (int): Penalty for a mismatch.
        gap (int): Penalty for a gap.

    Returns:
        int: Total alignment score.
    """
    score = 0
    for a, b in zip(aligned1, aligned2):
        if a == '-' or b == '-':
            score += gap
        elif a == b:
            score += match
        else:
            score += mismatch
    return score

def pairwise_distance_matrix(sequences, match=1, mismatch=-1, gap=-2):
    """
    Builds a pairwise distance matrix using Needleman-Wunsch alignment scores.

    Parameters:
        sequences (list[str]): List of nucleotide sequences.
        match (int): Score for a match.
        mismatch (int): Penalty for a mismatch.
        gap (int): Penalty for a gap.

    Returns:
        list[list[int]]: Symmetric matrix of alignment scores.
    """
    n = len(sequences)
    matrix = [[0] * n for _ in range(n)]

    #compute upper triangle (and mirror to lower triangle)
    for i in range(n):
        for j in range(i + 1, n):
            m = generate_array(sequences[i], sequences[j], match, mismatch, gap)
            a1, a2 = analyze_alignment(m, sequences[i], sequences[j], match, mismatch, gap)
            score = alignment_score(a1, a2, match, mismatch, gap)
            matrix[i][j] = matrix[j][i] = score

    return matrix
