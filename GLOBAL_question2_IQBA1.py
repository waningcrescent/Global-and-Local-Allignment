import numpy as np

# YES , the result will change when the scoring scheme is changed, as the value selected for each cell is heavily
# dependent on the scoring scheme(s) and we can see the result in the output of this code

# Initializing the two sequences as arrays:
s1 = "GATGCGCAG"
s2 = "GGCAGTA"


# the scoring scheme is as follows:
Match = +2
Mismatch = -1
Gap = -3

# initializing the matrix :
n = len(s1)
m = len(s2)

# As per the Needleman-Wunsch algorithm, the initialization for the matrix is to assign gap penalty to first row & col
# matrix[i][0] = gap_penalty*i & matrix[0][j] = gap_penalty*j
matrix = np.zeros((len(s1) + 1, len(s2) + 1))
# initialize first row and column with gap penalties:
for i in range(1, n + 1):
    matrix[i][0] = Gap * i
for j in range(1, m + 1):
    matrix[0][j] = Gap * j

for i in range(n):
    for j in range(m):
        # match_score : check if the letter matches then assign match else mismatch
        if s1[i] == s2[j]:
            match_score = Match
        else:
            match_score = Mismatch

        # choosing the max value for each cell, thereby generating the matrix :
        matrix[i + 1][j + 1] = max(matrix[i][j] + match_score, matrix[i][j + 1] + Gap, matrix[i + 1][j] + Gap)

# Now, the matrix has been filled, we need to perform the traceback after identifying the score : There are more than
# one possibility of aligning, depending on the path taken during traceback, and we will be tracing all of those
# below :
aligned_seqs = []


def traceback(i, j, aligned_s1, aligned_s2):
    if i == 0 and j == 0:
        aligned_seqs.append((aligned_s1, aligned_s2))
        return
    if i > 0 and j > 0 and matrix[i][j] == matrix[i - 1][j - 1] + (Match if s1[i - 1] == s2[j - 1] else Mismatch):
        traceback(i - 1, j - 1, s1[i - 1] + aligned_s1, s2[j - 1] + aligned_s2)
    if i > 0 and matrix[i][j] == matrix[i - 1][j] + Gap:
        traceback(i - 1, j, s1[i - 1] + aligned_s1, '-' + aligned_s2)
    if j > 0 and matrix[i][j] == matrix[i][j - 1] + Gap:
        traceback(i, j - 1, '-' + aligned_s1, s2[j - 1] + aligned_s2)


traceback(n, m, '', '')
# Formatting and printing the output
print()
print('GLOBAL ALLIGNMENT MATRIX : ')
print()
matrix_print = matrix.transpose()
for i in range(m + 1):
    print(matrix_print[i])
print()
print("SCORE : ", matrix[n][m])
print()

print("Alignment(s) for the score generated above (optimal alignments): ")
for i, seq in enumerate(aligned_seqs):
    print(f'Alignment {i + 1}:')
    print(seq[0])
    print(seq[1])
    print('\n')
