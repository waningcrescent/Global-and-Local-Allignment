# Global-and-Local-Allignment
Dynamic Programming to find Global and Local allignments of two nucleotide sequences.

BIO213 – Introduction to Quantitative Biology
ASSIGNMENT-1
Shriya Verma 

1 Question. 1
There are two sequences given in the question : GATGCGCAG and GGCAGTA.
We are expected to perform global allignment on the given sequences with the
following scoring scheme:
Match = +2
Mismatch = -3
Gap = -1

# Can there be more than one optimal alignments?
It is possible to get more than one optimal alignment depending on the path
taken while tracing the matrix back. Multiple alignments can have the same
score. Therefore, it is important to list all the best possible alignments of a
given global alignment.


2 Question. 2
# Show Modification of the Results if the Scoring scheme is changed to:
Match = +2
Mismatch = -1
Gap = -3

2.1 Bidimensional Array obtained in the computation of
the Needleman-Wunsch Algorithm

The matrix obtained is different from the matrix obtained in Question 1, As
at iteration, the maximum value for a cell is selected. These values are heavily
dependant on the scoring scheme used. The values obtained by gap/matchmismatch
will change depending on the scoring scheme, and therefore, the maximum
value filled in each cell is also subject to change, implying that the matrix will
be different and so will the traceback.
Iteration: F(i, j)=max
{
F(i − 1, j) + gappenalty
F(i, j − 1) + gappenalty
F(i − 1, j − 1) + match (xi
, yj )

3
3 Question. 3
Now, we are required to perform local allignment of the same sequences :
GATGCGCAG and GGCAGTA using the scoring scheme defined below:
Match = +2
Mismatch = -1
Gap = -3

Local alignment is done using Smith-Waterman Algorithm, The initialization
step differs from the Needleman-Wunsch Algorithm here as the first row and
column are initialized with zeroes unlike the gap penalty used in NeedlemanWunsch Algorithm. The iteration step also differs as only non-negative values
are taken to fill the matrix :
Iteration: F(i, j)=max
{
0
F(i − 1, j) + gappenalty
F(i, j − 1) + gappenalty
F(i − 1, j − 1) + match (xi
, yj )
}


The score for a local alignment is the maximum value in the matrix and can
be present at any position.


4 Question. 4

# 4.1 Changes required in the program to perform local
rather than global pairwise sequence alignment
1. The first change was the initialization of the matrix:
For Global Alignment, The matrix is initialized by assigning gap penalties to all cells in the first row/column
For Local Alignment, The matrix is initialized by assigning value 0 to all
cells in the first row/column
2. Comparison with zero while filling matrix using :
max(0,value)

As in the Smith-Waterman Algorithm, the minimum value of any cell is 0.
3. The final change was the initialization of the traceback, While in global
alignment, the traceback always starts from the bottom-right-most cell.
In local allignment, the traceback starts from the cell having the maximum value. The code for local is required to iterate through the entire
matrix to find the maximum value in the matrix and it’s position before
calling the traceback function on it. I also had to stop the traceback as
soon as a cell value of ’0’ is reached, contrary to stopping on reaching
matrix[0][0].
4. Fundamentally the differences in Global and local alignment are :
Global: Initialization : F(i, 0) = gappenalty ∗ i
F(0, j) = gappenalty ∗ j
Iteration: F(i, j)=max
{
F(i − 1, j) + gappenalty
F(i, j − 1) + gappenalty
F(i − 1, j − 1) + match (xi
, yj ) }
Termination : Bottom Right
Local: Initialization : F(i, 0) = 0
F(0, j) = 0
Iteration: F(i, j)=max
{
0
F(i − 1, j) + gappenalty
F(i, j − 1) + gappenalty
F(i − 1, j − 1) + match (xi
, yj ) }
Termination : max. value in the matrix
