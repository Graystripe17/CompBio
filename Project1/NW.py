#!/usr/local/bin/python3

def alphabetize(l1, l2):
    return min(l1, l2), max(l1, l2)

# Load BLOSUM matrix
def loadBlosum(filename="blosum.txt"):
    axes = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']
    weightsDict = {}
    with open(filename, "r+") as oh:
        for line in oh.readlines():
            l1, _, weightStr = line.partition(' ')
            weights = [int(x) for x in weightStr.split()]
            for l2, w in zip(axes, weights):
                weightsDict[alphabetize(l1, l2)] = w
    assert(len(weightsDict) == 300)
    return weightsDict

weightsDict = loadBlosum()

def loadSequence(filename):
    with open(filename, "r+") as oh:
        first = oh.readline()[1:].strip()
        seq   = map(str.strip, oh.readlines())
    return ''.join(seq)

seq1 = loadSequence("1dzl_HPV.fasta")
seq2 = loadSequence("3mge_HIV.fasta")

import numpy as np
import seaborn as sns; sns.set()
def printMatrix(grid):
    print('\n'.join([''.join(['{:4}'.format(item) for item in row])
      for row in grid]))
    ax = sns.heatmap(np.array(grid))

seq1 = "GCATGCAAAAA"
seq2 = "GATTACAX"

seq1 = "EEEEEKKKKK"
seq2 = "EEEEE"

#seq1 = "AAAAAFFF"
#seq2 = "BBBBBFFF"
# Build empty grid
grid = [[0 for i in range(len(seq1) + 1)] for j in range(len(seq2) + 1)]
# Populate grid with gaps
try:
    for x in range(len(seq1) + 1):
        grid[0][x] = -4 * x
    for x in range(len(seq2) + 1):
        grid[x][0] = -4 * x
except:
    assert(False)

#   seq 1
# s
# e
# q
# 2

def getBestScore(i, j):
    # Get diag neighbor
    if j > 0 and i > 0:
        row = i - 1
        col = j - 1
        c = weightsDict[alphabetize(seq1[col-1], seq2[row-1])]
        d_score = grid[row][col] + c
    # Get left neighbor
    if j > 0:
        row = i
        col = j - 1
        # Gapping seq1
        c = weightsDict[alphabetize('*', seq2[row-1])]
        l_score = grid[row][col] + c
    # Get up neighbor
    if i > 0:
        row = i - 1
        col = j
        print("Row, Col, l1, l2", row, col, len(seq1), len(seq2))
        print(seq1, seq2)
        # Gapping seq2
        c = weightsDict[alphabetize(seq1[col-1], '*')]
        u_score = grid[row][col] + c
    return max(d_score, l_score, u_score)

maxScore = 0, 0, 0

for i in range(1, len(seq2) + 1):
    for j in range(1, len(seq1) + 1):
        bestScore = getBestScore(i, j)
        grid[i][j] = bestScore 
        if bestScore > maxScore[2]:
            maxScore = i, j, bestScore 

printMatrix(grid)
print(maxScore[2])

backtrackSteps = []

def backtrack(i, j):
    if i == 0 and j == 0:
        return
    # Get diag neighbor
    if j > 0 and i > 0:
        row = i - 1
        col = j - 1
        d_neighbor = grid[row][col]
    # Get left neighbor
    if j > 0:
        row = i
        col = j - 1
        l_neighbor = grid[row][col]
    # Get up neighbor
    if i > 0:
        row = i - 1
        col = j
        u_neighbor = grid[row][col]
    # As long as it's the largest, push onto the stack
    # Invariant: max backwards always converges
    bestNeighbor = max(d_neighbor, l_neighbor, u_neighbor)
    if bestNeighbor == d_neighbor:
        backtrackSteps.append('D')
        backtrack(i - 1, j - 1)
    elif bestNeighbor == l_neighbor:
        # seq1 has been gapped
        backtrackSteps.append('1')
        backtrack(i, j - 1)
    elif bestNeighbor == u_neighbor:
        # seq2 has been gapped
        backtrackSteps.append('2')
        backtrack(i - 1, j)

backtrack(maxScore[0], maxScore[1])
# backtrack(len(seq2), len(seq1))
print("BTSTEPS", backtrackSteps)
steps = backtrackSteps[::-1]
print("STEPS", steps)

def seqGenerator(seq):
    i = 0
    while i < len(seq):
        yield seq[i]
        i += 1

seq1Generator = seqGenerator(seq1)
seq2Generator = seqGenerator(seq2)

# Printing line 1
seq1List = []
midList  = []
seq2List = []
for s in steps:
    if s == 'D':
        f = next(seq1Generator)
        s = next(seq2Generator)
        seq1List.append(f)
        if f == s:
            midList.append('|')
        else:
            midList.append('*')
        seq2List.append(s)
    elif s == '1':
        seq1List.append('-')
        midList.append(' ')
        seq2List.append(next(seq2Generator))
    elif s == '2':
        seq1List.append(next(seq1Generator))
        midList.append(' ')
        seq2List.append('-')
# The max has been reached
try:
    while True:
        seq1List.append(next(seq1Generator))
        midList.append(' ')
        seq2List.append('-')
except:
    pass
try:
    while True:
        seq2List.append(next(seq2Generator))
        seq1List.append('-')
        midList.append(' ')
except:
    pass
print(seq1List, midList, seq2List)
print(''.join(seq1List))
print(''.join(midList))
print(''.join(seq2List))
