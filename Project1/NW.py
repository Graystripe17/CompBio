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

