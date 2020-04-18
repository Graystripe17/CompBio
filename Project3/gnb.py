import os
import math
import random
import pickle

def padMatrix(matrix):
    return [([-1] * 20, '*')] * 2 + matrix + [([-1] * 20, '*')] * 2

# Returns P(xi|classification)
def gaussian(x, mean, sigma):
    if sigma == 0:
        if x == mean:
            return 1.0
        else:
            return 0.0
    exp = math.exp(-1 * ((x - mean) ** 2) / (2 * sigma ** 2))
    flat = 1 / (2 * math.pi * sigma ** 2) ** 0.5
    probability = exp * flat
    return probability

def processFeatures(matrix):
    processedMatrix = []
    for i in range(2, len(matrix) - 2):
        Y = matrix[i][-1]
        features = []
        for m in matrix[i-2:i+3]:
            features += m[0]
        processedMatrix.append((features, Y))
    assert len(matrix) == len(processedMatrix) + 4
    return processedMatrix

def avg(x):
    return sum(x) / len(x)

def stddev(X):
    mean = avg(X)
    variance = sum([((x - mean) ** 2) for x in X]) / len(X)
    stddev = variance ** 0.5
    return stddev

def populateModel(posterior, targetFeatures, allFeatures):
    M = []
    S = []
    for i in range(100):
        values = [scan[i] for scan in targetFeatures]
        M.append(avg(values))
        S.append(stddev(values))
    P = len(targetFeatures) / len(allFeatures)
    return M, S, P 

def train(records):
    model = {}
    allFeatures = [record[1] for record in records]
    for posterior in {'H', 'E', 'C'}:
        targetFeatures = [record[0] for record in records if record[1] == posterior]
        model[posterior] = populateModel(posterior, targetFeatures, allFeatures)
    return model

def classify(record, model):
    results = {}
    features = record[0]
    for posterior in {'H', 'E', 'C'}:
        M, S, P = model[posterior] 
        featureProduct = 1
        for i in range(100):
            featureProduct *= gaussian(features[i], M[i], S[i])
        results[posterior] = featureProduct * P
    return results

def selectWinner(record, model):
    results = classify(record, model)
    winner = max(results, key=results.get)
    return winner
    

def test(records, model):
    correct, incorrect = 0, 0
    for record in records:
        actual = record[1]
        predict = selectWinner(record, model)
        print(actual, predict)
        if actual == predict:
            correct += 1
        else:
            incorrect += 1
    #print(correct / (correct + incorrect))
    print(predict)

def getMatrix(ss, pssmPath):
    with open(pssmPath, 'r+') as oh:
        _ = oh.readline() # blank line
        _ = oh.readline() # instruction
        labels = oh.readline().split()[:20]
        rawMatrix = []
        pssm = oh.readline().strip()
        residueNumber = 0
        while pssm:
            pssm = pssm.split()
            residue = pssm[2]
            Y = ss[residueNumber]
            feature = [int(x) for x in pssm[3:23]]
            rawMatrix.append((feature, Y))
            residueNumber += 1
            pssm = oh.readline().strip()
        MATRIX = padMatrix(rawMatrix)
    return MATRIX
    
def getRecordsByProtein(stem):
    fastaPath = os.path.join(path, "fasta", stem + ".fasta")
    ssPath = os.path.join(path, "ss", stem + ".ss")
    pssmPath = os.path.join(path, "pssm", stem + '.pssm')
    with open(fastaPath, 'r+') as oh:
        stem = oh.readline()[1:]
        fasta = oh.readline().strip()
    with open(ssPath, 'r+') as oh:
        stem = oh.readline()[1:]
        ss = oh.readline().strip()
    assert(len(fasta) == len(ss))
    MATRIX = getMatrix(ss, pssmPath)
    records = processFeatures(MATRIX)
    return records

def trainModel():
    train_data = []
    test_data = []
    path = './'
    for filename in os.listdir(os.path.join(path, "fasta")):
        stem, _, fasta = filename.partition('.')
        records = getRecordsByProtein(stem)
        ri = random.randint(0, 3)
        if ri == 0:
            test_data += records
        else:
            train_data += records
    model = train(train_data)
    print("Training complete")
    return model

def getHEC(stem):
    records = getRecordsByProtein(stem)
    model = trainModel()
    return classify(records, model)

def writePickle(hecDict, filename):
    with open(filename, 'wb') as oh:
        pickle.dump(hecDict, oh)

if __name__ == '__main__':
    train_data = []
    test_data = []
    path = './'
    for filename in os.listdir(os.path.join(path, "fasta")):
        stem, _, fasta = filename.partition('.')
        records = getRecordsByProtein(stem)
        ri = random.randint(0, 3)
        if ri == 0:
            test_data += records
        else:
            train_data += records
    model = train(train_data)
    test(test_data, model)
    input("Part 1 complete")
    hecDict = {}
    for filename in os.listdir(os.path.join(path, "fasta")):
        stem, _, fasta = filename.partition('.')
        hecCount = {'H': 0, 'E': 0, 'C': 0}
        records = getRecordsByProtein(stem)
        for record in records:
            winner = selectWinner(record, model)
            hecCount[winner] += 1
        hecTotal = hecCount['H'] + hecCount['E'] + hecCount['C']
        hecH = hecCount['H'] / hecTotal
        hecE = hecCount['E'] / hecTotal
        hecC = hecCount['C'] / hecTotal
        hecDict[stem] = hecH, hecE, hecC
    pickleName = "HEC.pickle"
    writePickle(hecDict, pickleName)
    print(hecDict)
    print("Pickle written to", pickleName)
