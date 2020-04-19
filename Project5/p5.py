import os
import math
import pickle
import random

def dotProd(x1, x2):
    return sum([a * b for a, b in zip(x1, x2)])

def predict(B, B_0, features):
    return B_0 + dotProd(B, features) 

def getTM(p1, p2):
    tmPath = os.path.join("tmalign", p1 + '_' + p2 + '_' + 'tmalign')
    with open(tmPath, 'r+') as oh:
        for i in range(17):
            oh.readline()
        TM1 = float(oh.readline().split()[1])
        TM2 = float(oh.readline().split()[1])
        return (TM1 + TM2) / 2

if __name__ == '__main__':
    path = './'
    percentagesDict = {} # 20 features per protein
    hecDict = {} # 3 features per protein
    ebDict = {} # 2 features per protein
    for filename in os.listdir(os.path.join(path, "fasta")):
        stem, _, fasta = filename.partition('.')
        fastaPath = os.path.join(path, "fasta", stem + ".fasta")
        pssmPath = os.path.join(path, "pssm", stem + ".pssm")
        with open(fastaPath, 'r+') as oh:
            stem = oh.readline()[1:].strip()
            fasta = oh.readline().strip()
        with open(pssmPath, 'r+') as oh:
            _ = oh.readline() # blank line
            _ = oh.readline() # instruction
            labels = oh.readline().split()[20:]
            assert(len(labels) == 20)
            percentages = []
            line = oh.readline().strip()
            while line:
                pline = line.split()[22:42]
                percentages.append(pline)
                line = oh.readline().strip()

            # Calc column percent avg
            meanp = []
            for i in range(20):
                total = 0
                for line in percentages:
                    total += int(line[i])
                final = total / len(percentages) / 100.0
                meanp.append(final)
            print(meanp)
            percentagesDict[stem] = meanp

    with open("hec.pickle", 'rb') as oh:
        hecDict = pickle.load(oh)
    
    with open("eb.pickle", 'rb') as oh:
        ebDict = pickle.load(oh)

    allProteins = [p.split('.')[0] for p in os.listdir(os.path.join(path, "fasta"))]
    train = []
    test = []
    for i in range(len(allProteins) - 1):
        for j in range(i + 1, len(allProteins)):
            p1 = allProteins[i]
            p2 = allProteins[j]
            TM = getTM(p1, p2)
            features = []
            # Feature order: 40 percents, 6 hec, 4 eb
            hecList = hecDict[p1] + hecDict[p2]
            assert len(hecList) == 6
            ebList = [ebDict[p1]['E'], ebDict[p1]['B'], ebDict[p2]['E'], ebDict[p2]['B']]
            features = tuple(percentagesDict[p1]) + tuple(percentagesDict[p2]) + (hecList) + tuple(ebList)
            assert len(features) == 50
            ri = random.randint(0, 3)
            if ri == 0:
                test.append((TM, p1, p2, features))
            else:
                train.append((TM, p1, p2, features))
    # Initialize weights
    B = [0 for _ in range(50)]
    B_0 = random.random()
    alpha = 0.001
    epsilon = 0.01
    iterations = 0
    while True:
        TM, p1, p2, features = random.choice(train)
        pred = predict(B, B_0, features)
        squared_error = (TM - pred) ** 2
        partial_wrt_b0 = (1 / 50) * ((TM - (B_0 + dotProd(B, features))))
        partial_wrt_b = [(1 / 50) * (features[i]) * (TM - B_0 + dotProd(B, features)) for i in range(50)]
        # Update weights
        B_0 = B_0 + alpha * partial_wrt_b0
        B = [B[i] + alpha * partial_wrt_b[i] for i in range(50)]
        print(str(squared_error) + '\t' + str(pred))
        if squared_error < epsilon:
            break
        iterations += 1
        if iterations % 10000 == 0:
            print("FEATURES", features)
            print("B", B, B_0)
            print("pred", pred)


