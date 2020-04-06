import os
import math
import random

def padMatrix(matrix):
    return [[-1] * 20] * 2 + matrix + [[-1] * 20] * 2

def extract200Features(matrix, i1, i2):
    features = []
    pi1 = i1 + 2 # RR INDICES START AT 1
    pi2 = i2 + 2
    print(i1, i2)
    print(pi1, pi2, len(matrix))
    for m in matrix[pi1-2:pi1+3]:
        print("M", m)
        features += m
    print("moving to pi2")
    for m in matrix[pi2-2:pi2+3]:
        print("M", m)
        features += m
    print(features, len(features))
    assert len(features) == 200
    return features

def sigmoid(x):
    return 1 / (1 + math.exp(-x))

def cost(Y, pred):
    cost = (Y * math.log(pred) + (1 - Y) * (math.log(1 - pred)))
    return cost

def gradient(X, Y, cost):
    dw = X * (cost - Y)
    db = cost - Y
    return dw, db

def initializeWeights():
    w = []
    for i in range(200):
        w.append(random.random())
    b = random.random()
    return w, b

def scalarMul(a, l):
    return [a * x for x in l]

def addLists(a, b):
    return [i + j for i, j in zip(a, b)]

def subtractLists(a, b):
    return [i - j for i, j in zip(a, b)]

if __name__ == '__main__':
    path = './'
    for filename in os.listdir(os.path.join(path, "rr")):
        train_data = []
        test_data = []
        stem, _, rr = filename.partition('.')
        rrPath = os.path.join(path, "rr", stem + ".rr")
        pssmPath = os.path.join(path, "pssm", stem + '.pssm')
        allContacts = []
        with open(rrPath, 'r+') as oh:
            seq = oh.readline()
            contacts = oh.readline().strip()
            while contacts:
                i, j, d1, d2, d = contacts.split()[:5]
                allContacts.append((int(i) - 1, int(j) - 1)) # RR index start 0
                contacts = oh.readline().strip()
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
                Y = seq[residueNumber]
                feature = [int(x) for x in pssm[3:23]]
                rawMatrix.append(feature)
                residueNumber += 1
                pssm = oh.readline().strip()
            MATRIX = padMatrix(rawMatrix)
        contactFeatures = [(extract200Features(MATRIX, *contact), 'C') for contact in allContacts]
        # Generate case Y = 0
        noncontactFeatures = []
        while len(noncontactFeatures) < len(contactFeatures):
            L = len(seq)
            r = random.randint(0, L - 2) # Inclusive to L
            s = random.randint(0, L - 2)
            a, b = min(r, s), max(r, s)
            while abs(a - b) < 5 or (a, b) in contactFeatures:
                b = random.randint(0, L - 2) # Reselect until farther away
            noncontact = (min(a, b), max(a, b))
            noncontactFeatures.append((extract200Features(MATRIX, *noncontact), 'N'))
        print("apples2oranges", contactFeatures[0], noncontactFeatures[0])
        for c in contactFeatures + noncontactFeatures:
            ri = random.randint(0, 3)
            if ri == 0:
                test_data += c
            else:
                train_data += c
        learning_rate = 0.01
        w, b = initializeWeights()
        print(len(train_data), train_data[0])
        for case in train_data:
            # Update weights
            X = case[0]
            print("X", X)
            assert len(X) == 200
            if case[1] == 'C':
                Y = 1
            elif case[1] == 'N':
                Y = 0
            pred = map(sigmoid, [i * j for i, j in zip(X, w)])
            dw, db = gradient(X, Y, cost(Y, pred))
            w = subtractLists(w, scalarMul(learning_rate, w))
            b = subtractLists(b, scalarMul(learning_rate, b))
            print("W B", w, b)
        for case in test_data:
            X = case[0]
            if case[1] == 'C':
                Y = 1
            elif case[1] == 'N':
                Y = 0
            pred = map(sigmoid, [i * j for i, j in zip(X, w)])
            print(pred)

