import os
import math
import random

def padMatrix(matrix):
    return [[-1] * 20] * 2 + matrix + [[-1] * 20] * 2

def extract200Features(matrix, i1, i2):
    features = []
    pi1 = i1 + 2 # RR INDICES START AT 1
    pi2 = i2 + 2
    for m in matrix[pi1-2:pi1+3]:
        features += m
    for m in matrix[pi2-2:pi2+3]:
        features += m
    assert len(features) == 200
    return features

def sigmoid(x):
    if x > 100:
        return 1
    if x < -100:
        return 0
    return 1 / (1 + math.exp(-x))

def cost(Y, pred):
    p = float(pred)
    assert p <= 1 and p >= 0
    cost = -1 * (Y * math.log(float(p)) + (1 - Y) * (math.log(1 - float(p))))
    return cost

def gradient(X, Y, cost):
    dw = scalarMul(cost - Y, X)
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

def dotProduct(a, b):
    return [i * j for i, j in zip(a, b)]

if __name__ == '__main__':
    path = './'
    predictions = []
    for filename in os.listdir(os.path.join(path, "rr")):
        if filename in {".", "..", ""}:
            continue
        stem, _, rr = filename.partition('.')
        rrPath = os.path.join(path, "rr", stem + ".rr")
        pssmPath = os.path.join(path, "pssm", stem + '.pssm')
        allContacts = []
        try:
            with open(rrPath, 'r+') as oh:
                seq = oh.readline()
                contacts = oh.readline().strip()
                while contacts:
                    i, j, d1, d2, d = contacts.split()[:5]
                    allContacts.append((int(i) - 1, int(j) - 1)) # RR index start 0
                    contacts = oh.readline().strip()
        except:
            continue
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
        contactFeatures = [(extract200Features(MATRIX, *contact), 'C', contact) for contact in allContacts]
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
            noncontactFeatures.append((extract200Features(MATRIX, *noncontact), 'N', noncontact))
        train_data = []
        test_data = []
        for c in contactFeatures + noncontactFeatures:
            ri = random.randint(0, 3)
            if ri == 0:
                test_data.append(c)
            else:
                train_data.append(c)
        learning_rate = 0.01
        eta = 0.25
        w, b = initializeWeights()
        iteration = 0
        while True:
            case = train_data[random.randint(0, len(train_data) - 1)]
            X = case[0]
            assert len(X) == 200
            if case[1] == 'C':
                Y = 1
            elif case[1] == 'N':
                Y = 0
            pred = (1 / 200) * (sum(list(map(sigmoid, [i * j for i, j in zip(X, w)]))) + sigmoid(b))
            error = cost(Y, pred)
            dw, db = gradient(X, Y, error)
            w = subtractLists(w, scalarMul(learning_rate, dw))
            b = b - learning_rate * db
            iteration += 1
            if iteration % 10000 == 0 or error < eta:
                print("Iteration: ", iteration)
                print("W b", w[:5], b)
                print("pred", pred)
                print("error", error)
                print("dw, db", sum(dw) / len(dw), db)
            if error < eta:
                break
        # Done with SGD
        resultsDict = {}
        for case in test_data:
            X = case[0]
            if case[1] == 'C':
                Y = 1
            elif case[1] == 'N':
                Y = 0
            i, j = case[2]
            pred = (1 / 200) * (sum(list(map(sigmoid, [s * t for s, t in zip(X, w)]))) + sigmoid(b))
            print(pred, Y)
            predictions.append((abs(Y - pred), j - i, L))
    print("Top L / 10: ", [x[0] for x in sorted(predictions)[::-1] if x[1] == x[2] // 10])
    print("Top L / 5: ", [x[0] for x in sorted(predictions)[::-1] if x[1] == x[2] // 5])
    print("Top L / 2: ", [x[0] for x in sorted(predictions)[::-1] if x[1] == x[2] // 2])
