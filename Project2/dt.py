import os
import math
import random

path = './'
    
attributes = { 
        "Hydrophobic",
        "Polar",
        "Small",
        "Proline",
        "Tiny",
        "Aliphatic",
        "Aromatic",
        "Positive",
        "Negative",
        "Charged"
        }
residues = {
    "A" : { "Hydrophobic", "Small", "Tiny" },
    "C" : { "Hydrophobic", "Small" },
    "D" : { "Polar", "Small", "Negative", "Charged" },
    "E" : { "Polar", "Negative", "Charged" },
    "F" : { "Hydrophobic", "Aromatic" },
    "G" : { "Hydrophobic", "Small", "Tiny" },
    "H" : { "Polar", "Aromatic", "Positive", "Charged" },
    "I" : { "Hydrophobic", "Aliphatic" },
    "K" : { "Polar", "Positive", "Charged" },
    "L" : { "Hydrophobic", "Aliphatic" },
    "M" : { "Hydrophobic" },
    "N" : { "Polar", "Small" },
    "P" : { "Hydrophobic", "Small", "Proline" },
    "Q" : { "Polar" },
    "R" : { "Polar", "Positive", "Charged" },
    "S" : { "Polar", "Small", "Tiny" },
    "T" : { "Hydrophobic", "Polar", "Small" },
    "V" : { "Hydrophobic", "Small", "Aliphatic" },
    "W" : { "Hydrophobic", "Aromatic" },
    "Y" : { "Hydrophobic", "Polar", "Aromatic" },
}

class Question:
    def __init__(self, trait):
        self.trait = trait

    def match(self, record):
        return self.trait in residues[record]

    def __str__(self):
        return self.trait

class Node:
    def __init__(self, question, true_branch, false_branch):
        self.question = question
        self.true_branch = true_branch
        self.false_branch = false_branch

class Leaf:
    def __init__(self, records):
        E = 0
        B = 0
        for r in records:
            if r[1] == 'E':
                E += 1
            elif r[1] == 'B':
                B += 1
        if E > B:
            self.predictions = 'E'
        else:
            self.predictions = 'B'

def print_tree(node, indent=""):
    if isinstance(node, Leaf):
        print(indent + "Leaf prediction", node.predictions)
        return
    print(str(node.question))
    print(indent + "--> True:")
    print_tree(node.true_branch, indent + "  ")
    print(indent + "--> False:")
    print_tree(node.false_branch, indent + "  ")

def current_uncertainty(records):
    E = 0
    for r in records:
        if r[1] == 'E':
            E += 1
    return max(E, len(records) - E) / len(records)

def build_tree(records):
    max_gain, best_question = find_best_partition(records)
    # Reached the end of optimization
    if max_gain == 0:
        return Leaf(records)
    true_records, false_records = partition(records, best_question)
    true_branch = build_tree(true_records)
    false_branch = build_tree(false_records)
    return Node(best_question, true_branch, false_branch)

def shannon(left, right):
    p = float(len(left)) / (len(left) + len(right))
    return -1 * (p * math.log2(p) + (1 - p) * math.log2(1 - p))
    
def info_gain(left, right, current_uncertainty):
    print(current_uncertainty, shannon(left, right), current_uncertainty-shannon(left, right))
    return shannon(left, right)

def partition(records, question):
    true_records, false_records = [], []
    for record in records:
        if question.match(record[0]):
            true_records.append(record)
        else:
            false_records.append(record)
    return true_records, false_records

def find_best_partition(records):
    max_gain = 0
    best_question = None
    uncertainty = current_uncertainty(records)
    for attribute in attributes:
        question = Question(attribute)
        true_records, false_records = partition(records, question)
        if len(true_records) == 0 or len(false_records) == 0:
            continue
        gain = info_gain(true_records, false_records, uncertainty)
        if gain > max_gain:
            max_gain, best_question = gain, question
    return max_gain, best_question

def classify(record, node):
    if isinstance(node, Leaf):
        return node.predictions
    if node.question.match(record[0]):
        return classify(record, node.true_branch)
    else:
        return classify(record, node.false_branch)

if __name__ == '__main__':
    train = []
    test = []
    for filename in os.listdir(os.path.join(path, "fasta")):
        stem, _, fasta = filename.partition('.')
        aminoPath = os.path.join(path, "fasta", stem + '.fasta')
        sa = os.path.join(path, "sa", stem + '.sa')
        with open(aminoPath, 'r+') as oh:
            stem = oh.readline()[1:]
            amino = oh.readline().strip()
        with open(sa, 'r+') as oh:
            stem = oh.readline()[1:]
            EB = oh.readline().strip()
        assert(len(amino) == len(EB))
        ri = random.randint(0, 3)
        # 25%
        if ri == 0:
            test += list(zip(amino, EB))
        # 75%
        else:
            train += list(zip(amino, EB))
    tree = build_tree(train)
    E_realE = 0
    B_realE = 0
    E_realB = 0
    B_realB = 0
    for key, val in test:
        pred = classify(key, tree)
        if val == 'E':
            if pred == 'E':
                E_realE += 1
            elif pred == 'B':
                B_realE += 1
        elif val == 'B':
            if pred == 'E':
                E_realB += 1
            elif pred == 'B':
                B_realB += 1
    correct = E_realE + B_realB
    incorrect = E_realB + B_realE
    print("Accuracy", correct / (correct + incorrect))
    precision = E_realE / (E_realE + E_realB)
    print("Precision", precision)
    recall = E_realE / (E_realE + B_realE)
    print("Recall", recall)
    print("F1", 2 * (precision * recall) / (precision + recall))
    print_tree(tree)
    for key, val in residues.items():
        print(key, ": ", classify(key, tree))
