import pickle
from src import Distribution


with open('replacements.pickle', 'rb') as file:
    replacements = pickle.load(file)

with open('amino_acids.pickle', 'rb') as file:
    amino_acids = pickle.load(file)

MFFinder = Distribution()
MFFinder.fit(replacements, amino_acids)
