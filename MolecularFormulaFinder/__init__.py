import pickle

from .src import Distribution


with open('MolecularFormulaFinder/replacements.pickle', 'rb') as file:
    replacements = pickle.load(file)

with open('MolecularFormulaFinder/amino_acids.pickle', 'rb') as file:
    amino_acids = pickle.load(file)

mffinder = Distribution()
mffinder.fit(replacements, amino_acids)
