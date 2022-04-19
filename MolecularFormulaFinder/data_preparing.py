import pickle

import networkx as nx

from src import MAX_MASS_DIFF, decompose_formula, PNP


LIB_SIZE = 5021
path_pnpdb = '../../side_sources/pnpdatabase/'

pnps = []
replacements = []
amino_acids = set()

with open(path_pnpdb + 'pnp_AAGraphs.txt', 'r') as file:
    for ind in range(LIB_SIZE):
        _ = file.readline()

        _, _, formula, _, _, mass = file.readline().split()
        mass = float(mass)

        _ = file.readline()
        *_, num_aa = file.readline().split()
        num_aa = int(num_aa)

        structure = nx.MultiDiGraph()

        for i in range(num_aa):
            AA_id, AA_formula, AA_mass = file.readline().split()
            AA_id = int(AA_id)
            AA_mass = float(AA_mass)
            structure.add_node(AA_id, formula=AA_formula, mass=AA_mass)

        *_, num_edges = file.readline().split()
        num_edges = int(num_edges)

        multiedges = False

        for i in range(num_edges):
            v1, bound, v2 = file.readline().split()
            assert bound == '-NC>'

            v1, v2 = int(v1), int(v2)
            if (v1, v2) in structure.edges:
                multiedges = True
            structure.add_edge(v1, v2)

        pnp_type = file.readline()[:-1]

        pnps.append(PNP(formula, mass, multiedges, structure))

for ind1 in range(LIB_SIZE):
    amino_acids.update(list(pnps[ind1].aa_decomp().keys()))
    for ind2 in range(LIB_SIZE):
        p1 = pnps[ind1]
        p2 = pnps[ind2]
        if ind1 != ind2 and abs(p1.mass - p2.mass) < MAX_MASS_DIFF:
            AA_set1 = p1.aa_decomp()
            AA_set2 = p2.aa_decomp()
            gone = AA_set1 - AA_set2
            gain = AA_set2 - AA_set1
            if (sum(gone.values()) == 1) and (sum(gain.values()) == 1):
                diff = decompose_formula(list(gone)[0]) - decompose_formula(list(gain)[0])
                replacements.append(diff)

with open('replacements.pickle', 'wb') as file:
    pickle.dump(replacements, file)

with open('amino_acids.pickle', 'wb') as file:
    pickle.dump(amino_acids, file)
