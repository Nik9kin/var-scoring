# -*- coding: utf-8 -*-

"""
    Molecular Formula Finder
"""

import time


def decomposition_printing():
    for ind in range(len(atoms)):
        if decomposition[ind] != 0:
            print(atoms[ind], decomposition[ind], sep="", end=" ")
    print()


def decompose_residue(remaining_mass: float, cur_atom: int):
    global num_of_solutions
    atom = atoms[cur_atom]
    if cur_atom == len(atoms) - 1:
        num_of_atom = round(remaining_mass / atomic_masses[atom])
        if abs(num_of_atom) <= coefficient_bounds[atom]:
            decomposition.append(num_of_atom)
            if abs(remaining_mass - atomic_masses[atom] * num_of_atom) < tolerance:
                decomposition_printing()
                num_of_solutions += 1
            decomposition.pop()
    else:
        for num_of_atom in range(-coefficient_bounds[atom], coefficient_bounds[atom] + 1):
            decomposition.append(num_of_atom)
            decompose_residue(remaining_mass - atomic_masses[atom] * num_of_atom, cur_atom + 1)
            decomposition.pop()


atoms = ("I", "Br", "Cl", "S", "P", "F", "O", "N", "C", "H")
atomic_masses = {"I": 126.904468, "Br": 78.9183376, "Cl": 34.96885271,
                 "S": 31.97207069, "P": 30.97376151, "F": 18.99840320,
                 "O": 15.9949146221, "N": 14.0030740052, "C": 12.0,
                 "H": 1.0078250321
                 }
coefficient_bounds = {"I": 0, "Br": 0, "Cl": 0, "S": 6, "P": 0, "F": 0, "O": 28, "N": 10, "C": 59, "H": 86}

print("Input total mass: ", end="")
total_mass = float(input())
tolerance = 0.02
decomposition = []
num_of_solutions = 0

start = time.time()
decompose_residue(total_mass, 0)
end = time.time()
print("Total number of solutions: ", num_of_solutions)
print("Program execution time: ", end - start, "s")
