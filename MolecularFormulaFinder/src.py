import numpy as np

from collections import Counter

N_ATOMS = 10
MAX_MASS_DIFF = 300.0

atoms = ("I", "Br", "Cl", "S", "P", "F", "O", "N", "C", "H")
atom2ind = {"I": 0, "Br": 1, "Cl": 2, "S": 3, "P": 4, "F": 5, "O": 6, "N": 7, "C": 8, "H": 9}
atom_masses = np.array([126.904468, 78.9183376, 34.96885271, 31.97207069, 30.97376151,
                        18.99840320, 15.9949146221, 14.0030740052, 12.0, 1.0078250321])


def decompose_formula(formula):
    res = np.zeros(N_ATOMS, dtype=int)
    i, j, k = 0, 0, 0
    while i < len(formula):
        j = i + 1
        while j < len(formula) and formula[j].islower():
            j += 1
        k = j
        while k < len(formula) and formula[k].isdigit():
            k += 1
        atom = formula[i:j]
        if j == k:
            count = 1
        else:
            count = int(formula[j:k])
        res[atom2ind[atom]] += count
        i = k
    return res


def formula2mass(formula):
    return np.dot(decompose_formula(formula), atom_masses)


class PNP:

    def __init__(self, formula, mass, multiedges, structure):
        self.formula = formula
        self.mass = mass
        self.multiedges = multiedges
        self.structure = structure

    def size(self):
        return self.structure.number_of_nodes()

    def atom_decomp(self, aa=-1):
        if aa == -1:
            formula = self.formula
        else:
            formula = self.structure.nodes[aa]['formula']
        return decompose_formula(formula)

    def aa_decomp(self):
        return Counter([self.structure.nodes[i]['formula'] for i in self.structure.nodes()])


class Node:

    def __init__(self, pred):
        self.childs = {}
        self.leaf = False
        self.pred = pred
        self.max_val = 0

    def __to_leaf(self):
        self.leaf = True
        self.vals = {}

    def insert(self, row):
        if row.shape[0] == 1:
            if not self.leaf:
                self.__to_leaf()
            self.vals[row[0]] = self.vals.get(row[0], 0) + 1
            self.max_val = max(self.max_val, self.vals[row[0]])
        else:
            if not row[0] in self.childs.keys():
                self.childs[row[0]] = Node(pred=self)
            self.childs[row[0]].insert(row[1:])
            self.max_val = max(self.max_val, self.childs[row[0]].max_val)

    def delete(self, row):
        if row.shape[0] == 1:
            self.vals[row[0]] -= 1
            self.max_val = max(self.vals.values())
        else:
            self.childs[row[0]].delete(row[1:])
            self.max_val = max([child.max_val for child in self.childs.values()])


class Distribution:

    def __init__(self, smooth_coef=10):
        self.smooth = smooth_coef
        self.cov = None
        self.__cov_inv = None
        self.__aprior_dist_coef = None
        self.bound = None
        self.data = None
        self.aa_list = None
        self.k = None
        self.top_k = None

    def fit(self, data, amino_acids):

        # aprior normal distribution params
        self.cov = np.mean(np.square(np.append(data, [[1] * N_ATOMS], axis=0)), axis=0)
        self.__cov_inv = 1 / self.cov
        self.__aprior_dist_coef = np.sum(np.log(self.cov)) + N_ATOMS * np.log(2 * np.pi)

        # discrete distribution from data
        min_bound = np.abs(np.amin(data, axis=0))
        max_bound = np.abs(np.amax(data, axis=0))
        self.bound = np.amax([min_bound, max_bound], axis=0)
        self.data = Node(pred=None)
        for row in data:
            self.data.insert(row)

        # list of amino_acids
        self.aa_list = [(formula2mass(formula), formula) for formula in amino_acids]
        self.aa_list.sort()

    def add_row(self, row):
        self.data.insert(row)

    def del_row(self, row):
        self.data.delete(row)

    def aprior_prob(self, x):
        return np.exp(-0.5 * (np.dot(np.square(x), self.__cov_inv) + self.__aprior_dist_coef))

    @staticmethod
    def next_ind(i):
        if i > 0:
            return -i
        else:
            return -i + 1

    def search(self, mass, ind_no_H, cur_node, is_node_real):

        if ind_no_H.shape[0] == N_ATOMS - 1:
            num_H = round((mass - np.dot(ind_no_H, atom_masses[:-1])) / atom_masses[-1])
            ind = np.append(ind_no_H, num_H)
            p = self.smooth * self.aprior_prob(ind)
            if is_node_real:
                p += cur_node.vals.get(num_H, 0)
            if abs(np.dot(ind, atom_masses) - mass) < 0.02 and p > self.top_k[self.k - 1][0]:
                self.top_k[self.k - 1] = (p, tuple(ind))
                self.top_k.sort(reverse=True)
        else:
            i = 0
            while i <= self.bound[ind_no_H.shape[0]]:
                best_estimated_ind = np.append(ind_no_H, [i] + [0] * (N_ATOMS - 1 - ind_no_H.shape[0]))
                if is_node_real and (i in cur_node.childs.keys()):
                    estimated_p = cur_node.childs[i].max_val + self.smooth * self.aprior_prob(best_estimated_ind)
                    if estimated_p >= self.top_k[self.k - 1][0]:
                        self.search(mass, np.append(ind_no_H, [i]), cur_node.childs[i], True)
                else:
                    estimated_p = self.smooth * self.aprior_prob(best_estimated_ind)
                    if estimated_p >= self.top_k[self.k - 1][0]:
                        self.search(mass, np.append(ind_no_H, [i]), None, False)
                i = Distribution.next_ind(i)

    def predict(self, mass, k=5):
        self.k = k
        self.top_k = [(0.0, tuple([0] * N_ATOMS))] * k
        self.search(mass, np.array([], dtype=int), self.data, True)
        indels = []
        mass = abs(mass)
        if len(self.aa_list) and self.aa_list[0][0] <= mass + 0.02 and mass - 0.02 <= self.aa_list[-1][0]:
            i_begin = 0
            i_end = len(self.aa_list)
            while i_end - i_begin > 1:
                i_mid = (i_begin + i_end) // 2
                if self.aa_list[i_mid][0] < mass - 0.02:
                    i_begin = i_mid
                else:
                    i_end = i_mid
            i1 = i_begin
            i_begin = 0
            i_end = len(self.aa_list)
            while i_end - i_begin > 1:
                i_mid = (i_begin + i_end) // 2
                if self.aa_list[i_mid][0] <= mass + 0.02:
                    i_begin = i_mid
                else:
                    i_end = i_mid
            i2 = i_end
            for i in range(i1, i2):
                if abs(self.aa_list[i][0] - mass) < 0.02:
                    indels.append(self.aa_list[i])
        return self.top_k, indels
