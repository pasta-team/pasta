# coding: utf-8
# Copyright © 2016 YunXing Zuo, WeiJi Hsiao

from math import fabs
import numpy as np
import itertools
import random
from collections import OrderedDict
from pasta.element import Element
import spglib
import forstr


class Structure(object):

    def __init__(self, lattice, positions, atoms, is_cartesian=False, name=None):
        self.__lattice = np.array(lattice).reshape((3, 3))
        if is_cartesian:
            self.__cart_positions = np.array(positions)
            self.__frac_positions = np.dot(
                self.__cart_positions, np.linalg.inv(self.__lattice))
        else:
            self.__frac_positions = np.array(positions)
            self.__cart_positions = np.dot(
                self.__frac_positions, self.__lattice)
        self.__atoms = []
        for atom in atoms:
            self.__atoms.append(atom)

        self.__species = OrderedDict()
        for species, atoms in itertools.groupby(self.__atoms):
            self.__species[species] = len(list(atoms))
        self.__rdf = OrderedDict()
        if name is None:
            name = ' '.join('{:s}{:d}'.format(
                element, self.__species[element]) for element in self.__species)
        self.__name = name

    def __repr__(self):
        strings = []
        strings.append('Structure of {name}'.format(name=self.__name))
        strings.append(str(self.__lattice))
        strings.append('Atoms')
        for element, num in zip(self.element_species, self.num_per_species):
            strings.append('{:>4s} : {:<4d}'.format(element, num))
        strings.append('Fractional Coordinates')
        for atom, position in zip(self.__atoms, self.__frac_positions):
            strings.append(
                '{:>4s} : {:>16.8f}{:>16.8f}{:>16.8f}'.format(atom, *position))
        return '\n'.join(strings)

    @staticmethod
    def import_from_vasp(filename):
        lattice = []
        positions = []
        atoms = []
        with open(filename) as f:
            name = f.readline().strip()
            scale = float(f.readline())
            for i in range(3):
                lattice += [float(coord) *
                            scale for coord in f.readline().split()]
            species = f.readline().split()
            number = [int(i) for i in f.readline().split()]
            for ele, num in zip(species, number):
                for i in range(num):
                    atoms.append(ele)
            line = f.readline()
            if line[0] in ('S', 's'):
                line = f.readline()
            if line[0] in ('C', 'c', 'K', 'k'):
                is_cartesian = True
            else:
                is_cartesian = False
            while True:
                line = f.readline()
                if not line:
                    break
                positions.append([float(i) for i in line.split()[0:3]])
        return Structure(lattice, positions, atoms, is_cartesian, name)

    @staticmethod
    def import_from_pwmat(filename):
        lattice = []
        positions = []
        atoms = []
        with open(filename) as f:
            line_num = 0
            while True:
                line = f.readline()
                if not line:
                    break
                line_num += 1
                if line_num in (3, 4, 5):
                    lattice += [float(i) for i in line.split()]
                elif line_num > 6:
                    if len(line.split()) == 1:
                        break
                    data = [float(i) for i in line.split()]
                    atoms.append(Element(int(data[0])).atomic_symbol)
                    positions.append(data[1:4])
        return Structure(lattice, positions, atoms)

    def export_to_vasp(self, filename, is_cartesian=False, scale=1.0):
        lines = []
        lines.append(self.__name)
        lines.append(str(scale))
        lattice = self.__lattice / scale
        for i in (0, 1, 2):
            lines.append('{:>16.8f}{:>16.8f}{:>16.8f}'.format(
                lattice[i][0], lattice[i][1], lattice[i][2]))
        lines.append(''.join('{:>6s}'.format(element)
                             for element in self.element_species))
        lines.append(''.join('{:>6d}'.format(num)
                             for num in self.num_per_species))
        if is_cartesian:
            positions = self.cart_positions / scale
            lines.append('Cartesian')
        else:
            positions = self.frac_positions
            lines.append('Direct')
        for position in positions:
            lines.append('{:>16.8f}{:>16.8f}{:>16.8f}'.format(*position))

        with open(filename, 'w') as f:
            f.write('\n'.join(lines))

    def export_to_pwmat(self, filename):
        with open(filename, 'w') as f:
            f.write('{:d}\n'.format(self.number_of_atoms))
            f.write('Lattice vector\n')
            for i in (0, 1, 2):
                f.write('{:>16.8f}{:>16.8f}{:>16.8f}\n'.format(
                    self.lattice[i][0], self.lattice[i][1], self.lattice[i][2]))
            f.write('Position\n')
            atnum = []
            for ele in self.__species:
                num = Element(ele).atomic_number
                for i in range(self.__species[ele]):
                    atnum.append(num)
            for (number, position) in zip(atnum, self.frac_positions):
                f.write('{:>4d}{:>16.8f}{:>16.8f}{:>16.8f}{:>6d}{:>6d}{:>6d}\n'.format(
                    number, position[0], position[1], position[2], 1, 1, 1))

    def as_tuple(self):
        return (self.lattice, self.frac_positions, np.array([Element(atom).atomic_number for atom in self.__atoms]))

    def multiple(self, trans_matrix):
        lattice = self.lattice.copy()
        positions = self.cart_positions.copy()
        pc = positions.copy()
        n1 = int(trans_matrix[0])
        n2 = int(trans_matrix[1])
        n3 = int(trans_matrix[2])
        atoms = self.atoms
        for i in range(n1):
            for j in range(n2):
                for k in range(n3):
                    if i + j + k > 0:
                        inc = i * lattice[0] + j * lattice[1] + k * lattice[2]
                        p = pc + inc
                        positions = np.vstack((positions, p))
        atoms *= n1 * n2 * n3
        return Structure(lattice, positions, atoms, is_cartesian=True)

    def get_crystal_system(self, symprec=1e-1, angle_tolerance=5):
        n = spglib.get_symmetry_dataset(
            self.as_tuple(), symprec=symprec, angle_tolerance=angle_tolerance)['number']

        def f(i, j): return n - i + 1 if i <= n <= j else 0
        cs = {"triclinic": (1, 2), "monoclinic": (3, 15), "orthorhombic": (16, 74),
              "tetragonal": (75, 142), "trigonal": (143, 167), "hexagonal": (168, 194),
              "cubic": (195, 230)}
        crystal_system = None
        for k, v in cs.items():
            if f(*v):
                crystal_system = (k, f(*v))
                break
        return crystal_system

    def find_rdf(self, max_R=20, step=0.2,):
        volume = np.dot(self.__lattice[0], np.cross(
            self.__lattice[1], self.__lattice[2]))
        q = [0] * 3
        for i in (0, 1, 2):
            q[i] = int(max_R * np.linalg.norm(np.cross(
                self.lattice[(i + 2) % 3], self.lattice[(i + 1) % 3])) / volume + 1)

        n = len(self.__species)
        nr = int(max_R // step + 1)
        rdf = np.zeros((int(n * (n + 1) / 2), nr))
        p = self.cart_positions
        l = np.array(list(self.num_per_species))
        lat = self.lattice

        forstr.mul(rdf, p, l, n, int(self.number_of_atoms),
                   nr, q[0], q[1], q[2], lat, max_R, step)

        k = 0
        lsp = list(self.element_species)
        llsp = len(lsp)
        for i in range(llsp):
            for j in range(i, llsp):
                self.__rdf[(lsp[i], lsp[j])] = rdf[k].copy()
                k += 1

    def rdf_for_np(self, uplimit=4):
        if len(list(self.element_species)) > uplimit:
            raise ValueError
        row = uplimit * (uplimit - 1) / 2 + uplimit
        res = np.zeros((row, 102))
        li = list(range(row))
        random.shuffle(li)
        for num, i in enumerate(self.rdf):
            res[li[num]] = [elements_dict[i[0]],
                            elements_dict[i[1]]] + self.rdf[i]
        return res

    def write_rdf(self, directory):
        for i in self.rdf.keys():
            f = open('{:s}/{:s}-{:s}_zhu'.format(directory, i[0], i[1]), 'w')
            for j in range(len(self.rdf[i])):
                f.write('{:f}\t{:f}\n'.format(0.02 * (j + 1), self.rdf[i][j]))
            f.close()

    def dif(self, s, max_R=20, step=0.02):
        self.find_rdf()
        s.find_rdf()
        rdf1 = self.rdf
        rdf2 = s.rdf
        value = 0
        length = int(max_R // step + 1)
        for i in rdf1:
            datom = 0
            atom1 = 0
            atom2 = 0
            for j in range(length):
                atom1 += rdf1[i][j] * step * step * step * length * length
                atom2 += rdf2[i][j] * step * step * step * length * length
                datom += fabs(atom1 - atom2) * step
            if atom1 < atom2:
                q = datom / atom1
                if value < q:
                    value = q
            else:
                q = datom / atom2
                if value < q:
                    value = q
        return value

    def create_supercell(self, trans_matrix):
        lat = np.dot(trans_matrix, self.lattice)
        s = self.multiple((fabs(np.linalg.det(trans_matrix)),
                           fabs(np.linalg.det(trans_matrix)), 1))
        frac = np.dot(s.cart_positions, np.linalg.inv(lat))
        for i in range(len(frac)):
            for j in range(3):
                while frac[i][j] < 0:
                    frac[i][j] += 1
                while frac[i][j] >= 1:
                    frac[i][j] -= 1
        return Structure(lat, frac, s.atoms)

    @property
    def rdf(self):
        return self.__rdf

    @property
    def atoms(self):
        return [atom for atom in self.__atoms]

    @property
    def lattice(self):
        return self.__lattice

    @property
    def frac_positions(self):
        return self.__frac_positions

    @property
    def cart_positions(self):
        return self.__cart_positions

    @property
    def number_of_atoms(self):
        return len(self.__atoms)

    @property
    def element_species(self):
        return self.__species.keys()

    @property
    def num_per_species(self):
        return self.__species.values()

    @property
    def species(self):
        return self.__species

    @property
    def min_bond_length(self):
        pos = self.cart_positions.copy()
        pc = pos.copy()
        n = 8 * self.number_of_atoms
        for i in (0, 1):
            for j in (0, 1):
                for k in (0, 1):
                    if i + j + k > 0:
                        inc = i * self.lattice[0] + j * \
                            self.lattice[1] + k * self.lattice[2]
                        p = pc + inc
                        pos = np.vstack((pos, p))

        return forstr.min_bond_length(pos, n)

    @property
    def volume(self):
        return fabs(np.dot(self.lattice[0], np.cross(self.lattice[1], self.lattice[2])))

    @property
    def atom_density(self):
        return self.number_of_atoms / self.volume

    @property
    def rdf_for_lx(self):
        li = self.rdf[(self.atoms[0], self.atoms[0])][:3]
        for i in range(3):
            li[i] *= self.min_bond_length**3 * (i + 1)**2 * 1.02**3
            li[i] = int(li[i])
        return li

    @property
    def name(self):
        return self.__name
