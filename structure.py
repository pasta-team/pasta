# coding: utf-8
# Copyright Â© 2016 YunXing Zuo, WeiJi Hsiao

from math import fabs,sqrt
import numpy as np
from scipy.spatial import Voronoi
import itertools
import ctypes
import random
from collections import OrderedDict
from numpy.ctypeslib import ndpointer
from pasta.element import Element
import spglib
import re
lib = ctypes.cdll.LoadLibrary('/home/jjs/zhupylab/pasta/forstr.so')

def gcd(a,b):
    while b>0:
        a,b = b, a%b
    return a

vesta_dic = {}
with open('/home/jjs/machinel/hse_gap/temp/bond_vesta') as f:
    while True:
        line = f.readline()
        if not line:
            break
        data = line.split()
        vesta_dic[tuple(sorted(data[:2]))] = float(data[-1])

class Structure(object):

    def __init__(self, lattice, positions, atoms, is_cartesian=False, name=None,):
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

        #atoms.sort()
        self.__species = OrderedDict()
        self.__sorted = True
        for species, at_species in itertools.groupby(atoms):
            if species not in self.__species:
                self.__species[species] = len(list(at_species))
            else:
                self.__species[species] += len(list(at_species))
                self.__sorted = False
        self.__rdf = OrderedDict()
        if name is None:
            name = ' '.join('{:s}{:d}'.format(element, self.__species[element]) for element in self.__species)
        self.__name = name
        self.__lengths = [np.linalg.norm(i) for i in self.__lattice]
        self.__angles = np.zeros(3)
        for i in range(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            cos_angle = np.dot(self.__lattice[j], self.__lattice[k]) / (np.linalg.norm(self.__lattice[j]) * np.linalg.norm(self.__lattice[k]))
            self.__angles[i] = np.arccos(cos_angle) * 180.0 / np.pi

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
                lattice += [float(coord) * scale for coord in f.readline().split()]
            species = f.readline().split()
            number = [int(i) for i in f.readline().split()]
            for ele,num in zip(species,number):
                for i in range(num):
                    atoms.append(ele)
            line = f.readline()
            if line[0] in ('S','s'):
                line = f.readline()
            if line[0] in ('C','c','K','k'):
                is_cartesian = True
            else:
                is_cartesian = False
            for i in range(len(atoms)):
                line = f.readline()
                positions.append([float(i) for i in line.split()[0:3]])
        return Structure(lattice, positions, atoms, is_cartesian, name)
    
    @staticmethod
    def import_from_pwmat(filename):
        lattice = []
        positions = []
        atoms = []
        with open(filename) as f:
            line = f.readline()
            num_atoms = int(line.split()[0])
            f.readline()
            for i in range(3):
                lattice.append([float(coord) for coord in f.readline().split()])
            f.readline()
            for i in range(num_atoms):
                data = f.readline().split()
                atoms.append(Element(int(data[0])).atomic_symbol)
                positions.append(list(map(float,data[1:4])))
        return Structure(lattice,positions,atoms)
    
    def export_to_vasp(self, filename, is_cartesian=False, scale=1.0):
        if not self.__sorted:
            self.sort_()
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
        if not self.__sorted:
            self.sort_()
        with open(filename, 'w') as f:
            f.write('{:d}\n'.format(self.number_of_atoms))
            f.write('Lattice vector\n')
            for i in (0,1,2):
                f.write('{:>16.8f}{:>16.8f}{:>16.8f}\n'.format(self.lattice[i][0],self.lattice[i][1],self.lattice[i][2]))
            f.write('Position\n')
            for (ele,position) in zip (self.atoms,self.frac_positions):
                number = Element(ele).atomic_number
                f.write('{:>4d}{:>16.8f}{:>16.8f}{:>16.8f}{:>6d}{:>6d}{:>6d}\n'.format(number,position[0],position[1],position[2],1,1,1))

    def sort_(self):
        new_fracs = []
        new_atoms = []
        for i in self.species:
            for j in range(self.number_of_atoms):
                if self.atoms[j] == i:
                    new_atoms.append(self.atoms[j])
                    new_fracs.append(self.frac_positions[j])
        self.__frac_positions = np.array(new_fracs)
        self.__atoms = new_atoms
        self.__cart_positions = np.dot(self.__frac_positions, self.__lattice)
        self.__sorted = True
    
    def as_tuple(self):
        return (self.lattice, self.frac_positions,np.array([Element(atom).atomic_number for atom in self.__atoms]))
    
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
        lattice[0] *= n1
        lattice[1] *= n2
        lattice[2] *= n3
        return Structure(lattice, positions, atoms, is_cartesian=True)
    
    def get_crystal_system(self, symprec = 1e-1, angle_tolerance = 5):
        n = spglib.get_symmetry_dataset(self.as_tuple(), symprec = symprec, angle_tolerance = angle_tolerance)['number']
        f = lambda i, j: n-i+1 if i <= n <= j else 0
        cs = {"triclinic": (1, 2), "monoclinic": (3, 15), "orthorhombic": (16, 74), \
        "tetragonal": (75, 142), "trigonal": (143, 167), "hexagonal": (168, 194), \
        "cubic": (195, 230)}
        crystal_system = None
        for k, v in cs.items():
            if f(*v):
                crystal_system = (k,f(*v))
                break
        return crystal_system
    
    @property
    def spg_number(self, symprec = 1e-1, angle_tolerance = 5):
        return spglib.get_symmetry_dataset(self.as_tuple(), symprec = symprec, angle_tolerance = angle_tolerance)['number']

    def find_conventional_cell(self, symprec = 1e-4):
        lattice, frac_position, ele_num = spglib.refine_cell(self.as_tuple(), symprec = symprec)
        atoms = [Element(int(i)).atomic_symbol for i in ele_num]
        return Structure(lattice, frac_position, atoms)


    def find_rdf(self, max_R=15, step=0.15,):
        volume = np.dot(self.__lattice[0], np.cross(
            self.__lattice[1], self.__lattice[2]))
        q = [0] * 3
        for i in (0, 1, 2):
            q[i] = int(max_R * np.linalg.norm(np.cross(
                self.lattice[(i + 2) % 3], self.lattice[(i + 1) % 3])) / volume + 1)
        # print(q)

        n = len(self.__species)
        nr = int(max_R // step + 1)
        rdf = np.zeros((int(n * (n + 1) / 2), nr))
        rpp = (rdf.__array_interface__[
               'data'][0] + np.arange(rdf.shape[0]) * rdf.strides[0]).astype(np.uintp)
        p = self.cart_positions
        ppp = (p.__array_interface__[
               'data'][0] + np.arange(p.shape[0]) * p.strides[0]).astype(np.uintp)
        l = list(self.num_per_species)
        nps = (ctypes.c_int * n)(*l)
        lat = self.lattice
        lpp = (lat.__array_interface__[
               'data'][0] + np.arange(lat.shape[0]) * lat.strides[0]).astype(np.uintp)

        _doublepp = ndpointer(dtype=np.uintp, ndim=1, flags='C')
        lib.mul.argtypes = [_doublepp, _doublepp, ctypes.POINTER(
            ctypes.c_int), ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, _doublepp, ctypes.c_double, ctypes.c_double]

        lib.mul(rpp, ppp, nps, n, self.number_of_atoms,
                nr, q[0], q[1], q[2], lpp, max_R, step)

        k = 0
        lsp = list(self.element_species)
        llsp = len(lsp)
        for i in range(llsp):
            for j in range(i, llsp):
                self.__rdf[(lsp[i], lsp[j])] = rdf[k].copy()
                k += 1

    def rdf_for_np(self, uplimit=4, sample=6):
        self.find_rdf()
        if len(list(self.element_species)) > uplimit:
            raise ValueError
        row = int(uplimit * (uplimit - 1) / 2 + uplimit)
        res = []
        for i in range(sample):
            m = np.zeros((row, 106))
            li = list(range(row))
            random.shuffle(li)
            for num, i in enumerate(self.rdf):
                ele1 = Element(i[0])
                ele2 = Element(i[1])
                if ele1.electronegativity is None or ele2.electronegativity is None:
                    raise ValueError('electronegativity not exist')
                li2 = list(ele1.atomic_coordination) + [ele1.electronegativity] + list(
                    ele2.atomic_coordination) + [ele2.electronegativity] + list(self.rdf[i])
                m[li[num]] = np.array(li2)
            res.extend(list(m.flatten()))
        return res

    def write_rdf(self, directory):
        for i in self.rdf.keys():
            f = open('{:s}/{:s}-{:s}_zhu'.format(directory, i[0], i[1]), 'w')
            for j in range(len(self.rdf[i])):
                f.write('{:f}\t{:f}\n'.format(0.02 * (j + 1), self.rdf[i][j]))
            f.close()

    def dif(self, s, max_R=15, step=0.15):
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
    
    def pair_bond_length(self, element1, element2, filename=None):
        s = self.multiple((3,3,3))
        bond = []
        for i,p,f in zip(self.atoms,self.cart_positions,self.frac_positions):
            if i == element1:
                li = []
                for j,p2 in zip(self.atoms,self.cart_positions):
                    if j == element2:
                        li.append(sqrt((p2[0]-p[0])**2+(p2[1]-p[1])**2+(p2[2]-p[2])**2))
                ma = max(li)
                p += self.lattice[0]+self.lattice[1]+self.lattice[2]
                for n,(k,p3) in enumerate(zip(s.atoms,s.cart_positions)):
                    if (n<13*self.number_of_atoms or n>=14*self.number_of_atoms) and k == element2:
                        dis = sqrt((p3[0]-p[0])**2+(p3[1]-p[1])**2+(p3[2]-p[2])**2)
                        if dis < max(li) + 0.01:
                            li.append(dis)
                li.sort()
                bond.append((f,li))
        if filename is not None:
            with open(filename,'w') as f:
                for i in bond:
                    f.write('{:s}\n'.format(str(i)))

    

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
        n = 8*self.number_of_atoms
        for i in (0,1):
            for j in (0,1):
                for k in (0,1):
                    if i+j+k>0:
                        inc = i * self.lattice[0] + j * self.lattice[1] + k * self.lattice[2]
                        p = pc + inc
                        pos = np.vstack((pos,p))
        ppp = (pos.__array_interface__['data'][0] + np.arange(pos.shape[0]) * pos.strides[0]).astype(np.uintp)
        _doublepp = ndpointer(dtype=np.uintp, ndim=1, flags='C')
        lib.min_bond_length.argtypes = [_doublepp, ctypes.c_int]
        lib.min_bond_length.restype = ctypes.c_double
        return lib.min_bond_length(ppp,n)
    
    @property
    def volume(self):
        return fabs(np.dot(self.lattice[0], np.cross(self.lattice[1],self.lattice[2])))
    
    @property
    def atom_density(self):
        return self.number_of_atoms/self.volume

    @property
    def adjacency_matrix(s):
        '''calculating galvez matrix for a given structure'''
        def equi(i,amp,cal,na):
            no_cell = i//na
            x = no_cell//(amp*amp)
            y = (no_cell%(amp*amp))//amp
            z = no_cell%amp
            mid_l = int((amp-cal)/2)
            mid_h = int((amp+cal)/2)-1
            while x < mid_l:
                x += cal
            while x > mid_h:
                x -= cal
            while y < mid_l:
                y += cal
            while y > mid_h:
                y -= cal
            while z < mid_l:
                z += cal
            while z > mid_h:
                z -= cal
            return i%na + na * (z-mid_l + cal * (y-mid_l) + cal * cal * (x-mid_l))

        na = len(s.atoms)
        amp = 1
        li = [0]*na
        if na < 3:
            amp += 2
        while True:
            #print('place 1:{:d}'.format(amp))
            index = amp ** 3 // 2
            beg = index * na
            reg = range(beg, beg + na)
            ss = s.multiple((amp, amp, amp))
            try:
                vor = Voronoi(ss.cart_positions)
                break
            except:
                amp += 2
                if amp > 10:
                    raise ValueError
                continue
        amp -= 2
        while True:
            amp += 2
            #print('place 2:{:d}'.format(amp))
            if amp==1:
                amp += 2
            if amp > 10:
                raise ValueError
            s_new = s.multiple((amp,amp,amp))
            atoms = s_new.atoms
            li_new = [0] * na
            vor = Voronoi(s_new.cart_positions)
            middle = set(range(na*int(amp*amp*amp//2), na*int(amp*amp*amp//2+1)))
            for i in vor.ridge_points:
                if i[0] in middle:
                    at_tup = tuple(sorted([atoms[i[0]], atoms[i[1]]]))
                    dis = np.linalg.norm(vor.points[i[0]] - vor.points[i[1]])
                    #print(i[0]%na, i[1]%na, atoms[i[0]], atoms[i[1]], dis)
                    if at_tup in vesta_dic and vesta_dic[at_tup] >= dis:
                        li_new[i[0]%na] += 1
                if i[1] in middle:
                    at_tup = tuple(sorted([atoms[i[0]], atoms[i[1]]]))
                    dis = np.linalg.norm(vor.points[i[0]] - vor.points[i[1]])
                    if at_tup in vesta_dic and vesta_dic[at_tup] >= dis:
                        li_new[i[1] % na] += 1
            #print(li)
            #print(li_new)
            if li_new == li:
                break
            li = li_new
        
        amp -= 2
        cal = -1
        dif = int((amp-cal)/2)
        while True:
            amp += 2
            #print('place 3:{:d}'.format(amp))
            if amp > 10:
                raise ValueError
            cal += 2
            s_find = s.multiple((amp,amp,amp))
            vor = Voronoi(s_find.cart_positions)
            mid = list(i for i in range(amp*amp*amp*na) if dif<=(i//na)%amp<(amp+cal)/2 and dif*amp<=(i//na)%(amp*amp)<(amp+cal)*amp/2 and dif<=(i//na)//(amp*amp)<(amp+cal)/2)
            na_new = na*cal*cal*cal
            ad = np.zeros((na_new,na_new))
            rmat = np.zeros((na_new, na_new))
            atoms = s_find.atoms
            conn = [[] for i in range(na_new)]
            for i in vor.ridge_points:
                if i[0] in mid:
                    at_tup = tuple(sorted([atoms[i[0]], atoms[i[1]]]))
                    dis = np.linalg.norm(vor.points[i[0]] - vor.points[i[1]])
                    #print(p0,p1,dis)
                    if at_tup in vesta_dic and vesta_dic[at_tup] >= dis:
                        p0 = equi(i[0],amp,cal,na)
                        p1 = equi(i[1],amp,cal,na)
                        conn[p0].append(p1)
                        ad[p0][p1]= ad[p1][p0]= 1
                        rmat[p0][p1] = rmat[p1][p0] = dis
                if i[1] in mid:
                    at_tup = tuple(sorted([atoms[i[0]], atoms[i[1]]]))
                    dis = np.linalg.norm(vor.points[i[0]] - vor.points[i[1]])
                    #print(p0, p1, dis)
                    if at_tup in vesta_dic and vesta_dic[at_tup] >= dis:
                        p0 = equi(i[0], amp, cal, na)
                        p1 = equi(i[1], amp, cal, na)
                        conn[p1].append(p0)
                        ad[p0][p1] = ad[p1][p0] = 1
                        rmat[p0][p1] = rmat[p1][p0] = dis
            #print(conn)
            for i in conn:
                if len(set(i))<len(i):
                    break
            else:
                break
        
        return ad, rmat

    @property
    def a(self):
        return self.__lengths[0]

    @property
    def b(self):
        return self.__lengths[1]

    @property
    def c(self):
        return self.__lengths[2]

    @property
    def alpha(self):
        return self.__angles[0]

    @property
    def beta(self):
        return self.__angles[1]

    @property
    def gamma(self):
        return self.__angles[2]

    @property
    def simplified_formula(self):
        eles = list(self.element_species)
        nums = list(self.num_per_species)
        if len(nums) > 1:
            tmp = nums[0]
            for i in range(1, len(nums)):
                tmp = gcd(tmp, nums[i])
            for i in range(len(nums)):
                nums[i] = int(nums[i] / tmp)
            formula_li = sorted(list(zip(eles, nums)), key=lambda x: x[0])
            for i in range(len(formula_li)):
                formula_li[i] = formula_li[i][0] + '{:d}'.format(formula_li[i][1])
            return ''.join(formula_li)
        else:
            return '{:s}1'.format(eles[0])

    @property
    def alloic(self):
        metal = False
        non_metal = False
        for i in self.element_species:
            if Element(i).metal is True:
                metal = True
            else:
                non_metal = True
        if metal and non_metal:
            return 1
        if metal:
            return 0
        return 2

    @property
    def exist_covalent_bond(self):
        count = 0
        for i in self.atoms:
            ele = Element(i)
            if ele.metal is True:
                count += ele.max_val
            else:
                count -= (32-ele.atomic_coordination[1])
        return count < 0

    @property
    def has_halo(self):
        for i in self.atoms:
            ele = Element(i)
            if ele.atomic_coordination[1] == 31:
                return True
        return False

    @property
    def metal_elements(self):
        metal = 0
        for i in self.element_species:
            if Element(i).metal is True:
                metal += 1
        return metal