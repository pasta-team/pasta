# coding: utf-8
# Copyright © 2016 YunXing Zuo, WeiJi Hsiao

from .element_property import *

__author__ = 'YunXing Zuo, WeiJi Hsiao'
__email__ = 'weiji.hsiao@gmail.com'
__date__ = 'Oct. 25, 2016'


class Element(object):
    def __init__(self, arg):
        """
        Create a Element object.
        Args:
            arg (str/int): An atomic symbol string or an atomic number.
        """
        if isinstance(arg, int):
            self.__atomic_symbol = NUMBER[arg]
            self.__atomic_number = arg
        elif isinstance(arg, str):
            self.__atomic_symbol = arg
            self.__atomic_number = SYMBOL_TO_NUM[arg]

    def __repr__(self):
        return self.__atomic_symbol

    def __eq__(self, other):
        if isinstance(other, Element):
            if self.atomic_number == other.atomic_number:
                return True
        return False

    def __hash__(self):
        return self.atomic_number

    def __gt__(self, other):
        if isinstance(other, Element):
            if self.atomic_number > other.atomic_number:
                return True
        return False

    def __lt__(self, other):
        if isinstance(other, Element):
            if self.atomic_number < other.atomic_number:
                return True
        return False

    @property
    def atomic_symbol(self):
        return self.__atomic_symbol

    @property
    def atomic_number(self):
        return self.__atomic_number

    @property
    def metal(self):
        return self.atomic_symbol in METAL

    @property
    def transitional_metal(self):
        return self.atomic_symbol in TRANSITION_METAL

    @property
    def max_val(self):
        if self.metal is True:
            return METAL[self.atomic_symbol]

    @property
    def atomic_coordination(self):
        return ATOMIC_COORDINATION[self.atomic_number-1]

    @property
    def electronegativity(self):
        return ELECTRONEGATIVITY[self.atomic_number-1]

    @property
    def valence_electron_descriptor(self):
        return VALENCE_ELECTRON_DESCRIPTOR[self.atomic_number-1]

    @property
    def covalent_radii(self):
        """
        Cordero, B. et al. Covalent radii revisited. Dalton Trans. 2008, 2832–2838 (2008)
        average is taken when there are different values for hs&ls or sp&sp2&sp3
        """
        return COVALENT_RADII[self.atomic_number-1]

    @property
    def heat_capacity(self):
        return HEAT_CAPACITY[self.atomic_number-1]

    @property
    def electron_affinity(self):
        """
        CRC, 92nd edition, 10-147
        simply set to 0 when 'not stable'
        ignoring '>' and '<'
        missing value set to 0
        """
        return ELECTRON_AFFINITY[self.atomic_number-1]

    @property
    def polarizability(self):
        """
        CRC, 92nd edition, 10-186
        average is taken when multivalues available
        """
        return POLARIZABILITY[self.atomic_number-1]

    @property
    def ionization_potential(self):
        """
        CRC, 92nd edition,10-196
        non-exist values are set to 10000
        """
        return IONIZATION_POTENTIAL[self.atomic_number-1]

    @property
    def standard_molar_entropy(self):
        """
        CRC, 92nd edition,5-4
        missing value for Tl?
        """
        return STANDARD_MOLAR_ENTROPY[self.atomic_number-1]

    #@property
    def effective_nuclear_charge(self):
        """
        wikipedia, 2017/9/14
        """
        #return elements_data[self.atomic_symbol][9]
        pass

    @property
    def thermal_conductivity(self):
        """
        wikipedia, 2017/9/14
        WEL
        """
        return THERMAL_CONDUCTIVITY[self.atomic_number-1]

    @property
    def molar_volume(self):
        """
        http://periodictable.com/Properties/A/MolarVolume.html
        """
        return MOLAR_VOLUME[self.atomic_number-1]

    @property
    def atomic_mass(self):
        return ATOMIC_MASS[self.atomic_number-1]

    @property
    def enthalpy_of_fusion(self):
        return ENTHALPY_OF_FUSION[self.atomic_number-1]

    @property
    def enthalpy_of_atomization(self):
        return ENTHALPY_OF_ATOMIZATION[self.atomic_number-1]

    @property
    def enthalpy_of_vaporization(self):
        return ENTHALPY_OF_VAPORIZATION[self.atomic_number-1]

    def get_property(self, s):
        if s == 'row_number':
            return self.atomic_coordination[0]
        if s == 'group_number':
            return self.atomic_coordination[1]
        if s == 'electronegativity':
            return self.electronegativity
        if s == 'valence_electron_descriptor':
            return self.valence_electron_descriptor
        if s == 'covalent_radii':
            return self.covalent_radii
        if s == 'heat_capacity':
            return self.heat_capacity
        if s == 'electron_affinity':
            return self.electron_affinity
        if s == 'polarizability':
            return self.polarizability
        if s == 'IP1':
            return self.ionization_potential[0]
        if s == 'IP2':
            return self.ionization_potential[1]
        if s == 'IP3':
            return self.ionization_potential[2]
        if s == 'standard_molar_entropy':
            return self.standard_molar_entropy
        '''if s == 'effective_nuclear_charge':
            return self.effective_nuclear_charge'''
        if s == 'thermal_conductivity':
            return self.thermal_conductivity
        if s == 'molar_volume':
            return self.molar_volume
        if s == 'atomic_mass':
            return self.atomic_mass
        if s == 'enthalpy_of_fusion':
            return self.enthalpy_of_fusion
        if s == 'enthalpy_of_atomization':
            return self.enthalpy_of_atomization
        if s == 'enthalpy_of_vaporization':
            return self.enthalpy_of_vaporization


'''elements_dict = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Kr": 36,
    "Rb": 37,
    "Sr": 38,
    "Y": 39,
    "Zr": 40,
    "Nb": 41,
    "Mo": 42,
    "Tc": 43,
    "Ru": 44,
    "Rh": 45,
    "Pd": 46,
    "Ag": 47,
    "Cd": 48,
    "In": 49,
    "Sn": 50,
    "Sb": 51,
    "Te": 52,
    "I": 53,
    "Xe": 54,
    "Cs": 55,
    "Ba": 56,
    "La": 57,
    "Ce": 58,
    "Pr": 59,
    "Nd": 60,
    "Pm": 61,
    "Sm": 62,
    "Eu": 63,
    "Gd": 64,
    "Tb": 65,
    "Dy": 66,
    "Ho": 67,
    "Er": 68,
    "Tm": 69,
    "Yb": 70,
    "Lu": 71,
    "Hf": 72,
    "Ta": 73,
    "W": 74,
    "Re": 75,
    "Os": 76,
    "Ir": 77,
    "Pt": 78,
    "Au": 79,
    "Hg": 80,
    "Tl": 81,
    "Pb": 82,
    "Bi": 83,
    "Po": 84,
    "At": 85,
    "Rn": 86,
    "Fr": 87,
    "Ra": 88,
    "Ac": 89,
    "Th": 90,
    "Pa": 91,
    "U": 92,
    "Np": 93,
    "Pu": 94,
    "Am": 95,
    "Cm": 96,
    "Bk": 97,
    "Cf": 98,
    "Es": 99,
    "Fm": 100,
    "Md": 101,
    "No": 102,
    "Lr": 103,
    "Rf": 104,
    "Db": 105,
    "Sg": 106,
    "Bh": 107,
    "Hs": 108,
    "Mt": 109,
    "Ds": 110,
    "Rg": 111,
    "Cn": 112,
    "Uut": 113,
    "Uuq": 114,
    "Uup": 115,
    "Uuh": 116,
    "Uus": 117,
    "Uuo": 118,
}

elements_data = {
    "H": [(1, 1), 2.20, 1, 0.31, 28.836, 0.754195, 0.666793, (13.598443, 10000, 10000), 130.7, 1.000, 0.1805, 0.01121],
    "He": [(1, 32), None, 2, 0.28, 20.786, 0, 0.2050520, (24.587387, 54.417760, 10000), 126.2, 1.688, 0.1513, 0.022424],
    "Li": [(2, 1), 0.98, 3, 1.28, 24.860, 0.618049, 24.33, (5.391719, 75.6400, 122.45429), 29.1, 1.279, 85, 1.297e-5],
    "Be": [(2, 2), 1.57, 2, 0.96, 16.443, 0, 5.60, (9.32270, 18.21114, 153.89661), 9.5, 1.912, 190, 4.8767e-6],
    "B": [(2, 27), 2.04, 3, 0.84, 11.087, 0.279723, 3.03, (8.29802, 25.1548, 37.93064), 5.9, 2.421, 27, 4.3947e-6],
    "C": [(2, 28), 2.55, 4, 0.73, 8.517, 1.262119, 1.67, (11.26030, 24.3833, 47.8878), 5.7, 3.136, 140, 5.3145e-6],
    "N": [(2, 29), 3.04, 5, 0.71, 29.124, 0, 1.10, (14.5341, 29.6013, 47.44924), 191.6, 3.834, 0.02583, 0.011196],
    "O": [(2, 30), 3.44, 6, 0.66, 29.378, 1.4611135, 0.802, (13.61805, 35.1211, 54.9355), 205.2, 4.453, 0.02658, 0.011196],
    "F": [(2, 31), 3.98, 7, 0.57, 31.304, 3.4011895, 0.557, (17.4228, 34.9708, 62.7084), 202.8, 5.100, 0.0277, 0.011202],
    "Ne": [(2, 32), None, 8, 0.58, 20.786, 0, 0.39432, (21.56454, 40.96296, 63.45), 146.3, 5.758, 0.0491, 0.02242],
    "Na": [(3, 1), 0.93, 7, 1.66, 28.230, 0.547926, 24.11, (5.139076, 47.2864, 71.6200), 51.3, 2.507, 140, 2.375e-5],
    "Mg": [(3, 2), 1.31, 2, 1.41, 24.869, 0, 10.8, (7.646235, 15.03527, 80.1437), 32.7, 3.308, 160, 1.3984e-5],
    "Al": [(3, 27), 1.61, 3, 1.21, 24.20, 0.43283, 6.8, (5.985768, 18.82855, 28.44765), 28.3, 4.066, 235, 9.99e-6],
    "Si": [(3, 28), 1.90, 4, 1.11, 19.99, 1.3895213, 5.53, (8.15168, 16.34584, 33.49302), 18.8, 4.285, 150, 1.2054e-5],
    "P": [(3, 29), 2.19, 5, 1.07, 23.824, 0.7465, 3.63, (10.48669, 19.7695, 30.2027), 41.1, 4.886, 0.236, 1.6991e-5],
    "S": [(3, 30), 2.58, 6, 1.05, 22.7, 2.07710418, 2.90, (10.36001, 23.33788, 34.79), 32.1, 5.482, 0.205, 1.636e-5],
    "Cl": [(3, 31), 3.16, 7, 1.02, 33.949, 3.612724, 2.18, (12.96763, 23.8136, 39.61), 223.1, 6.116, 0.0089, 0.011031],
    "Ar": [(3, 32), None, 8, 1.06, 20.786, 0, 1.6411, (15.759610, 27.62966, 40.74), 154.8, 6.764, 0.01772, 0.022392],
    "K": [(4, 1), 0.82, 9, 2.03, 29.6, 0.50147, 43.4, (4.3406633, 31.63, 45.806), 64.7, 3.495, 100, 4.568e-5],
    "Ca": [(4, 2), 1.00, 10, 1.76, 25.929, 0.02455, 25.7, (6.11316, 11.87172, 50.9131), 41.6, 4.398, 200, 2.5857e-5],
    "Sc": [(4, 17), 1.36, 11, 1.70, 25.52, 0.188, 17.8, (6.56149, 12.79977, 24.75666), 34.6, 7.120, 16, 1.5061e-5],
    "Ti": [(4, 18), 1.54, 12, 1.60, 25.06, 0.079, 14.6, (6.82812, 13.5755, 27.4917), 30.7, 8.141, 22, 1.0621e-5],
    "V": [(4, 19), 1.63, 13, 1.53, 24.89, 0.525, 12.4, (6.74619, 14.618, 29.311), 28.9, 8.983, 31, 8.3374e-6],
    "Cr": [(4, 20), 1.66, 12, 1.39, 23.35, 0.666, 11.6, (6.76651, 16.4857, 30.96), 23.8, 9.757, 94, 7.2824e-6],
    "Mn": [(4, 21), 1.55, 13, 1.50, 26.32, 0, 9.4, (7.43402, 15.6400, 33.668), 32.0, 10.528, 7.8, 7.3545e-6],
    "Fe": [(4, 22), 1.86, 8, 1.42, 25.1, 0.151, 8.4, (7.9024, 16.1877, 30.652), 27.3, 11.180, 80, 7.0923e-6],
    "Co": [(4, 23), 1.88, 9, 1.38, 24.81, 0.662, 7.5, (7.88101, 17.084, 33.50), 30.0, 11.855, 100, 6.62e-6],
    "Ni": [(4, 24), 1.91, 10, 1.24, 26.07, 1.156, 6.8, (7.6398, 18.16884, 35.19), 29.9, 12.530, 91, 6.5888e-6],
    "Cu": [(4, 25), 1.90, 11, 1.32, 24.440, 1.235, 6.2, (7.72638, 20.2924, 36.841), 33.2, 13.201, 400, 7.124e-6],
    "Zn": [(4, 26), 1.65, 12, 1.22, 25.39, 0, 5.8, (9.394199, 17.96439, 39.723), 41.6, 13.878, 120, 9.161e-6],
    "Ga": [(4, 27), 1.81, 13, 1.22, 26.03, 0.43, 8.12, (5.999301, 20.51515, 30.7258), 40.8, 6.222, 29, 1.1809e-5],
    "Ge": [(4, 28), 2.01, 14, 1.20, 23.222, 1.232712, 5.84, (7.89943, 15.93461, 34.2241), 31.1, 6.780, 60, 1.3646e-5],
    "As": [(4, 29), 2.18, 5, 1.19, 24.64, 0.804, 4.31, (9.7886, 18.5892, 28.351), 35.1, 7.449, 50, 1.3082e-5],
    "Se": [(4, 30), 2.55, 6, 1.20, 25.363, 2.020670, 3.77, (9.75239, 21.19, 30.8204), 42.4, 8.287, 0.52, 1.6385e-5],
    "Br": [(4, 31), 2.96, 7, 1.20, 75.69, 3.362588, 3.05, (11.8138, 21.591, 36), 152.2, 9.028, 0.12, 2.561e-5],
    "Kr": [(4, 32), 3.00, 8, 1.16, 20.786, 0, 2.4844, (13.99961, 24.35984, 36.950), 164.1, 9.338, 0.00943, 0.02235],
    "Rb": [(5, 1), 0.82, 9, 2.20, 31.06, 0.48592, 47.3, (4.177128, 27.2895, 40), 76.8, 4.985, 58, 5.5788e-5],
    "Sr": [(5, 2), 0.95, 10, 1.95, 26.79, 0.048, 25.5, (5.69485, 11.0301, 42.89),55.0, 6.071, 35, 3.3316e-5],
    "Y": [(5, 17), 1.22, 11, 1.90, 26.53, 0.307, 22.7, (6.2173, 12.224, 20.52), 44.4, 15.958, 17, 1.9881e-5],
    "Zr": [(5, 18), 1.33, 12, 1.75, 25.36, 0.426, 17.9, (6.63390, 13.1, 22.99), 39.0, 13.072, 23, 1.4011e-5],
    "Nb": [(5, 19), 1.60, 13, 1.64, 24.6, 0.916, 15.7, (6.75885, 14.0, 25.04), 36.4, 11.238, 54, 1.0841e-5],
    "Mo": [(5, 20), 2.16, 14, 1.54, 24.06, 0.748, 12.8, (7.09243, 16.16, 27.13), 28.7, 11.392, 139, 9.333e-6],
    "Tc": [(5, 21), 1.90, 13, 1.47, None, 0.55, 11.4, (7.28, 15.26, 29.54), None, 12.882, 51, 8.522e-6],
    "Ru": [(5, 22), 2.20, 14, 1.46, 24.06, 1.05, 9.6, (7.36050, 16.76, 28.47), 28.5, 12.813, 120, 8.1706e-6],
    "Rh": [(5, 23), 2.28, 15, 1.42, 24.98, 1.137, 8.6, (7.45890, 18.08, 31.06), 31.5, 13.442, 150, 8.2655e-6],
    "Pd": [(5, 24), 2.20, 10, 1.39, 25.98, 0.562, 4.8, (8.3369, 19.43, 32.93), 37.6, 13.618, 72, 8.8514e-6],
    "Ag": [(5, 25), 1.93, 11, 1.45, 25.350, 1.302, 7.2, (7.57623, 21.47746, 34.83), 42.6, 14.763, 430, 1.0283e-5],
    "Cd": [(5, 26), 1.69, 12, 1.44, 26.020, 0, 7.3, (8.99382, 16.90831, 37.48), 51.8, 15.877, 97, 1.2995e-5],
    "In": [(5, 27), 1.78, 13, 1.42, 26.74, 0.3, 9.6, (5.78636, 18.8703, 28.03), 57.8, 8.470, 82, 1.5707e-5],
    "Sn": [(5, 28), 1.96, 14, 1.39, 26.99, 1.112067, 7.06, (7.34392, 14.6322, 30.50260), 51.2, 9.102, 67, 1.6239e-5],
    "Sb": [(5, 29), 2.05, 5, 1.39, 25.23, 1.046, 6.6, (8.60839, 16.63, 25.3), 45.7, 9.995, 24, 1.8181e-5],
    "Te": [(5, 30), 2.10, 6, 1.38, 25.73, 1.970876, 5.5, (9.0096, 18.6, 27.96), 49.7, 10.809, 3, 2.0449e-5],
    "I": [(5, 31), 2.66, 7, 1.39, 54.43, 3.0590463, 5.0, (10.45126, 19.1313, 33), 116.1, 11.612, 0.449, 2.5689e-5],
    "Xe": [(5, 32), 2.60, 8, 1.40, 20.786, 0, 4.044, (12.12984, 20.9750, 32.1230), 169.7, 12.425, 0.00565, 0.0223],
    "Cs": [(6, 1), 0.79, 9, 2.44, 32.210, 0.471626, 59.42, (3.893905, 23.15744, None), 85.2, None, 36, 7.0732e-5],
    "Ba": [(6, 2), 0.89, 10, 2.15, 28.07, 0.14462, 39.7, (5.211664, 10.00383, None), 62.5, None, 18, 3.9125e-5],
    "La": [(6, 3), 1.10, 11, 2.07, 27.11, 0.47, 31.1, (5.5769, 11.059, 19.1773), 56.9, None, 13, 2.2601e-5],
    "Ce": [(6, 4), 1.12, 12, 2.04, 26.94, 0.65, 29.6, (5.5387, 10.85, 20.198), 72.0, None, 11, 2.0947e-5],
    "Pr": [(6, 5), 1.13, 11, 2.03, 27.20, 0.962, 28.2, (5.473, 10.55, 21.624), 73.2, None, 13, 2.1221e-5],
    "Nd": [(6, 6), 1.14, 11, 2.01, 27.45, 1.916, 31.4, (5.5250, 10.72, 22.1), 71.5, None, 17, 2.0576e-5],
    "Pm": [(6, 7), 1.13, 11, 1.99, None, None, 30.1, (5.582, 10.90, 22.3), None, None, 15, 1.9961e-5],
    "Sm": [(6, 8), 1.17, 11, 1.98, 29.54, None, 28.8, (5.6437, 11.07, 23.4), 69.6, None, 13, 2.0449e-5],
    "Eu": [(6, 9), 1.20, 8, 1.98, 27.66, 0.864, 27.7, (5.67038, 11.25, 24.92), 77.8, None, 14, 2.8979e-5],
    "Gd": [(6, 10), 1.20, 9, 1.96, 37.03, None, 23.5, (6.14980, 12.09, 20.63), 68.1, None, 11, 1.9903e-5],
    "Tb": [(6, 11), 1.10, 9, 1.94, 28.91, 1.165, 25.5, (5.8638, 11.52, 21.91), 73.2, None, 11, 1.9336e-5],
    "Dy": [(6, 12), 1.22, 9, 1.92, 28.16, 0, 24.5, (5.9389, 11.67, 22.8), 75.6, None, 11, 1.9004e-5],
    "Ho": [(6, 13), 1.23, 9, 1.92, 27.15, None, 23.6, (6.0215, 11.80, 22.84), 75.3, None, 16, 1.8753e-5],
    "Er": [(6, 14), 1.24, 9, 1.89, 28.12, None, 22.7, (6.1077, 11.93, 22.74), 73.2, None, 15, 1.8449e-5],
    "Tm": [(6, 15), 1.25, 9, 1.90, 27.03, 1.029, 21.8, (6.18431, 12.05, 23.68), 74.0, None, 17, 1.8124e-5],
    "Yb": [(6, 16), 1.10, 8, 1.87, 26.74, -0.020, 20.9, (6.25416, 12.176, 25.05), 59.9, None, 39, 2.6338e-5],
    "Lu": [(6, 17), 1.27, 9, 1.87, 26.86, 0.34, 21.9, (5.42586, 13.9, 20.9594), 51.0, None, 16, 1.7779e-5],
    "Hf": [(6, 18), 1.30, 10, 1.75, 25.73, 0.017, 16.2, (6.82507, 15, 23.3), 43.6, None, 23, 1.34102e-5],
    "Ta": [(6, 19), 1.50, 11, 1.70, 25.36, 0.322, 13.1, (7.54957, None, None), 41.5, None, 57, 1.08677e-5],
    "W": [(6, 20), 2.36, 12, 1.62, 24.27, 0.81626, 11.1, (7.86403, 16.1, None), 32.6, None, 170, 9.5501e-6],
    "Re": [(6, 21), 1.90, 7, 1.51, 25.48, 0.15, 9.7, (7.83352, None, None), 36.9, None, 48, 8.85856e-6],
    "Os": [(6, 22), 2.20, 8, 1.44, 24.7, 1.1, 8.5, (8.43823, None, None), 32.6, None, 88, 8.4135e-6],
    "Ir": [(6, 23), 2.20, 9, 1.41, 25.1, 1.5638, 7.6, (8.96702, None, None), 35.5, None, 150, 8.4864e-6],
    "Pt": [(6, 24), 2.28, 10, 1.36, 25.86, 2.128, 6.5, (8.9588, 18.563, None), 41.6, None, 72, 9.2498e-6],
    "Au": [(6, 25), 2.54, 11, 1.36, 25.418, 2.30863, 5.8, (9.22553, 20.20, None), 47.4, None, 320, 1.021e-5],
    "Hg": [(6, 26), 2.00, 12, 1.32, 27.983, 0, 5.4, (10.4375, 18.7568, 34.2), 75.9, None, 8.3, 1.48212e-5],
    "Tl": [(6, 27), 1.62, 13, 1.45, 26.32, 0.377, 7.6, (6.108194, 20.4283, 29.83), None, None, 46, 1.72475e-5],
    "Pb": [(6, 28), 1.87, 14, 1.46, 26.84, 0.364, 7.00, (7.41663, 15.03248, 31.9373), 64.8, None, 35, 1.8272e-5],
    "Bi": [(6, 29), 2.02, 15, 1.48, 25.52, 0.942362, 7.4, (7.2855, 16.703, 25.56), 56.7, None, 8, 2.1368e-5],
    "Po": [(6, 30), 2.00, 16, 1.40, None, 1.9, 6.8, (8.414, None, None), None, None, 0.2, 2.2727e-5],
    "At": [(6, 31), 2.20, 17, 1.50, None, 2.8, 6.0, (None, None, None), None, None, 2, None],
    "Rn": [(6, 32), 2.20, 8, 1.50, 20.786, 0, 5.3, (10.7485, None, None), 176.2, None, 0.00361, 0.02282],
    "Fr": [(7, 1), 0.70, 9, 2.60, None, 0.486, 47.8, (4.072741, None, None), 95.4, None, None, None],
    "Ra": [(7, 2), 0.90, 10, 2.21, None, 0.10, 38.3, (5.278423, 10.14715, None), 71.0, None, 19, 4.52e-5],
    "Ac": [(7, 3), 1.10, 11, 2.15, 27.2, 0.35, 32.1, (5.17, 11.75, None), 56.5, None, 12, 2.25422e-5],
    "Th": [(7, 4), 1.30, 12, 2.06, 27.32, None, 32.1, (6.3067, 11.9, 20.0), 51.8, None, 54, 1.97917e-5],
    "Pa": [(7, 5), 1.50, 13, 2.00, None, None, 25.4, (5.89, None, None), 51.9, None, 47, 1.50316e-5],
    "U": [(7, 6), 1.38, 14, 1.96, 27.665, None, 24.9, (6.1941, 10.6, None), 50.2, None, 27, 1.2495e-5],
    "Np": [(7, 7), 1.36, 15, 1.90, None, None, 24.8, (6.2657, None, None), None, None, 6, 1.15892e-5],
    "Pu": [(7, 8), 1.28, 16, 1.87, None, None, 24.5, (6.0260, 11.2, None), None, None, 6, 1.23133e-5],
    "Am": [(7, 9), 1.13, 17, 1.80, None, None, 23.3, (5.9738, None, None), None, None, 10, None],
    "Cm": [(7, 10), 1.28, 18, 1.69, None, None, 23.0, (5.9914, None, None), None, None, None, 1.82828e-5],
    "Bk": [(7, 11), 1.30, None, None, None, None, 22.7, (6.1979, None, None), None, None, 10, 1.67118e-5],
    "Cf": [(7, 12), 1.30, None, None, None, None, 20.5, (6.2817, 11.8, None), None, None, None, 1.662e-5],
    "Es": [(7, 13), 1.30, None, None, None, None, 19.7, (6.42, 12.0, None), None, None, None, None],
    "Fm": [(7, 14), 1.30, None, None, None, None, 23.8, (6.50, None, None), None, None, None, None],
    "Md": [(7, 15), 1.30, None, None, None, None, 18.2, (6.58, None, None), None, None, None, None],
    "No": [(7, 16), 1.30, None, None, None, None, 16.4, (6.65, None, None), None, None, None, None],
    "Lr": [(7, 17), 1.30, None, None, None, None, None, (4.9, None, None), None, None, None, None],
    "Rf": [(7, 18), None, None, None, None, None, None, (6.0, None, None), None, None, None, None],
    "Db": [(7, 19), None, None, None, None, None, None, (None, None, None), None, None, None, None],
    "Sg": [(7, 20), None, None, None, None, None, None, (None, None, None), None, None, None, None],
    "Bh": [(7, 21), None, None, None, None, None, None, (None, None, None), None, None, None, None],
    "Hs": [(7, 22), None, None, None, None, None, None, (None, None, None), None, None, None, None],
    "Mt": [(7, 23), None, None, None, None, None, None, (None, None, None), None, None, None, None],
    "Ds": [(7, 24), None, None, None, None, None, None, (None, None, None), None, None, None, None],
    "Rg": [(7, 25), None, None, None, None, None, None, (None, None, None), None, None, None, None],
    "Cn": [(7, 26), None, None, None, None, None, 4.06, (None, None, None), None, None, None, None],
    "Uut": [(7, 27), None, None, None, None, None, None, (None, None, None), None, None, None, None],
    "Uuq": [(7, 28), None, None, None, None, None, 4.48, (None, None, None), None, None, None, None],
    "Uup": [(7, 29), None, None, None, None, None, None, (None, None, None), None, None, None, None],
    "Uuh": [(7, 30), None, None, None, None, None, None, (None, None, None), None, None, None, None],
    "Uus": [(7, 31), None, None, None, None, None, None, (None, None, None), None, None, None, None],
    "Uuo": [(7, 32), None, None, None, None, 0.056, None, (None, None, None), None, None, None, None],
}

metal = {
    "H":1,
    "Li":1,
    "Be":2,
    "Na":1,
    "Mg":2,
    "Al":3,
    "K":1,
    "Ca":2,
    "Sc":3,
    "Ti":4,
    "V":5,
    "Cr":6,
    "Mn":7,
    "Fe":7,
    "Co":5,
    "Ni":4,
    "Cu":4,
    "Zn":2,
    "Ga":3,
    "Ge":4,
    "Rb":1,
    "Sr":2,
    "Y":3,
    "Zr":4,
    "Nb":5,
    "Mo":6,
    "Tc":7,
    "Ru":8,
    "Rh":6,
    "Pd":6,
    "Ag":4,
    "Cd":2,
    "In":3,
    "Sn":4,
    "Sb":5,
    "Cs":1,
    "Ba":2,
    "La":3,
    "Ce":4,
    "Pr":5,
    "Nd":4,
    "Pm":3,
    "Sm":3,
    "Eu":3,
    "Gd":3,
    "Tb":4,
    "Dy":4,
    "Ho":3,
    "Er":3,
    "Tm":3,
    "Yb":3,
    "Lu":3,
    "Hf":4,
    "Ta":5,
    "W":6,
    "Re":7,
    "Os":8,
    "Ir":5,
    "Pt":6,
    "Au":5,
    "Hg":2,
    "Tl":3,
    "Pb":4,
    "Bi":5,
    "Po":6,
    "Fr":1,
    "Ra":2,
    "Ac":3,
    "Th":4,
    "Pa":5,
    "U":6,
    "Np":7,
    "Pu":7,
    "Am":7,
    "Cm":6,
    "Bk":4,
    "Cf":4,
    "Es":4,
    "Fm":3,
    "Md":3,
    "No":3,
    "Lr":3,
    "Rf":4,
    "Db":5,
    "Sg":6,
    "Bh":7,
    "Hs":8,
    "Mt":3,
    "Ds":3,
    "Rg":3,
    "Cn":3,
    "Uut":3,
    "Uuq":3,
    "Uup":3,
    "Uuh":3,
    "Uus":3,
    "Uuo":3,}

transitional_metal = {'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
                      'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
                      'La','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg'}'''
