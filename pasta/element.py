# coding: utf-8
# Copyright © 2016 YunXing Zuo, WeiJi Hsiao

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
        if arg in elements_dict:
            self.__atomic_symbol = arg
            self.__atomic_number = elements_dict[arg][0]
        else:
            self.__atomic_number = arg
            self.__atomic_symbol = elements_list[arg]

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
    def atomic_coordination(self):
        return elements_dict[self.atomic_symbol][1]
    
    @property
    def atomic_mass(self):
        return elements_dict[self.atomic_symbol][2]


elements_list = [
    'None',
    'H',
	'He',
	'Li',
	'Be',
	'B',
	'C',
	'N',
	'O',
	'F',
	'Ne',
	'Na',
	'Mg',
	'Al',
	'Si',
	'P',
	'S',
	'Cl',
	'Ar',
	'K',
	'Ca',
	'Sc',
	'Ti',
	'V',
	'Cr',
	'Mn',
	'Fe',
	'Co',
	'Ni',
	'Cu',
	'Zn',
	'Ga',
	'Ge',
	'As',
	'Se',
	'Br',
	'Kr',
	'Rb',
	'Sr',
	'Y',
	'Zr',
	'Nb',
	'Mo',
	'Tc',
	'Ru',
	'Rh',
	'Pd',
	'Ag',
	'Cd',
	'In',
	'Sn',
	'Sb',
	'Te',
	'I',
	'Xe',
	'Cs',
	'Ba',
	'La',
	'Ce',
	'Pr',
	'Nd',
	'Pm',
	'Sm',
	'Eu',
	'Gd',
	'Tb',
	'Dy',
	'Ho',
	'Er',
	'Tm',
	'Yb',
	'Lu',
	'Hf',
	'Ta',
	'W',
	'Re',
	'Os',
	'Ir',
	'Pt',
	'Au',
	'Hg',
	'Tl',
	'Pb',
	'Bi',
	'Po',
	'At',
	'Rn',
	'Fr',
	'Ra',
	'Ac',
	'Th',
	'Pa',
	'U',
	'Np',
	'Pu',
	'Am',
	'Cm',
	'Bk',
	'Cf',
	'Es',
	'Fm',
	'Md',
	'No',
	'Lr',
	'Rf',
	'Db',
	'Sg',
	'Bh',
	'Hs',
	'Mt',
	'Ds',
	'Rg',
	'Cn',
	'Uut',
	'Uuq',
	'Uup',
	'Uuh',
	'Uus',
	'Uuo',
]

elements_dict = {
    'H':(1, (1, 1), 1.00794),
	'He':(2, (1, 32), 4.002602),
	'Li':(3, (2, 1), 6.941),
	'Be':(4, (2, 2), 9.012182),
	'B':(5, (2, 27), 10.811),
	'C':(6, (2, 28), 12.0107),
	'N':(7, (2, 29), 14.0067),
	'O':(8, (2, 30), 15.9994),
	'F':(9, (2, 31), 18.9984032),
	'Ne':(10, (2, 32), 20.1797),
	'Na':(11, (3, 1), 22.98976928),
	'Mg':(12, (3, 2), 24.305),
	'Al':(13, (3, 27), 26.9815386),
	'Si':(14, (3, 28), 28.0855),
	'P':(15, (3, 29), 30.973762),
	'S':(16, (3, 30), 32.065),
	'Cl':(17, (3, 31), 35.453),
	'Ar':(18, (3, 32), 39.948),
	'K':(19, (4, 1), 39.0983),
	'Ca':(20, (4, 2), 40.078),
	'Sc':(21, (4, 17), 44.955912),
	'Ti':(22, (4, 18), 47.867),
	'V':(23, (4, 19), 50.9415),
	'Cr':(24, (4, 20), 51.9961),
	'Mn':(25, (4, 21), 54.938045),
	'Fe':(26, (4, 22), 55.845),
	'Co':(27, (4, 23), 58.933195),
	'Ni':(28, (4, 24), 58.6934),
	'Cu':(29, (4, 25), 63.546),
	'Zn':(30, (4, 26), 65.38),
	'Ga':(31, (4, 27), 69.723),
	'Ge':(32, (4, 28), 72.64),
	'As':(33, (4, 29), 74.9216),
	'Se':(34, (4, 30), 78.96),
	'Br':(35, (4, 31), 79.904),
	'Kr':(36, (4, 32), 83.798),
	'Rb':(37, (5, 1), 85.4678),
	'Sr':(38, (5, 2), 87.62),
	'Y':(39, (5, 17), 88.90585),
	'Zr':(40, (5, 18), 91.224),
	'Nb':(41, (5, 19), 92.90638),
	'Mo':(42, (5, 20), 95.96),
	'Tc':(43, (5, 21), None),
	'Ru':(44, (5, 22), 101.07),
	'Rh':(45, (5, 23), 102.9055),
	'Pd':(46, (5, 24), 106.42),
	'Ag':(47, (5, 25), 107.8682),
	'Cd':(48, (5, 26), 112.411),
	'In':(49, (5, 27), 114.818),
	'Sn':(50, (5, 28), 118.71),
	'Sb':(51, (5, 29), 121.76),
	'Te':(52, (5, 30), 127.6),
	'I':(53, (5, 31), 126.90447),
	'Xe':(54, (5, 32), 131.293),
	'Cs':(55, (6, 1), 132.9054519),
	'Ba':(56, (6, 2), 137.327),
	'La':(57, (6, 3), 138.90547),
	'Ce':(58, (6, 4), 140.116),
	'Pr':(59, (6, 5), 140.90765),
	'Nd':(60, (6, 6), 144.242),
	'Pm':(61, (6, 7), None),
	'Sm':(62, (6, 8), 150.36),
	'Eu':(63, (6, 9), 151.964),
	'Gd':(64, (6, 10), 157.25),
	'Tb':(65, (6, 11), 158.92535),
	'Dy':(66, (6, 12), 162.5),
	'Ho':(67, (6, 13), 164.93032),
	'Er':(68, (6, 14), 167.259),
	'Tm':(69, (6, 15), 168.93421),
	'Yb':(70, (6, 16), 173.054),
	'Lu':(71, (6, 17), 174.9668),
	'Hf':(72, (6, 18), 178.49),
	'Ta':(73, (6, 19), 180.94788),
	'W':(74, (6, 20), 183.84),
	'Re':(75, (6, 21), 186.207),
	'Os':(76, (6, 22), 190.23),
	'Ir':(77, (6, 23), 192.217),
	'Pt':(78, (6, 24), 195.084),
	'Au':(79, (6, 25), 196.966569),
	'Hg':(80, (6, 26), 200.59),
	'Tl':(81, (6, 27), 204.3833),
	'Pb':(82, (6, 28), 207.2),
	'Bi':(83, (6, 29), 208.9804),
	'Po':(84, (6, 30), None),
	'At':(85, (6, 31), None),
	'Rn':(86, (6, 32), None),
	'Fr':(87, (7, 1), None),
	'Ra':(88, (7, 2), None),
	'Ac':(89, (7, 3), None),
	'Th':(90, (7, 4), 232.03806),
	'Pa':(91, (7, 5), 231.03588),
	'U':(92, (7, 6), 238.02891),
	'Np':(93, (7, 7), None),
	'Pu':(94, (7, 8), None),
	'Am':(95, (7, 9), None),
	'Cm':(96, (7, 10), None),
	'Bk':(97, (7, 11), None),
	'Cf':(98, (7, 12), None),
	'Es':(99, (7, 13), None),
	'Fm':(100, (7, 14), None),
	'Md':(101, (7, 15), None),
	'No':(102, (7, 16), None),
	'Lr':(103, (7, 17), None),
	'Rf':(104, (7, 18), None),
	'Db':(105, (7, 19), None),
	'Sg':(106, (7, 20), None),
	'Bh':(107, (7, 21), None),
	'Hs':(108, (7, 22), None),
	'Mt':(109, (7, 23), None),
	'Ds':(110, (7, 24), None),
	'Rg':(111, (7, 25), None),
	'Cn':(112, (7, 26), None),
	'Uut':(113, (7, 27), None),
	'Uuq':(114, (7, 28), None),
	'Uup':(115, (7, 29), None),
	'Uuh':(116, (7, 30), None),
	'Uus':(117, (7, 31), None),
	'Uuo':(118, (7, 32), None),
}