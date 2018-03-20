#!/usr/bin/python
'''
# Construct the crystal structure js file for the 3Dmodel
# include the data of atom,bond,polyhedra
# auto scale the size
'''

from numpy import random
import os
import math
# message of element
elementtable = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn']

# vesta value of radius
element_radious = [0.46, 1.22, 1.57, 1.12, 0.81, 0.77, 0.74, 0.74, 0.72, 1.6, 1.91, 1.6, 1.43, 1.18, 1.1, 1.04, 0.99, 1.92, 2.35, 1.97, 1.64, 1.47, 1.35, 1.29, 1.37, 1.26, 1.25, 1.25, 1.28, 1.37, 1.53, 1.22, 1.21, 1.04, 1.14, 1.98, 2.5, 2.15, 1.82, 1.6, 1.47, 1.4, 1.35, 1.34, 1.34, 1.37, 1.44, 1.52, 1.67, 1.58, 1.41, 1.37, 1.33, 2.18, 2.72, 2.24, 1.88, 1.82, 1.82, 1.82, 1.81, 1.81, 2.06, 1.79, 1.77, 1.77, 1.76, 1.75, 1.0, 1.94, 1.72, 1.59, 1.47, 1.41, 1.37, 1.35, 1.36, 1.39, 1.44, 1.55, 1.71, 1.75, 1.82, 1.77, 0.62, 0.8, 1.0, 2.35, 2.03, 1.8, 1.63, 1.56, 1.56, 1.64, 1.73, 0.8]
element_valence = [1,1,1,2,3,4,5,2,1,1,1,2,3,4,5,6,5,1,1,2,3,4,5,6,7,3,3,2,2,2,3,4,5,6,5,1,1,2,3,4,5,6,7,4,4,4,1,2,3,4,5,6,6,1,1,2,3,4,3,3,3,3,2,3,4,3,3,3,3,3,3,4,5,6,7,6,4,4,3,2,3,4,5,6,5,1,1,2,3,4,5,6,5,5,3,3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,1,1,1]

# vesta color setting
element_color = [[1.0, 0.8, 0.8],[0.98907, 0.91312, 0.81091],[0.52731, 0.87953, 0.4567],[0.37147, 0.8459, 0.48292],[0.1249, 0.63612, 0.05948],[0.5043, 0.28659, 0.16236],[0.69139, 0.72934, 0.9028],[0.99997, 0.01328, 0.0],[0.69139, 0.72934, 0.9028],[0.99954, 0.21788, 0.71035],[0.97955, 0.86618, 0.23787],[0.98773, 0.48452, 0.0847],[0.50718, 0.70056, 0.84062],[0.10596, 0.23226, 0.98096],[0.75557, 0.61256, 0.76425],[1.0, 0.98071, 0.0],[0.19583, 0.98828, 0.01167],[0.81349, 0.99731, 0.77075],[0.63255, 0.13281, 0.96858],[0.35642, 0.58863, 0.74498],[0.71209, 0.3893, 0.67279],[0.47237, 0.79393, 1.0],[0.9, 0.1, 0.0],[0.0, 0.0, 0.62],[0.66148, 0.03412, 0.62036],[0.71051, 0.44662, 0.00136],[0.0, 0.0, 0.68666],[0.72032, 0.73631, 0.74339],[0.1339, 0.28022, 0.86606],[0.56123, 0.56445, 0.50799],[0.62292, 0.89293, 0.45486],[0.49557, 0.43499, 0.65193],[0.45814, 0.81694, 0.34249],[0.6042, 0.93874, 0.06122],[0.49645, 0.19333, 0.01076],[0.98102, 0.75805, 0.95413],[1.0, 0.0, 0.6],[0.0, 1.0, 0.15259],[0.40259, 0.59739, 0.55813],[0.0, 1.0, 0.0],[0.29992, 0.70007, 0.46459],[0.70584, 0.52602, 0.68925],[0.80574, 0.68699, 0.79478],[0.81184, 0.72113, 0.68089],[0.80748, 0.82205, 0.67068],[0.75978, 0.76818, 0.72454],[0.72032, 0.73631, 0.74339],[0.95145, 0.12102, 0.86354],[0.84378, 0.50401, 0.73483],[0.60764, 0.56052, 0.72926],[0.84627, 0.51498, 0.31315],[0.67958, 0.63586, 0.32038],[0.55914, 0.122, 0.54453],[0.60662, 0.63218, 0.97305],[0.05872, 0.99922, 0.72578],[0.11835, 0.93959, 0.17565],[0.3534, 0.77057, 0.28737],[0.82055, 0.99071, 0.02374],[0.9913, 0.88559, 0.02315],[0.98701, 0.5556, 0.02744],[0.0, 0.0, 0.96],[0.99042, 0.02403, 0.49195],[0.98367, 0.03078, 0.83615],[0.75325, 0.01445, 1.0],[0.44315, 0.01663, 0.99782],[0.1939, 0.02374, 0.99071],[0.02837, 0.25876, 0.98608],[0.28688, 0.45071, 0.23043],[0.0, 0.0, 0.88],[0.15323, 0.99165, 0.95836],[0.15097, 0.99391, 0.71032],[0.70704, 0.70552, 0.3509],[0.71952, 0.60694, 0.33841],[0.55616, 0.54257, 0.50178],[0.70294, 0.69401, 0.55789],[0.78703, 0.69512, 0.47379],[0.78975, 0.81033, 0.45049],[0.79997, 0.77511, 0.75068],[0.99628, 0.70149, 0.22106],[0.8294, 0.72125, 0.79823],[0.58798, 0.53854, 0.42649],[0.32386, 0.32592, 0.35729],[0.82428, 0.18732, 0.97211],[0.0, 0.0, 1.0],[0.0, 0.0, 1.0],[1.0, 1.0, 0.0],[0.0, 0.0, 0.0],[0.42959, 0.66659, 0.34786],[0.39344, 0.62101, 0.45034],[0.14893, 0.99596, 0.47106],[0.16101, 0.98387, 0.20855],[0.47774, 0.63362, 0.66714],[0.3, 0.3, 0.3],[0.3, 0.3, 0.3],[0.3, 0.3, 0.3],[0.3, 0.3, 0.3]]

# element negativity message from OQMD
element_negativity = [2.2, 0.0, 0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 0.0, 0.93, 1.31, 1.61, 1.9, 2.19, 2.58, 3.16, 0.0, 0.82, 1.0, 1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.9, 1.65, 1.81, 2.01, 2.18, 2.55, 2.96, 3.0, 0.82, 0.95, 1.22, 1.33, 1.6, 2.16, 1.9, 2.2, 2.28, 2.2, 1.93, 1.69, 1.78, 1.96, 2.05, 2.1, 2.66, 2.6, 0.79, 0.89, 1.1, 1.12, 1.13, 1.14, 1.13, 1.17, 1.2, 1.2, 1.2, 1.22, 1.23, 1.24, 1.25,1.1, 1.27, 1.3, 1.5, 2.36, 1.9, 2.2, 2.2, 2.28, 2.54, 2.0, 1.62, 2.33, 2.02, 2.0, 2.2, 2.2, 0.7, 0.9, 1.1, 1.3, 1.5, 1.38, 1.36, 1.28, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,]

# vesta bond setting
with open('/data/control/style.ini','r') as tmp:
    bond_list_ini = tmp.readlines()

class element(object):
    def __init__(self,index):
        if type(index) == int:
            self.index = index
            self.element = elementtable[index - 1]
            self.radius = element_radious[index - 1]
            self.valence = element_valence[index - 1]
            self.color = element_color[index - 1]
        elif index in elementtable:
            index = elementtable.index(index) + 1
            self.index = index
            self.element = elementtable[index - 1]
            self.radius = element_radious[index - 1]
            self.valence = element_valence[index - 1]
            self.color = element_color[index - 1]


# create bond list message

bond_list = {}
for i in range(13,927):
    file_context = bond_list_ini[i].split()
    p1 = file_context[1]
    if p1 in elementtable:
        index1 = elementtable.index(p1) + 1
    p2 = file_context[2]
    if p2 in elementtable:
        index2 = elementtable.index(p2) + 1
    distance = file_context[4]
    bond_pair = [index1,index2]
    bond_list[(index1,index2)] = float(distance)

# calculation realposition function

def calc_coordinate(position, vector):
    new_position = []
    for i in range(3):
        x = 0
        for j in range(3):
            x += float(position[j]) * float(vector[j][i])
        new_position.append(x)
    return new_position

# translattion function
def translation(positionA, positionB):
    return [positionA[i] + positionB[i] for i in range(3)]

def translation_del(positionA, positionB):
    return [positionA[i] - positionB[i] for i in range(3)]

# define ternary expansion
def resolve_3(num):
    tmp = [-1] * 3
    index = 0
    while(num > 2):
        tmp[index] = (num % 3) - 1
        num = num // 3
        index += 1
    tmp[index] = num - 1
    return tmp

# define function calculation distance between A and B :
def calc_distance(atomA, atomB = [0,0,0]):
    distance = 0
    for i in range(3):
        distance += (atomA[i] - atomB[i]) ** 2
    return distance ** 0.5

# judge the position in the crystal:
def judge_vector(position):
    for i in position:
        if i < 0 or i > 1:
            return False
    return True

# construct the function of two vector:
def vdot(vector1, vector2):
    product = 0
    for i in range(3):
        product += vector1[i] * vector2[i]
    return product

def vcross(vector1, vector2):
    x = vector1[1] * vector2[2] - vector1[2] * vector2[1]
    y = vector1[2] * vector2[0] - vector1[0] * vector2[2]
    z = vector1[0] * vector2[1] - vector1[1] * vector2[0]
    return [x, y, z]

def conclude_cover(i, j, k, index_list):
    a = translation_del(pertubation_list[j]['position'], pertubation_list[i]['position'])
    b = translation_del(pertubation_list[k]['position'], pertubation_list[i]['position'])
    normal_vector = vcross(a,b)
    plane = [i, j, k]
    state = 0
    for index in index_list:
        if index not in [i, j, k]:
            c = translation_del(pertubation_list[index]['position'], pertubation_list[i]['position'])
            result = vdot(c, normal_vector)
            if state != 0 and state * result <0:
                return False
            elif result != 0:
                state = result
            else:
                plane.append(index)
    return plane

# read the log
with open('/data/control/add_position.log','r') as file:
    path_list_old = [i.replace('\n','') for i in file.readlines()]

path_list_total = os.listdir('/data/datafile/')

path_list = list(set(path_list_total) ^ set(path_list_old))

for path in path_list:
    # only one atom.config read
    # atom.config,POSCAR
    print(path)
    os.chdir('/data/datafile/{0}'.format(path))
    with open('/data/datafile/{0}/atom.config'.format(path),'r') as tmp:
        atom_config = tmp.readlines()

    atom_amount = int(atom_config[0].split()[0])
    lattice_vector = [[float(x) for x in atom_config[i+2].split()] for i in range(3)]

    # define scale ceofficient
    lattice_vector_length = [calc_distance(i) for i in lattice_vector]
    alpha = math.acos(vdot(lattice_vector[1],lattice_vector[2]) / lattice_vector_length[1] / lattice_vector_length[2]) * 180 / math.pi
    beta = math.acos(vdot(lattice_vector[2],lattice_vector[0]) / lattice_vector_length[2] / lattice_vector_length[0]) * 180 / math.pi
    gama = math.acos(vdot(lattice_vector[0],lattice_vector[1]) / lattice_vector_length[0] / lattice_vector_length[1]) * 180 / math.pi
    angle = [alpha, beta, gama]
    max_lattice_vector = max(lattice_vector_length)
    scale_position = 500 / max_lattice_vector
    scale_bond_radius = scale_position / 9
    scale_atom_radius = scale_position / 3

    # claculation the atom's message
    atom_list = []
    atom_species = []
    for i in range(atom_amount):
        relative_position = atom_config[6+i].split()
        unit_atom = {}
        unit_atom['index'] = int(relative_position[0])
        atom_species.append(int(relative_position[0]))
        unit_atom['relative_position'] = [float(x) for x in relative_position[1:4]]
        unit_atom['position'] = calc_coordinate(unit_atom['relative_position'], lattice_vector)
        unit_atom['data'] = 'initial'
        atom_list.append(unit_atom)

    atom_kind = list(set(atom_species))
    # extend the crystal structure
    extend_atom_list = []

    for i in range(len(atom_list)):
        atom = atom_list[i]
        for k in range(27):
            if k != 13:
                translation_ceof = resolve_3(k)
                unit_atom = {}
                unit_atom['index'] = atom['index']
                unit_atom['relative_position'] = translation(atom['relative_position'],translation_ceof)
                unit_atom['position'] = calc_coordinate(unit_atom['relative_position'], lattice_vector)
                if judge_vector(unit_atom['relative_position']):
                    unit_atom['data'] = 'extend'
                    atom_list.append(unit_atom)
                else:
                    unit_atom['data'] = 'need'
                    extend_atom_list.append(unit_atom)

    #calculation the initial crystal's real position
    total_atom_amount = len(atom_list)
    super_crystal = atom_list + extend_atom_list

    # calculation bond list:
    crystal_bond_list = []
    polyhedra = {}
    for i in range(len(atom_list)):
        atomA = atom_list[i]
        for j in range(len(super_crystal)):
            atomB = super_crystal[j]
            bond_key = (atomA['index'],atomB['index'])
            if (bond_key in bond_list) and (i != j):
                bond_length = calc_distance(atomA['position'],atomB['position'])
                if bond_length <= bond_list[bond_key]:
                    if atomB in atom_list:
                        extend_index = atom_list.index(atomB)
                    else:
                        atom_list.append(atomB)
                        extend_index = atom_list.index(atomB)
                    crystal_bond_list.append([i,extend_index])
                    if i in polyhedra:
                        polyhedra[i].append(extend_index)
                    else:
                        polyhedra[i] = [extend_index]

    # position-data:atom_list[index]['position']
    # construct the polyhedra data:

    # pertubation aton list
    pertubation_list = atom_list[:]
    for i in pertubation_list:
        pertubation = [0.05 * random.rand() for j in range(3)]
        i['position'] = translation(i['position'], pertubation)

    # don't care len(plane) > 3, we add a small perturbation to the crystal
    surface_list = {}
    for center_index in polyhedra:
        cover_list = polyhedra[center_index]
        num = len(cover_list)
        faces = []
        for i in range(num - 2):
            for j in range(i+1, num - 1):
                for k in range(j+1, num):
                    plane = conclude_cover(cover_list[i], cover_list[j], cover_list[k], cover_list)
                    if plane:
                        faces.append(plane)
        surface_list[center_index] = faces
    # define the vector of translation:
    a = lattice_vector[0]
    b = lattice_vector[1]
    c = lattice_vector[2]
    total_vector = translation(translation(a,b),c)
    translation_vector = [x * scale_position/2 for x in total_vector]

    with open('position.js', 'w') as f:
        f.write('var element = new Array();\n')
        for index in atom_kind:
            atom = element(index)
            f.write('element[{0}] = {1}\n'.format(atom.index,'{'))
            f.write("\telement: '{0}',\n".format(atom.element))
            f.write('\tradius: {0:.6f},\n'.format(atom.radius * scale_atom_radius))
            f.write('\tcolor: [{0}, {1}, {2}],\n{3};\n'.format(atom.color[0], atom.color[1], atom.color[2], '}'))
        f.write('var vector_a_ini = [{0:.6f}, {1:.6f}, {2:.6f}],\n'.format(a[0],a[1],a[2]))
        f.write('vector_b_ini = [{0:.6f}, {1:.6f}, {2:.6f}],\n'.format(b[0],b[1],b[2]))
        f.write('vector_c_ini = [{0:.6f}, {1:.6f}, {2:.6f}],\n'.format(c[0],c[1],c[2]))
        f.write('vector_a = [{0:.6f}, {1:.6f}, {2:.6f}],\n'.format(a[0]*scale_position,a[1]*scale_position,a[2]*scale_position))
        f.write('vector_b = [{0:.6f}, {1:.6f}, {2:.6f}],\n'.format(b[0]*scale_position,b[1]*scale_position,b[2]*scale_position))
        f.write('vector_c = [{0:.6f}, {1:.6f}, {2:.6f}];\n'.format(c[0]*scale_position,c[1]*scale_position,c[2]*scale_position))
        f.write('var atom_list = [\n')
        for atom in atom_list:
            x = atom['position'][0] * scale_position - translation_vector[0]
            y = atom['position'][1] * scale_position - translation_vector[1]
            z = atom['position'][2] * scale_position - translation_vector[2]
            f.write('\t[ {0}, [{1:.6f}, {2:.6f}, {3:.6f}] ],\n'.format(atom['index'],x,y,z))
        f.write(']\nvar bond_list = [\n')
        for bond in crystal_bond_list:
            f.write('[{0}, {1}],'.format(bond[0],bond[1]))
        f.write('];\nvar face_list = [\n')
        for i in polyhedra:
            f.write('{0} {1}, {2} {3}, \n'.format('[',i,surface_list[i],']'))
        f.write('];\n')
        f.write('var bond_radius = {0:.6f};\n'.format(scale_bond_radius))
        f.write('var translation_vector = [{0:.6f}, {1:.6f}, {2:.6f}];\n'.format(translation_vector[0],translation_vector[1],translation_vector[2]))
        f.write('var angle = [{0:.2f}, {1:.2f}, {2:.2f}],\n'.format(angle[0],angle[1],angle[2]))
        f.write('vector_length = {};\n'.format(lattice_vector_length))
        f.write('var scale_position = {:.2f}'.format(scale_position))

    os.system('echo "{}" >> /data/control/add_position.log'.format(path))