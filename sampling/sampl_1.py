from Bio.SeqUtils import *
from geometry.transform import *
from utils.io import *
from utils.calc import distance
from params.params_sr3 import *

def dist(coordA, coordB):
    # длина вектора
    return sqrt(pow(coordA[0] - coordB[0], 2) + pow(coordA[1] - coordB[1], 2) + pow(coordA[2] - coordB[2], 2))

def generate(s, angles):
    structure = read('../files/aminos_out.pdb', "test")
    arr = {}
    for residue in structure:
        arr[seq1(residue.get_resname())] = residue
    result=[]
    for i in range(len(s)):
        c=s[i]
        if len(result) != 0:
            last = result[-1]
            res = moveTo(last,arr[c], angles[i-1])
            result.append(res)
        else:
            result.append(arr[c])
    return result

def grad(angle):
    return angle*180/pi

def rot1(amino1,angle, O, C, N, H):
    """
        Поворот структуры вогрук оси С-N
        Параметры:
            amino - структура
            angle - угол поворота [pi,pi]
    """ # -2 из за заглушки
    amino=amino1.copy()
    DIangle = angle-calc_dihedral(O,C,N,H)
    extra =  DIangle
    for atom in amino:
        vecALL = atom.get_vector()
        curcord = calcnewcord(C, N, vecALL, extra)
        atom.set_coord(curcord)
    return amino

def rad(angle):
    return angle*pi/180

def move_dist(amino, dis, vec, A, B):
    coord = A - B
    for atom in amino:
        atom.set_coord(atom.get_vector() + coord)
    coord = vec
    coord.normalize()
    coord[0] *= dis
    coord[1] *= dis
    coord[2] *= dis
    for atom in amino:
        atom.set_coord(atom.get_vector() + coord)
    return amino

def rotate_amino(amino, A, B, C, D, angle):
    """выставляет angle между A,B,D; C лежит в этой же плоскости"""
    perp = cross(A - B, C - B)
    angle -= calc_angle(A, B, D)
    perp[0] += B[0]
    perp[1] += B[1]
    perp[2] += B[2]
    for atom in amino:
        atom.set_coord(rotate_vector(B.get_array(), perp, atom.get_vector().get_array(), -angle))
    return amino

def get_H(amino):
    if(amino.get_resname()=="PRO"):
        return amino['CD'].get_vector()
    else:
        return amino['H'].get_vector()

def moveTo(last, amino1, angles):
    amino = amino1.copy()
    N = amino['N'].get_vector()
    CA = last['CA'].get_vector()
    C = last['C'].get_vector()
    amino=move_dist(amino, 1.4, C-CA, C, N)
    #до сюда выставляется расстояние СN=1.4, дигидральный угол CA,C,O,N=180
    N = amino['N'].get_vector()
    O = last['O'].get_vector()
    C = last['C'].get_vector()
    H = get_H(amino)
    amino = rot1(amino, pi, O, C, N, H)
    #выставилась плоскость OCNH=180
    H = get_H(amino)
    N = amino['N'].get_vector()
    CA2 = amino['CA'].get_vector()
    C = last['C'].get_vector()
    amino=rot1(amino, pi, C, N, H, CA2)
    #выставили 6 атомов в одну плоскость
    N = amino['N'].get_vector()
    CA = last['CA'].get_vector()
    O = last['O'].get_vector()
    C = last['C'].get_vector()
    amino = rotate_amino(amino, O, C, CA, N, rad(123))
    #выставился угол 123 градуса
    N = amino['N'].get_vector()
    H = get_H(amino)
    O = last['O'].get_vector()
    C = last['C'].get_vector()
    amino = rotate_amino(amino, C, N, O, H, rad(120))
    #выставляем 120 градусов
    CA = last['CA'].get_vector()
    C = last['C'].get_vector()
    N = last['N'].get_vector()
    N2 = amino['N'].get_vector()
    last_copy=last.copy()
    amino = rot1(amino, rad(angles[0]), N, CA, C, N2)
    last_copy = rot1(last_copy, rad(angles[0]), N, CA, C, N2)
    O = last_copy['O'].get_vector()
    last['O'].set_coord(O)
    #выставляем пси
    CA = last['CA'].get_vector()
    C = last['C'].get_vector()
    CA2 = amino['CA'].get_vector()
    N2 = amino['N'].get_vector()
    amino = rot1(amino, rad(angles[1]), CA, C, N2, CA2)
    #выставляем омегу
    C = last['C'].get_vector()
    CA2 = amino['CA'].get_vector()
    N2 = amino['N'].get_vector()
    C2 = amino['C'].get_vector()
    H = get_H(amino)
    amino = rot1(amino, rad(angles[2]), C, N2, CA2, C2)
    if (amino.get_resname() != "PRO"):
        amino['H'].set_coord(H)
    return amino

def culc_angle(a, b, c):
    x = distance(a, b)
    y = distance(b, c)
    z = distance(a, c)
    return arccos((x * x + y * y - z * z) / (2 * x * y))

def samples(amino_seq, cnt):
    angles=varrand(amino_seq, cnt)
    print(angles)
    result=[]
    for angle in angles:
        s=generate(amino_seq, angle)
        result.append(s)
    return result

"""if __name__ == "__main__":
    main()"""
