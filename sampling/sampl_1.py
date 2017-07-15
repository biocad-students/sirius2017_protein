from Bio.SeqUtils import *
from geometry.rotation import *
from geometry.transform import *
from utils.io import *
from params.params_sr3 import *

def dist(coordA, coordB):
    # длина вектора
    return sqrt(pow(coordA[0] - coordB[0], 2) + pow(coordA[1] - coordB[1], 2) + pow(coordA[2] - coordB[2], 2))

def generate(s):
    structure = read('../sirius_out/aminos_out.pdb', "test")
    arr = {}
    for residue in structure:
        arr[seq1(residue.get_resname())] = residue
    result=[]
    for c in s:
        if len(result) != 0:
            last = result[-1]
            res = moveTo(last,arr[c])
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
    perp = vec_mult_vec(A - B, C - B)
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

def moveTo(last, amino1):
    amino = amino1.copy()
    N = amino['N'].get_vector()
    H=get_H(amino)
    CA2 = amino['CA'].get_vector()
    CA = last['CA'].get_vector()
    O = last['O'].get_vector()
    C = last['C'].get_vector()
    C2 = amino['C'].get_vector()
    amino=move_dist(amino, 1.4, C-CA, C, N)
    #до сюда выставляется расстояние СN=1.4, дигидральный угол CA,C,O,N=180
    N = amino['N'].get_vector()
    CA = last['CA'].get_vector()
    O = last['O'].get_vector()
    C = last['C'].get_vector()
    amino = rotate_amino(amino, O, C, CA, N, rad(123))
    #выставился угол 123 градуса
    N = amino['N'].get_vector()
    O = last['O'].get_vector()
    C = last['C'].get_vector()
    H = get_H(amino)
    amino = rot1(amino, pi, O, C, N, H)
    #выставилась плоскость OCNH=180
    N = amino['N'].get_vector()
    H = get_H(amino)
    O = last['O'].get_vector()
    C = last['C'].get_vector()
    amino = rotate_amino(amino, C, N, O, H, rad(120))
    limit=2
    N = amino['N'].get_vector()
    H = get_H(amino)
    CA2 = amino['CA'].get_vector()
    if(dist(C, CA2)<limit):
        for atom in amino:
            atom.set_coord(rotate_vector(N.get_array(), H.get_array(), atom.get_vector().get_array(), pi))
    return amino

def culc_angle(a, b, c):
    x = dist(a, b)
    y = dist(b, c)
    z = dist(a, c)
    return arccos((x * x + y * y - z * z) / (2 * x * y))

def rotate_by_dih(residues, angles):
    for i in range(len(angles)):
        CA = residues[i]['CA'].get_vector()
        C = residues[i]['C'].get_vector()
        N = residues[i]['N'].get_vector()
        N2 = residues[i+1]['N'].get_vector()
        resi_copy=residues[i].copy()
        for j in range(i, len(residues)):
            residues[j] = rot1(residues[j], rad(angles[i][0]), N, CA, C, N2)
        O = residues[i]['O'].get_vector()
        resi_copy['O'].set_coord(O)
        residues[i]=resi_copy
        CA = residues[i]['CA'].get_vector()
        C = residues[i]['C'].get_vector()
        CA2 = residues[i+1]['CA'].get_vector()
        N2 = residues[i+1]['N'].get_vector()
        for j in range(i + 1, len(residues)):
            residues[j] = rot1(residues[j], rad(angles[i][1]), CA, C, N2, CA2)
        C = residues[i]['C'].get_vector()
        CA2 = residues[i+1]['CA'].get_vector()
        N2 = residues[i+1]['N'].get_vector()
        C2 = residues[i+1]['C'].get_vector()
        H = get_H(residues[i+1])
        for j in range(i + 1, len(residues)):
            residues[j] = rot1(residues[j], rad(angles[i][2]), C, N2, CA2, C2)
        if(residues[i+1].get_resname()!="PRO"):
            residues[i+1]['H'].set_coord(H)
    return residues

def samples(amino_seq, cnt):
    angles=varrand(amino_seq, cnt)
    result=[]
    for angle in angles:
        s=generate(amino_seq)
        result.append(rotate_by_dih(s, angle))
    return result

"""if __name__ == "__main__":
    main()"""
