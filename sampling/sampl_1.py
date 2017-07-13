from Bio.PDB import *
from Bio.SeqUtils import *
from Bio.PDB.Chain import Chain
import numpy as np
import math
import numbers
from Bio.SeqUtils import *
import math
from math import sin, cos, sqrt
from numpy import dot, linalg, arccos
import numbers
from utils.calc import *
from Bio.PDB.Superimposer import *

def imposer(_structure1,res11,res12):
    """
        Структурное выравнивание по 2 цепочкам атомов
        Параметры:
            structure1 - петелька
            res11,res12 - аминокислоты двух концов белка
    """
    first = 10000
    last = 0
    structure1 = _structure1.child_list
    #structure1 = _structure1.copy()
    for x in structure1:
        first = min(x.id[1],first)
        last = max(x.id[1],last)

    Nc = structure1[first-1]['N']
    Cc = structure1[first-1]['C']
    CAc = structure1[first-1]['CA']

    Ncf =  res11['N']
    Ccf =  res11['C']
    CAcf = res11['CA']

# IMPORTSER MODULE
    fixed_vectors = [Ncf,CAcf,Ccf]
    moving_vectors = [Nc,CAc,Cc]
    sup = Superimposer()
    sup.set_atoms(fixed_vectors,moving_vectors)
    for spx in range(first,last):
        for i in structure1[spx]:
            v = structure1[spx][i.get_name()].get_vector()
            curcord = dot(v._ar, sup.rotran[0])+sup.rotran[1]
            structure1[spx][i.get_name()].set_coord(curcord)
# 1 SELECT RES ALIGN
    _res11 = res11.copy()
    __res11 = res11.copy()
    _res11.__init__(structure1[first-1].id,__res11.resname,__res11.segid)
    structure1.__delitem__(first-1)
    structure1.insert(first-1,_res11)
    for x in res11:
        structure1[first-1].add(x)


    Nc = structure1[last-1]['N']
    Cc = structure1[last-1]['C']
    CAc = structure1[last-1]['CA']
    _res12 = res12.copy()
    __res12 = res12.copy()
    _res12.__init__(structure1[last-1].id,__res12.resname,__res12.segid)
    structure1.remove(structure1[last-1])
    structure1.insert(last,_res12)
    for x in res12:
        structure1[last-1].add(x)

    Ncf = res12['N']
    Ccf = res12['C']
    CAcf = res12['CA']


    fixed_vectors = [Nc,CAc,Cc]
    moving_vectors = [Ncf,CAcf,Ccf]
    sup = Superimposer()
    sup.set_atoms(fixed_vectors,moving_vectors)
    for i in structure1[last-1]:
        v = structure1[last-1][i.get_name()].get_vector()
        curcord = dot(v._ar, sup.rotran[0])+sup.rotran[1]
        structure1[last-1][i.get_name()].set_coord(curcord)

    return structure1
def rmsd(structure1,structure2):
    """
        Считает rmsd между двумя структурами
        Параметры:
            structure1,structure2 - цепочки атомов
    """
    fixed = []

    for x in structure1:
        fixed.append(x['CA'].get_vector())
    moving = []
    for x in structure2:
        moving.append(x['CA'].get_vector())
    # sup = Superimposer()
    # sup.set_atoms(fixed,moving)
    # return sup.rms

    N = len(moving)
    preval = 0
    for x in range(0,N):
        dist = distance(fixed[x],moving[x])
        preval += dist**2
    value = sqrt(preval/N)
    return value
def normalize(v):
    """
        Нормализует вектор v
        Параметры:
            v - ветор
    """
    n = linalg.norm(v)
    if n == 0:
        return v
    return div(v, linalg.norm(v))
def shift(structure,value,atomA = 'C',atomB = 'N'):
    """
        Увеличивает расстояние между атомами C N
        Параметры:
             structure - цепочка атомов
             atomA,atomB - атомы между которыми будет произведено изменение длины
    """
    for j in range(1,len(structure)):
        cordC = structure[j]['C'].get_vector()
        cordN = structure[j+1]['N'].get_vector()
        CN = cordC-cordN
        CN.normalize()
        CN[0]*=value
        CN[1]*=value
        CN[2]*=value
        Ni =  cordC - CN
        NNi = Ni - cordN
        for spx in range(j+1,len(structure)):
            for i in structure[spx]:
                    vect = structure[spx][i.get_name()].get_vector()+NNi
                    structure[spx][i.get_name()].set_coord(vect)
    return structure
def calcnewcord(_veca,_vecb,_vecc,angle):
    """
        Расчёт нового вектора, при повороте на заданный угол вокруг оси
        Параметры:
            _veca,vecb - векторы оси вектора PDBшные
            _vecc -  начальный вектор
            angle - угол поворота (в радианах)
    """
    veca=_veca.get_array()
    vecb=_vecb.get_array()
    vecc=_vecc.get_array()
    Eab = vecb - veca
    Eab = normalize(Eab)
    AC = vecc - veca
    O = Eab.copy()
    O[0] *= dot(AC, Eab)
    O[1] *= dot(AC, Eab)
    O[2] *= dot(AC, Eab)
    O += veca
    Eoc = normalize(vecc - O)
    S = vec_mult_vec(Eab,Eoc)
    OCd = distance(O,vecc)
    M = Eoc * OCd * math.cos(angle) + S * OCd * math.sin(angle) + O
    return M
def quat_invert(quat):
    '''Инверсия кватерниона'''
    res = quat * -1
    res[3] = quat[3]
    return res / np.linalg.norm(res)
def quat_mul_quat(quat1, quat2):
    """Произведение кватернионов"""
    res = vec_mult_vec(np.resize(quat1, 3), np.resize(quat2, 3)) + np.resize(quat1, 3) * quat2[3] + np.resize(quat2,
                                                                                                                  3) * \
                                                                                                        quat1[3]
    res = np.resize(res, 4)
    res[3] = quat1[3] * quat2[3] - np.dot(quat1, quat2)
    return res
def dist(coordA, coordB):
    # длина вектора
    return sqrt(pow(coordA[0] - coordB[0], 2) + pow(coordA[1] - coordB[1], 2) + pow(coordA[2] - coordB[2], 2))
def vec_mult_vec(vec1, vec2):
    '''Векторное произведение векторов'''
    vec = np.zeros(3)
    vec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1]
    vec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2]
    vec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0]
    return vec
def rotate_vector(point1, point2, vec, angle):
    '''Вращение вектора на заданный угол вокруг оси, заданной двумя точками
    point1, point2, vec, angle --> vec
    V' = qvq^(–1)
    '''
    new_vec = np.resize(vec - point1, 4)
    new_vec[3] = 0
    quat = make_quaternion(point1, point2, angle)
    t = quat_mul_quat(quat, new_vec)
    t = quat_mul_quat(t, quat_invert(quat))
    return np.resize(t, 3) + point1
def make_quaternion(point1, point2, angle):
    '''
    Создание кватерниона, заданного осью вращения и углом
    point1, point2, angle --> quat
    q = [cos(a/2), [x sin(a/2), y sin(a/2), z sin(a/2)]]
    '''
    startpoint = point2 - point1

    rotate_vector = startpoint / np.linalg.norm(startpoint)
    quat = np.zeros(4)
    quat[0] = rotate_vector[0] * sin(angle / 2)
    quat[1] = rotate_vector[1] * sin(angle / 2)
    quat[2] = rotate_vector[2] * sin(angle / 2)
    quat[3] = cos(angle / 2)
    return quat
def movechain(structure, coord, extra=0):
    for spx in range(extra, len(structure) + extra):
        for i in structure[spx]:
            vect = structure[spx][i.get_name()].get_vector() + coord
            structure[spx][i.get_name()].set_coord(vect)
    return structure
def getletter(structure):
    """
        Достаёт букву из структуры
        Параметры:
            structure - структура
    """
    return structure[0].child_list[0].id
def generate(s):
    structure = read('/home/ludmila/git/sirius2017_protein/sampling/aminos_out.pdb', "test")
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
    return [result]
def div(v, k):
    v[0]/=k
    v[1]/=k
    v[2]/=k
    return v
def grad(angle):
    return angle*180/pi
def rot(amino,angle, O, C, N, H):
    """
        Поворот структуры вогрук оси С-N
        Параметры:
            structure - структура
            angle - угол поворота [pi,pi]
    """ # -2 из за заглушки
    DIangle = angle-calc_dihedral(O,C,N,H)
    extra =  DIangle
    for atom in amino:
        vecALL = atom.get_vector()
        curcord = calcnewcord(C, N, vecALL, extra)
        atom.set_coord(curcord)
    return amino
def rad(angle):
    return angle*pi/180
"""def rot_in_plane(A, B, C, D, angle, amino):
    выставляет угол ABD=angle в плоскости ABC
    perp = vec_mult_vec(A - B, C - B)
    coord=[]
    coord[0] = perp[0]
    coord[1] = perp[1]
    coord[2] = perp[2]
    angle -= calc_angle(A, B, D)
    for atom in amino:
        atom.set_coord(rotate_vector(B.get_array(), (coord+B).get_array(), atom.get_vector().get_array(), -angle))"""
def move_dist(amino, dist, vec, A, B):
    coord = A - B
    for atom in amino:
        atom.set_coord(atom.get_vector() + coord)
    coord = vec
    coord.normalize()
    coord[0] *= 1.4
    coord[1] *= 1.4
    coord[2] *= 1.4
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
    amino = rot(amino, pi, O, C, N, H)
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


def getletter(structure):
    """
        Достаёт букву из структуры
        Параметры:
            structure - структура
    """
    return structure[0].child_list[0].id


def read(path,name = "test"):
    """
        Читает цепочку аминокислот
        Параметры:
            path - путь к файлу
            name - имя структуры
    """
    parser = PDBParser()
    structure = parser.get_structure(name,path)
    return structure[0][getletter(structure)]

def write(path,record):
    """
        Записывает цепочку в файл
        Параметры:
            path - путь к файлу
            name - имя структуры
    """
    io = PDBIO()
    io.set_structure(record)
    io.save(path)

def addDash(s,index):
    """
        Вставляет тире в строке по индексу
        Параметры:
            s - строка
            index - индекс
    """
    return s[:index] + '-' + s[index:]

def tta(tup,extra = 0):
    a = [' ']
    for x in tup:
        a.append(x)
    a.append(' ')
    return a

def compileres(residues,id = 1):
    """
        Собирает список из residue в chain
        Параметры:
            residues - список residue
            id - идентификатор цепочки
    """
    chain = Chain(id)
    startid = 0
    for x in residues:
        x.detach_parent()
        x.id = (' ',startid,' ')
        startid+=1
        chain.add(x)
    return chain

def writeres(path,residues,id = 1):
    #print(type(compileres(residues)))
    write(path,compileres(residues,id))

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
            residues[j] = rot(residues[j], rad(angles[i][0]), N, CA, C, N2)
        O = residues[i]['O'].get_vector()
        resi_copy['O'].set_coord(O)
        residues[i]=resi_copy
        CA = residues[i]['CA'].get_vector()
        C = residues[i]['C'].get_vector()
        CA2 = residues[i+1]['CA'].get_vector()
        N2 = residues[i+1]['N'].get_vector()
        for j in range(i + 1, len(residues)):
            residues[j] = rot(residues[j], rad(angles[i][1]), CA, C, N2, CA2)
        C = residues[i]['C'].get_vector()
        CA2 = residues[i+1]['CA'].get_vector()
        N2 = residues[i+1]['N'].get_vector()
        C2 = residues[i+1]['C'].get_vector()
        H = get_H(residues[i+1])
        for j in range(i + 1, len(residues)):
            residues[j] = rot(residues[j], rad(angles[i][2]), C, N2, CA2, C2)
        if(residues[i+1].get_resname()!="PRO"):
            residues[i+1]['H'].set_coord(H)
        """for i in range(len(residues)-1):
            C = residues[i]['C'].get_vector()
            O = residues[i]['O'].get_vector()
            N = residues[i+1]['N'].get_vector()
            H = residues[i+1]['H'].get_vector()
            perp = vec_mult_vec(C - N, O - N)
            angle = rad(120)-calc_angle(C, N, H)
            perp[0] += N[0]
            perp[1] += N[1]
            perp[2] += N[2]
            residues[i+1]['H'].set_coord(rotate_vector(N.get_array(), perp, residues[i+1]['H'].get_vector().get_array(), -angle))"""
    return residues
def main():
    s=generate("AQGP")
    angles=[(-53.005436, 166.00095, -107.3551), (137.84758, -172.47937, 90.52946), (-83.336555, -6.25399, -71.49093)]
    for i in range(len(s)):
        cur_s = s[i]
        writeres('/tmp/t' + str(i) + '.pdb', cur_s)
        result = rotate_by_dih(cur_s, angles)
        writeres('/tmp/'+str(i)+'.pdb', result)


if __name__ == "__main__":
    main()

    




