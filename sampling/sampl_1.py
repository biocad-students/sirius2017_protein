from Bio.PDB import *
from Bio.SeqUtils import *
from Bio.PDB.Chain import Chain
import numpy as np
import math
from math import sin, cos, sqrt
from numpy import dot, linalg, arccos
import numbers
from Bio.SeqUtils import *
import math
from math import sin, cos, sqrt
from numpy import dot, linalg, arccos
import numbers
from utils.calc import *

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
    veca = _veca.get_array()
    vecb = _vecb.get_array()
    vecc = _vecc.get_array()

    Eab = normalize(vecb - veca)
    AC = vecc - veca
    O = veca + Eab * dot(AC,Eab)
    Eoc = normalize(vecc - O)
    S = cross(Eab,Eoc)
    OCd = distance(O,vecc)
    M = Eoc * OCd * math.cos(angle) + S * OCd * math.sin(angle) + O
    return M
def rot(structure,point1, point2, angle):
    """
        Поворот структуры вогрук оси С-N
        Параметры:
            structure - структура
            angle - угол поворота [pi,pi]
    """
    for j in range(2,len(structure)): # -2 из за заглушки
        vecCA = structure[j]['CA'].get_vector()
        vecCAN = structure[j+1]['CA'].get_vector()
        DIangle = calc_dihedral(vecCA,point1,point2,vecCAN)
        extra = angle - DIangle
        for spx in range(j+1,len(structure)):
            for i in structure[spx]:
                vecALL = structure[spx][i.get_name()].get_vector()
                curcord = calcnewcord(point1,point2,vecALL, extra)
                structure[spx][i.get_name()].set_coord(curcord)
    return structure
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
def read(path, name):
    parser = PDBParser()
    structure = parser.get_structure(name, path)
    return structure[0][getletter(structure)]
def generate(s, cnt):
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
def moveTo(last, amino):
    N = amino['N'].get_vector()
    H=amino['H'].get_vector()
    CA2 = amino['CA'].get_vector()
    CA = last['CA'].get_vector()
    O = last['O'].get_vector()
    C = last['C'].get_vector()
    coord = C - N
    for atom in amino:
        atom.set_coord(atom.get_vector()+coord)
    coord = C - CA
    coord.normalize()
    coord[0] *= 1.4
    coord[1] *= 1.4
    coord[2] *= 1.4
    for atom in amino:
        atom.set_coord(atom.get_vector()+coord) #до сюда выставляется расстояние СN=1.4, дигидральный угол CA,C,O,N=180
    NH = dist(N, H)
    ON = N - O
    coord = N - H
    amino['H'].set_coord(atom.get_vector() + coord)
    ON.normalize()
    ON[0] *= NH
    ON[1] *= NH
    ON[2] *= NH
    coord=ON
    amino['H'].set_coord(atom.get_vector() + coord)
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
    for x in residues:
        chain.add(x)
    return chain
def culc_angle(a, b, c):
    x = dist(a, b)
    y = dist(b, c)
    z = dist(a, c)
    return arccos((x * x + y * y - z * z) / (2 * x * y))

def writeres(path,residues,id = 1):
    #print(type(compileres(residues)))
    write(path,compileres(residues,id))
def main():
    s=generate("AQG", 10)
    for i in range(len(s)):
        writeres('/tmp/'+str(i)+'.pdb', s[i])

if __name__ == "__main__":
    main()

    




