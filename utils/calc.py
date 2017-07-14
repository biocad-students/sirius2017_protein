from numpy import *
from math import *
from Bio.PDB.Superimposer import Superimposer

def distance(cordA,cordB):
    """
         Возвращает расстояние между точками в пространстве
         Параметры:
            cordA, cordB - координаты точек
    """
    return sqrt((cordA[0]-cordB[0])**2+(cordA[1]-cordB[1])**2+(cordA[2]-cordB[2])**2)

def normalize(v):
    """
        Нормализует вектор v
        Параметры:
            v - ветkор
    """
    n = linalg.norm(v)
    if n == 0:
        return v
    return v / n


def imposer(structure,_leftAmino,_rightAmino):
    """
        Структурное выравнивание по 2 цепочкам атомов
        Параметры:
            structure1 - петелька
            leftAmino,rightAmino - аминокислоты двух концов белка
    """
    leftAmino = _leftAmino.copy()
    rightAmino = _rightAmino.copy()
    leftAminoCDR = structure[0].copy()
    rightAminoCDR = structure[-1].copy()
    structure1 = [x.copy() for x in structure[1:-1]]

    Nc = leftAminoCDR['N']
    Cc = leftAminoCDR['C']
    CAc = leftAminoCDR['CA']

    Ncf =  leftAmino['N']
    Ccf =  leftAmino['C']
    CAcf = leftAmino['CA']

    fixed_vectors = [Ncf,CAcf,Ccf]
    moving_vectors = [Nc,CAc,Cc]
    sup = Superimposer()
    sup.set_atoms(fixed_vectors,moving_vectors)
    (rotation_m, transition_v) = sup.rotran
    for amino in (structure1 + [rightAminoCDR]):
        for atom in amino:
            v = atom.get_vector()
            curcord = dot(v._ar, rotation_m) + transition_v
            atom.set_coord(curcord)

    Nc = rightAminoCDR['N']
    Cc = rightAminoCDR['C']
    CAc = rightAminoCDR['CA']

    Ncf = rightAmino['N']
    Ccf = rightAmino['C']
    CAcf = rightAmino['CA']

    fixed_vectors = [Nc,CAc,Cc]
    moving_vectors = [Ncf,CAcf,Ccf]
    sup = Superimposer()
    sup.set_atoms(fixed_vectors,moving_vectors)
    (rotation_m, transition_v) = sup.rotran
    for atom in rightAmino:
        v = atom.get_vector()
        curcord = dot(v._ar, rotation_m) + transition_v
        atom.set_coord(curcord)

    result = [leftAmino] + structure1 + [rightAmino]
    return result

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

    N = len(moving)
    preval = 0
    for x in range(0,N):
        dist = distance(fixed[x],moving[x])
        preval += dist**2
    value = sqrt(preval/N)
    return value


def link(structure1,structure2):
    """
        Считает матрицу поворота между двумя аминокислотами
        Параметры:
            structure1 - стоит
            structure2 - двигается
    """
    Cf = structure1['C']
    CAf = structure1['CA']
    Nf = structure1['N']

    Cm = structure1['C']
    CAm = structure1['CA'].get_vector()
    Nm = structure1['N'].get_vector()

    fixed = [Cf,CAf,Nf]
    moving = [Cm,CAm,Nm]
    sup = Superimposer()
    sup.set_atoms(fixed,moving)
    return sup.rotran
