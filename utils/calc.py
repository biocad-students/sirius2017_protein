from numpy import linalg,dot
from math import sqrt
from Bio.PDB.Superimposer import *
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from utils.io import tta
def distance(cordA,cordB):
    """
         Возвращает расстояние между точками в пространстве
         Параметры:
            cordA, cordB - координаты точек
    """
    return sqrt(pow(cordA[0]-cordB[0],2)+pow(cordA[1]-cordB[1],2)+pow(cordA[2]-cordB[2],2))

def normalize(v):
    """
        Нормализует вектор v
        Параметры:
            v - ветор
    """
    n = linalg.norm(v)
    if n == 0:
        return v
    return v / linalg.norm(v)


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
