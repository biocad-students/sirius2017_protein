from numpy import linalg,multiply
from math import sqrt
from Bio.PDB.Superimposer import *
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
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

#

def imposer(structure1,res11,res12):
    """
        Структурное выравнивание по 2 цепочкам атомов
        Параметры:
            structure1 - петелька
            res11,res12 - аминокислоты двух концов белка
    """
    first = 10000
    last = 0

    for x in structure1:
        first = min(x.id[1],first)
        last = max(x.id[1],last)
    print(res11,res12,len(structure1))
    print(first,last)
    res11.id = structure1.child_list[0].id
    res12.id = structure1.child_list[last].id
    structure1.child_list[0] = res11
    structure1.child_list[last-1] = res12
    Nc = structure1[first]['N']
    Cc = structure1[first]['C']
    CAc = structure1[first]['CA']
    Ncf = res11['N']
    Ccf = res11['C']
    CAcf = res11['CA']
    fixed_vectors = [Ncf,Ccf,CAcf]
    moving_vectors = [Nc,Cc,CAc]
    sup = Superimposer()
    sup.set_atoms(fixed_vectors,moving_vectors)
    for spx in range(first,last):
        for i in structure1[spx]:
            vecALL = structure1[spx][i.get_name()].get_vector()
            curcord = vecALL._ar.dot(sup.rotran[0])+sup.rotran[1]
            #print(vecALL,curcord)
            structure1[spx][i.get_name()].set_coord(curcord)
    return structure1
    #return sup
