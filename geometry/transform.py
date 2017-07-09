from Bio.PDB import *

def move(structure,coord):
    """
    Сдвигает все атомы на cord
    Параметры:
        structure - цепочка
        coord - вектор сдвига
    """
    for spx in range(j,len(structure)):
        for i in structure[spx]:
                vect = structure[spx][i.get_name()].get_vector()+coord
                structure[spx][i.get_name()].set_coord(vect)
    return structure

def extend(structure,value):
    """
        Увеличивает расстояние между атомами C N
        Параметры:
             structure - цепочка атомов
    """
    for j in range(1,len(structure[0]['H'])):
        cordC = structure[0]['H'][j]['C'].get_vector()
        cordN = structure[0]['H'][j+1]['N'].get_vector()
        CN = cordC-cordN
        CN.normalize()
        CN[0]*=value
        CN[1]*=value
        CN[2]*=value
        Ni =  cordC - CN
        NNi = Ni - cordN
        print(NNi)
        for spx in range(j+1,len(structure)):
            for i in structure[spx]:
                    vect = structure[spx][i.get_name()].get_vector()+NNi
                    structure[spx][i.get_name()].set_coord(vect)
    return structure
