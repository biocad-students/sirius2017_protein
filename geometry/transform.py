from Bio.PDB import *
from utils.calc import *
from numpy import *

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

def rot(structure,angle,start = 0,isNewCalcReq = 1,FA = 'CA',FB = 'C',SA = 'N',SB = 'CA',):
    """
        Поворот структуры вогрук оси С-N
        Параметры:
            structure - структура
            angle - угол поворота радианы
            start - номер аминокислоты с которой начинать поворот
            stop - сдвиг конца
            FA,FB - атомы первой аминокислоты
            SA,SB - атомы второй аминокислоты
            isNewCalcReq - вычисление дигидрального угла

    """
    for j in range(start,len(structure)-1): # -2 из за заглушки
        vecCA = structure[j][FA].get_vector()
        vecC = structure[j][FB].get_vector()
        vecNN = structure[j+1][SA].get_vector()
        vecCAN = structure[j+1][SB].get_vector()
        if(isNewCalcReq == 1):
            DIangle = calc_dihedral(vecCA,vecC,vecNN,vecCAN)
            extra = angle - DIangle
        else:
            extra = angle
        for spx in range(j+1,len(structure)):
            for i in structure[spx]:
                vecALL = structure[spx][i.get_name()].get_vector()
                curcord = calcnewcord(vecC,vecNN,vecALL, extra)
                structure[spx][i.get_name()].set_coord(curcord)
    return structure
