# *-* encoding: utf-8 *-*

import numpy as np
from math import sin, cos
from utils.calc import *
from numpy import *


cs = cross
def calcnewcord(veca,vecb,_vecc, sinangle, cosangle, Eab):
    """
        Вращение вектора на угол с заданным синусом и косинусом вокруг оси, заданной двумя точками
        Параметры:
            veca,vecb - векторы оси вектора 
            _vecc -  начальный вектор (PDB)
            sinangle - синус угла поворота
            cosangle - косинус угла поворота
    """
    vecc = _vecc.get_array()

    AC = vecc - veca
    O = veca + Eab * dot(AC,Eab)
    Eoc = normalize(vecc - O)
    S = cs(Eab,Eoc)
    OCd = distance(O,vecc)
    M = Eoc * OCd * cosangle + S * OCd * sinangle + O
    return M

