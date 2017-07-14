import numpy as np
from Bio.PDB import *
import math
from math import sin, cos, sqrt
from numpy import dot
import numbers

def read(path):
    parser = PDBParser()
    structure = parser.get_structure('test', '15477.pdb')
    return structure


def div(self, scalar):
    #Деление вектора на число: v / c
    assert isinstance(scalar, numbers.Number)
    return Vector(self[0] / scalar, self[1] / scalar, self[2]/scalar)

def write(path, name):
    io=PDBIO()
    io.set_structure(name)
    io.save('test.pdb')

"""def dist(v1, v2):
    dv1v2=math.sqrt((v1[0]-v2[0])**2+(v1[1]-v2[1])**2+(v1[2]-v2[2])**2)
    return dv1v2"""

def dist(self):
    #длина вектора
    return (self[0] ** 2 + self[1] ** 2+self[2]**2) ** 0.5

def normal(v):
    #единичный вектор
    Ev1v2=div(v, dist(v))
    return Ev1v2

def rotate(a, b, c, angle):
    ac=c-a
    dao=dot(ac, normal(b-a))
    ao=div(normal(b-c), 1/dao)
    o=ao+a
    s=np.cross(c-o, b-a)
    m=div(normal(c-o), 1/(dist(o-c)*cos(angle)))+div(s, 1/(dist(o-c)*sin(angle)))
    print(m)
def main():
    structure=read('15477.pdb')
    angle=90
    n = list(structure[0]['H'])[1]['N'].get_vector()
    ca = list(structure[0]['H'])[1]['CA'].get_vector()
    for atom in structure.get_atoms():
       # c = list(structure[0]['H'])[1]['C'].get_vector()
        rotate(n, ca, atom.get_vector(), angle)
    write('15477.pdb', structure)

if __name__ == "__main__":
    main()
