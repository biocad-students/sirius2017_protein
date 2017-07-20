from Bio.PDB import *

from Bio.PDB.Structure import Structure
from utils.utilits import *
from geometry.transform import *
from sampling.sampl_1 import rot1
import numpy,math,sys
from utils.io import *
from sampling.read_k import *
from random import shuffle


def fastimpose(fixedRes,movingRes):
    fixed = [fixedRes[-1]['N'],fixedRes[-1]['CA'],fixedRes[-1]['C']]
    moving = [movingRes[0]['N'],movingRes[0]['CA'],movingRes[0]['C']]
    sp = Superimposer()
    sp.set_atoms(fixed,moving)
    for residue in movingRes:
        for atom in residue:
            v = atom.get_vector()
            cord = numpy.dot(v._ar,sp.rotran[0])+sp.rotran[1]
            atom.set_coord(cord)
    return movingRes

def smartsamp(structure):
    """
        Умный сэмплер
        Параметры:
            structure - список списков residue
    """
    struct = []
    index = 0
    for residueList in structure:
        if(residueList == structure[0]):
            for residue in residueList:
                _residue = residue.copy()
                _residue.id = (' ',index,' ')
                index+=1
                struct.append(_residue.copy())
        else:
            for residue in residueList:
                residue.detach_parent()
                residue.id = (' ',index,' ')
                index+=1
            toDel = len(struct)-1
            struct += fastimpose(struct,residueList)
            struct.pop(toDel)

    for x in range(1,len(struct)):
        if(letter(struct[x]) != 'P'):
            O = struct[x-1]['O'].get_vector()
            C = struct[x-1]['C'].get_vector()
            N = struct[x]['N'].get_vector()
            H = struct[x]['H'].get_vector()
            angle = calc_dihedral(O,C,N,H)
            if(abs(angle)<math.radians(160)):
                N = struct[x]['N'].get_vector()
                H = struct[x]['H'].get_vector()
                Vnc = O-N
                Vnh = C-N
                Dnh = distance(N,H)
                vertical = numpy.cross(Vnc.get_array(),Vnh.get_array())
                Ci = calcnewcord(N,N+vertical,C,math.radians(120))
                Hnew = normalize(Ci-N.get_array())*Dnh + N.get_array()
                struct[x]['H'].set_coord(Hnew)
    return struct

path = "../pdb/3-stor/"

def getresidue(name):
    if(name == "ASG"):
        return [read(path+name+"1.pdb").child_list,read(path+name+"2.pdb").child_list,read(path+name+".pdb").child_list]
    return [read(path+name+".pdb").child_list]



def cleversamp(var,mas,COUNT):
    """
        3 вариант сэмплирования
        Параметры:
            structure - список списков 3-меров
    """
    firstChar = var[0:3]
    firstList = mas.getValue(firstChar)
    shuffle(firstList)
    firstList = firstList[0:COUNT]
    for firstElement in firstList:
        todoLine = var[3:]
        tmp = generate_kmer(var,firstElement)
        while(len(todoLine)>0):
            angles = mas.getValue(todoLine[0:3])
            residue = []
            for angle in angles:
                print(todoLine[0:3])
                kmer = generate_kmer(todoLine[0:3],angle)

    exit()
