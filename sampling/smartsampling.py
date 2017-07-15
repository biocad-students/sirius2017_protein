from Bio.PDB import *
from Bio.PDB.Structure import Structure
from utils.utilits import *
from geometry.transform import *
from sampling.sampl_1 import rot1
import numpy,math
from utils.io import *

def smartimposer(ta,structure,splice):
    for x in splice:
        fixed = [structure[x-1]['N'],structure[x-1]['CA'],structure[x-1]['C']]
        moving = [structure[x]['N'],structure[x]['CA'],structure[x]['C']]
        sp = Superimposer()
        sp.set_atoms(fixed,moving)

        for y in range(x,len(structure)):
            for r in structure[y]:
                v = r.get_vector()
                cord = numpy.dot(v._ar,sp.rotran[0])+sp.rotran[1]
                structure[y][r.get_name()].set_coord(cord)
            writeres("RRR"+str(splice)+".pdb",structure)
    for x in splice:
        structure[x-1] = 0
    for x in splice:
        structure.remove(0)
    structure.pop()
    index = 0
    _structure = []
    for x in structure:
        tmp = x.copy()
        tmp.id = (' ',index,' ')
        index+=1
        _structure.append(tmp)
    structure = _structure
    for x in range(1,len(structure)):
        if(letter(structure[x]) != 'P'):
            O = structure[x-1]['O'].get_vector()
            C = structure[x-1]['C'].get_vector()
            N = structure[x]['N'].get_vector()
            H = structure[x]['H'].get_vector()
            angle = calc_dihedral(O,C,N,H)
            if(abs(angle)<math.radians(160)):
                N = structure[x]['N'].get_vector()
                H = structure[x]['H'].get_vector()
                Vnc = O-N
                Vnh = C-N
                Dnh = distance(N,H)
                vertical = numpy.cross(Vnc.get_array(),Vnh.get_array())
                Ci = calcnewcord(N,N+vertical,C,math.radians(120))
                Hnew = normalize(Ci-N.get_array())*Dnh + N.get_array()
                structure[x]['H'].set_coord(Hnew)

    return structure


def smartsamp(structure):
    """
        Умный сэмплер
        Параметры:
            structure - список списков residue
    """
    struct = list()
    splice = list()
    indexcount = 0
    for counter in range(len(structure)):
        for x in structure[counter]:
            x.detach_parent()
            x.id = (' ',indexcount,' ')
            indexcount+=1
            struct.append(x)
        splice.append(indexcount)
    splice.pop()

    struct = smartimposer(structure,struct,splice)
    return struct

path = "../pdb/3-stor/"

def getresidue(name):
    return read(path+name+".pdb")

def cleversamp(var):
    """
        3 вариант сэмплирования
        Параметры:
            structure - список списков 3-меров
    """
    finalresidue = getresidue(var[0:3])
    resindex = 0
    for x in finalresidue:
        x.id = (' ',resindex,' ')
        resindex+=1
    finalstr = var[0:3]
    varlen = len(var)
    index = 1
    while(varlen>index):
        for x in finalresidue:
            print(x)
        print("\n\n")
        req = finalstr[index]+finalstr[index+1]+var[index+2]
        index+=1
        print("looking for:",req)
        res = getresidue(req)
        for residue in res:
            residue.id = (' ',resindex,' ')
            resindex+=1
            finalresidue.add(residue)
        # fixed = [finalresidue[index-1]['N'],finalresidue[index-1]['CA'],finalresidue[index-1]['C']]
        # moving = [res[index-2]['N'],res[index-2]['CA'],res[index-2]['C']]
        # sp = Superimposer()
        # sp.set_atoms(fixed,moving)
        # for residue in res:
        #     for atom in residue:
        #         v = atom.get_vector()
        #         cord = numpy.dot(v._ar,sp.rotran[0])+sp.rotran[1]
        #         atom.set_coord(cord)

    return finalresidue
