from Bio.PDB import *
from Bio.PDB.Structure import Structure
from utils.utilits import *
from geometry.transform import *
from sampling.sampl_1 import rot1
import numpy,math
from utils.io import *


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


def cleversamp(var):
    """
        3 вариант сэмплирования
        Параметры:
            structure - список списков 3-меров
    """
    finalresidue = getresidue(var[0:3])[0]
    resindex = 0
    index = 1
    while(len(var)+1>index+3):
        rms = 1000
        isnew = True
        tmpresidue = finalresidue
        for residues in getresidue(var[index:index+3]):
            finalresidue = tmpresidue
            fixed = [finalresidue[-2]['N'],finalresidue[-2]['CA'],finalresidue[-2]['C'],finalresidue[-1]['N'],finalresidue[-1]['CA'],finalresidue[-1]['C']]
            finalresidue.pop()
            finalresidue.pop()
            print(residues)
            finalresidue.extend(residues)
            moving = [finalresidue[-3]['N'],finalresidue[-3]['CA'],finalresidue[-3]['C'],finalresidue[-2]['N'],finalresidue[-2]['CA'],finalresidue[-2]['C']]
            sp = Superimposer()
            sp.set_atoms(fixed,moving)
            for residue in range(index,len(finalresidue)):
                for atom in finalresidue[residue]:
                    v = atom.get_vector()
                    cord = numpy.dot(v._ar,sp.rotran[0])+sp.rotran[1]
                    atom.set_coord(cord)
            index+=1
            if(isnew):
                best = finalresidue
            else:
                trms = rmsd(best,finalresidue)
                if(trms<rms):
                    best = finalresidue
    for x in best:
        x.id = (' ',best,' ')
        best+=1
    return best
