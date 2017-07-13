from Bio.PDB import *
from Bio.PDB.Structure import Structure
from utils.calc import link
from utils.utilits import letter
from geometry.transform import *
import numpy



def smartimposer(ta,structure,splice):
    for x in splice:
        #print(x)
        fixed = [structure[x-1]['N'],structure[x-1]['CA'],structure[x-1]['C']]
        moving = [structure[x]['N'],structure[x]['CA'],structure[x]['C']]
        sp = Superimposer()
        sp.set_atoms(fixed,moving)

        for y in range(x,len(structure)):
            for r in structure[y]:
                v = r.get_vector()
                cord = numpy.dot(v._ar,sp.rotran[0])+sp.rotran[1]
                structure[y][r.get_name()].set_coord(cord)
        # for x in range(1,len(structure)):
        #     if((letter(structure[x-1]) != 'P') and (letter(structure[x]) != 'P')):
        #         O = structure[x-1]['O'].get_vector()
        #         C = structure[x-1]['C'].get_vector()
        #         N = structure[x]['N'].get_vector()
        #         H = structure[x]['H'].get_vector()
        #         angle = calc_dihedral(O,C,N,H)
        #         if(abs(angle)<math.radians(100)):
        #             print("turn")
        #             structure = rot(structure,math.pi-angle+0.3,x,0)
    for x in splice:
        structure[x-1] = 0
    for x in splice:
        structure.remove(0)
    structure.pop()

    for x in range(1,len(structure)):
        if(letter(structure[x]) != 'P'):
            O = structure[x-1]['O'].get_vector()
            C = structure[x-1]['C'].get_vector()
            N = structure[x]['N'].get_vector()
            H = structure[x]['H'].get_vector()
            angle = calc_dihedral(O,C,N,H)
            if(abs(angle)<math.radians(160)):
                structure = rot(structure,math.pi-angle,x,0)
                print("turn on angle ",math.degrees(angle))
    debugI(structure)
    return structure

def debugI(struct):
    for x in struct:
        print(x.id[1],x.resname)

def smartsamp(structure):
    """
        Умный сэмплер
        Параметры:
            structure - список списков residue
    """
    struct = list()
    splice = list()
    indexcount = 0 # индекс структуры
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
