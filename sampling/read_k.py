from Bio.SeqUtils import *
from geometry.transform import *
from utils.io import *
from sampl_1 import *

def copy_to(amino1, amino2, arr_2=[]):
    arr=['H', 'N', 'CA', 'C', 'O', 'HA']+arr_2
    for atom_name in arr:
       if(amino2.has_id(atom_name)):
            amino1[atom_name].set_coord(amino2[atom_name].get_vector())

def rotate_radical(amino1, y, z):
    """выставляет двугранный угол z=c-ca-n-cb и плоский угол n-ca-cb=y"""
    amino2=amino1.copy()
    amino3 = rot1(amino2, z,  amino2['C'].get_vector(), amino2['CA'].get_vector(), amino2['N'].get_vector(), amino2['CB'].get_vector())
    copy_to(amino3, amino2)
    amino3_copy=amino3.copy()
    amino = rotate_amino(amino3, amino3_copy['N'].get_vector(), amino3_copy['CA'].get_vector(), amino3_copy['CB'].get_vector(), amino3_copy['CB'].get_vector(), y)
    copy_to(amino, amino3_copy)
    return amino

def change_dist(amino, dis, A, B):
    coord = A - B
    for atom in amino:
        atom.set_coord(atom.get_vector() + coord)
    coord.normalize()
    coord[0] *= -dis
    coord[1] *= -dis
    coord[2] *= -dis
    for atom in amino:
        atom.set_coord(atom.get_vector() + coord)
    return amino

def rotate_main(_amino1, _amino2, main):
    amino1=_amino1.copy()
    amino2=_amino2.copy()
    rotated_amino1 = rotate_amino(amino1, amino1['N'].get_vector(), amino1['CA'].get_vector(), amino1['C'].get_vector(), amino1['C'].get_vector(), main[6])
    rotated_amino2 = rotate_amino(amino2, amino2['N'].get_vector(), amino2['CA'].get_vector(), amino2['C'].get_vector(), amino2['C'].get_vector(), main[7])
    rotated_amino2 = moveTo(rotated_amino1, rotated_amino2, [main[0], main[1], main[2]], main[5], main[4], main[3])
    return [rotated_amino1, rotated_amino2]

def rot_val(amino1, x1):
    amino = amino1.copy()
    amino_copy=amino.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(), amino['CG1'].get_vector())
    copy_to(rotated_amino, amino_copy)
    return rotated_amino

def rot_ile(amino1, x1, x2):
    amino = amino1.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(),amino['CG1'].get_vector())
    copy_to(rotated_amino, amino1.copy())
    rotated = rot1(rotated_amino, x2, rotated_amino['CA'].get_vector(), rotated_amino['CB'].get_vector(), rotated_amino['CG1'].get_vector(), rotated_amino['CD1'].get_vector())
    copy_to(rotated, rotated_amino.copy(), ['HB', 'CG2', 'HG21', 'HG22', 'HG23'])
    return rotated

def rot_leu(amino1, x1, x2):
    amino = amino1.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(),amino['CG'].get_vector())
    copy_to(rotated_amino, amino1.copy())
    rotated = rot1(rotated_amino, x2, rotated_amino['CA'].get_vector(), rotated_amino['CB'].get_vector(), rotated_amino['CG'].get_vector(), rotated_amino['CD1'].get_vector())
    copy_to(rotated, rotated_amino.copy(), ['HB1', 'HB2'])
    return rotated

def rot_met(amino1, x1, x2, x3):
    amino = amino1.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(), amino['CG'].get_vector())
    copy_to(rotated_amino, amino1.copy())
    rotated = rot1(rotated_amino, x2, rotated_amino['CA'].get_vector(), rotated_amino['CB'].get_vector(), rotated_amino['CG'].get_vector(), rotated_amino['SD'].get_vector())
    copy_to(rotated, rotated_amino.copy(), ['HB1', 'HB2'])
    rotated_2 = rot1(rotated, x3, rotated['CB'].get_vector(), rotated['CG'].get_vector(),rotated['SD'].get_vector(), rotated['CE'].get_vector())
    copy_to(rotated_2, rotated.copy(), ['HB1', 'HB2', 'HG1', 'HG2', 'CB'])
    return rotated_2

def rot_gln(amino1, x1, x2, x3):
    amino = amino1.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(), amino['CG'].get_vector())
    copy_to(rotated_amino, amino1.copy())
    rotated = rot1(rotated_amino, x2, rotated_amino['CA'].get_vector(), rotated_amino['CB'].get_vector(), rotated_amino['CG'].get_vector(), rotated_amino['CD'].get_vector())
    copy_to(rotated, rotated_amino.copy(), ['HB1', 'HB2'])
    rotated_2 = rot1(rotated, x3, rotated['CB'].get_vector(), rotated['CG'].get_vector(),rotated['CD'].get_vector(), rotated['OE1'].get_vector())
    copy_to(rotated_2, rotated.copy(), ['HB1', 'HB2', 'HG1', 'HG2', 'CB'])
    return rotated_2

def rot_glu(amino1, x1, x2, x3):
    amino = amino1.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(), amino['CG'].get_vector())
    copy_to(rotated_amino, amino1.copy())
    rotated = rot1(rotated_amino, x2, rotated_amino['CA'].get_vector(), rotated_amino['CB'].get_vector(), rotated_amino['CG'].get_vector(), rotated_amino['CD'].get_vector())
    copy_to(rotated, rotated_amino.copy(), ['HB1', 'HB2'])
    rotated_2 = rot1(rotated, x3, rotated['CB'].get_vector(), rotated['CG'].get_vector(),rotated['CD'].get_vector(), rotated['OE1'].get_vector())
    copy_to(rotated_2, rotated.copy(), ['HB1', 'HB2', 'HG1', 'HG2', 'CB'])
    return rotated_2

def rot_ser(amino1, x1):
    amino = amino1.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(),amino['OG'].get_vector())
    copy_to(rotated_amino, amino1.copy())
    return rotated_amino

def rot_thr(amino1, x1):
    amino = amino1.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(),amino['OG1'].get_vector())
    copy_to(rotated_amino, amino1.copy())
    return rotated_amino

def rot_cys(amino1, x1):
    amino = amino1.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(),amino['SG'].get_vector())
    copy_to(rotated_amino, amino1.copy())
    return rotated_amino

def rot_tyr_phe_trp(amino1, x1, x2):
    amino = amino1.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(),amino['CG'].get_vector())
    copy_to(rotated_amino, amino1.copy())
    rotated = rot1(rotated_amino, x2, rotated_amino['CA'].get_vector(), rotated_amino['CB'].get_vector(), rotated_amino['CG'].get_vector(), rotated_amino['CD1'].get_vector())
    copy_to(rotated, rotated_amino.copy(), ['HB1', 'HB2'])
    return rotated

def rot_asn(amino1, x1, x2):
    amino = amino1.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(),amino['CG'].get_vector())
    copy_to(rotated_amino, amino1.copy())
    rotated = rot1(rotated_amino, x2, rotated_amino['CA'].get_vector(), rotated_amino['CB'].get_vector(), rotated_amino['CG'].get_vector(), rotated_amino['OD1'].get_vector())
    copy_to(rotated, rotated_amino.copy(), ['HB1', 'HB2'])
    return rotated

def rot_his(amino1, x1, x2):
    amino = amino1.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(),amino['CG'].get_vector())
    copy_to(rotated_amino, amino1.copy())
    rotated = rot1(rotated_amino, x2, rotated_amino['CA'].get_vector(), rotated_amino['CB'].get_vector(), rotated_amino['CG'].get_vector(), rotated_amino['ND1'].get_vector())
    copy_to(rotated, rotated_amino.copy(), ['HB1', 'HB2'])
    return rotated

def rot_asp(amino1, x1, x2):
    amino=amino1.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(),amino['CG'].get_vector())
    copy_to(rotated_amino, amino1.copy())
    rotated = rot1(rotated_amino, x2, rotated_amino['CA'].get_vector(), rotated_amino['CB'].get_vector(), rotated_amino['CG'].get_vector(), rotated_amino['OD1'].get_vector())
    copy_to(rotated, rotated_amino, ['HB1', 'HB2'])
    return rotated

def rot_pro(amino1, x1, x2):
    amino = amino1.copy()
    CDH=distance(amino['CD'].get_vector(), amino['HD1'].get_vector())
    HCH=culc_angle(amino['HD1'].get_vector(), amino['CD'].get_vector(), amino['HD2'].get_vector())
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(),amino['CG'].get_vector())
    copy_to(rotated_amino, amino1.copy())
    rotated = rot1(rotated_amino, x2, rotated_amino['CA'].get_vector(), rotated_amino['CB'].get_vector(), rotated_amino['CG'].get_vector(), rotated_amino['CD'].get_vector())
    copy_to(rotated, rotated_amino, ['HB1', 'HB2', 'HB3'])
    rotated_CD = rotate_vector(rotated['N'].get_vector().get_array(), rotated['CG'].get_vector().get_array(), rotated['CD'].get_vector().get_array(), rad(180))
    CD = rotated['CD'].get_vector()
    H_coord=CD+normalize(rotated['CD'].get_vector().get_array()-rotated_CD)*CDH
    point_in_plane=rotate_vector(rotated['CD'].get_vector().get_array(), rotated_CD, amino['N'].get_vector().get_array(), rad(90))
    HD1 = rotate_vector(rotated['CD'].get_vector().get_array(), point_in_plane, H_coord.get_array(), HCH/2)
    HD2 = rotate_vector(rotated['CD'].get_vector().get_array(), point_in_plane, H_coord.get_array(), -HCH/2)
    HD1 = rotate_vector(rotated['CD'].get_vector().get_array(), rotated_CD, HD1, rad(90))
    HD2 = rotate_vector(rotated['CD'].get_vector().get_array(), rotated_CD, HD2, rad(90))
    rotated['HD1'].set_coord(HD1)
    rotated['HD2'].set_coord(HD2)
    return rotated

def rot_arg(amino1, x1, x2, x3, x4, x5):
    amino = amino1.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(), amino['CG'].get_vector())
    copy_to(rotated_amino, amino1.copy())
    rotated = rot1(rotated_amino, x2, rotated_amino['CA'].get_vector(), rotated_amino['CB'].get_vector(), rotated_amino['CG'].get_vector(), rotated_amino['CD'].get_vector())
    copy_to(rotated, rotated_amino.copy(), ['HB1', 'HB2'])
    rotated_2 = rot1(rotated, x3, rotated['CB'].get_vector(), rotated['CG'].get_vector(),rotated['CD'].get_vector(), rotated['NE'].get_vector())
    copy_to(rotated_2, rotated.copy(), ['HB1', 'HB2', 'HG1', 'HG2', 'CB'])
    rotated_3 = rot1(rotated_2, x4, rotated_2['CG'].get_vector(), rotated_2['CD'].get_vector(), rotated_2['NE'].get_vector(),rotated_2['CZ'].get_vector())
    copy_to(rotated_3, rotated_2.copy(), ['HB1', 'HB2', 'HG1', 'HG2', 'CB', 'CG', 'HD1', 'HD2'])
    rotated_4 = rot1(rotated_3, x5, rotated_3['CD'].get_vector(), rotated_3['NE'].get_vector(), rotated_3['CZ'].get_vector(),rotated_3['NH1'].get_vector())
    copy_to(rotated_4, rotated_3.copy(), ['HB1', 'HB2', 'HG1', 'HG2', 'CB', 'CG', 'HD1', 'HD2', 'HE', 'CD'])
    return rotated_4

def rot_lys(amino1, x1, x2, x3, x4):
    amino = amino1.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(), amino['CG'].get_vector())
    copy_to(rotated_amino, amino1.copy())
    rotated = rot1(rotated_amino, x2, rotated_amino['CA'].get_vector(), rotated_amino['CB'].get_vector(), rotated_amino['CG'].get_vector(), rotated_amino['CD'].get_vector())
    copy_to(rotated, rotated_amino.copy(), ['HB1', 'HB2'])
    rotated_2 = rot1(rotated, x3, rotated['CB'].get_vector(), rotated['CG'].get_vector(),rotated['CD'].get_vector(), rotated['CE'].get_vector())
    copy_to(rotated_2, rotated.copy(), ['HB1', 'HB2', 'HG1', 'HG2', 'CB'])
    rotated_3 = rot1(rotated_2, x4, rotated_2['CG'].get_vector(), rotated_2['CD'].get_vector(), rotated_2['CE'].get_vector(),rotated_2['NZ'].get_vector())
    copy_to(rotated_3, rotated_2.copy(), ['HB1', 'HB2', 'HG1', 'HG2', 'CB', 'CG', 'HD1', 'HD2'])
    return rotated_3

def rotRadicalInAmino(amino1, angles, letter):
    if(letter!='A' and letter!='G'):
        amino2 = amino1.copy()
        #print(letter, ' ', angles)
        amino = rotate_radical(amino2, rad(angles[0]), rad(angles[1]))
    else:
        amino = amino1
    if (letter == 'D'):
        amino = rot_asp(amino, rad(angles[2]), rad(angles[3]))
    if (letter == 'E'):
        amino = rot_glu(amino, rad(angles[2]), rad(angles[3]), rad(angles[4]))
    if (letter == 'C'):
        amino = rot_cys(amino, rad(angles[2]))
    if (letter == 'H'):
        amino = rot_his(amino, rad(angles[2]), rad(angles[3]))
    if (letter == 'I'):
        amino = rot_ile(amino, rad(angles[2]), rad(angles[3]))
    if (letter == 'K'):
        amino = rot_lys(amino, rad(angles[2]), rad(angles[3]), rad(angles[4]), rad(angles[5]))
    if (letter == 'L'):
        amino = rot_leu(amino, rad(angles[2]), rad(angles[3]))
    if (letter == 'M'):
        amino = rot_met(amino, rad(angles[2]), rad(angles[3]), rad(angles[4]))
    if (letter == 'N'):
        amino = rot_asn(amino, rad(angles[2]), rad(angles[3]))
    if (letter == 'P'):
        amino = rot_pro(amino, rad(angles[2]), rad(angles[3]))
    if (letter == 'Q'):
        amino = rot_gln(amino, rad(angles[2]), rad(angles[3]), rad(angles[4]))
    if (letter == 'R'):
        amino = rot_arg(amino, rad(angles[2]), rad(angles[3]), rad(angles[4]), rad(angles[5]), rad(angles[6]))
    if (letter == 'S'):
        amino = rot_ser(amino, rad(angles[2]))
    if (letter == 'T'):
        amino = rot_thr(amino, rad(angles[2]))
    if (letter == 'V'):
        amino = rot_asp(amino, rad(angles[2]), rad(angles[3]))
    if (letter == 'Y' or letter=='W' or letter=='F'):
        amino = rot_tyr_phe_trp(amino, rad(angles[2]), rad(angles[3]))
    return amino

def cnt_angles(letter):
    one_angle = ['V', 'S', 'T', 'C']
    two_angles = ['I', 'L', 'F', 'Y', 'W', 'N', 'P', 'H', 'D']
    three_angles = ['M', 'Q', 'E']
    four_angles = ['K']
    five_angles = ['R']
    if(letter in one_angle):
        return 1
    elif(letter in two_angles):
        return 2
    elif(letter in three_angles):
        return 3
    elif(letter in four_angles):
        return 4
    elif(letter in five_angles):
        return 5
    else:
        return -2

def main():
    structure = read('/home/ludmila/git/sirius2017_protein/sampling/aminos_out.pdb', "test")
    arr = {}
    for residue in structure:
        arr[seq1(residue.get_resname())] = residue
    am1='T'
    am2='A'
    am3='S'
    angles_and_dist=[-3.7548156, 178.12007, -106.46927, 119.63612, 123.63427, 1.3465911, 113.51516, 110.74306, 137.67319, 170.19287, -120.67031, 118.63731, 124.336075, 1.3409338, 110.74306, 110.09849, 173.68527, 179.59863, -154.52032, 120.319786, 123.57929, 1.3424668, 110.09849, 112.402664, 160.56654, 169.04048, -129.12004, 118.35881, 125.06439, 1.3427035, 112.402664, 107.13795, 112.36464, 121.535805, 73.16762, 110.11092, 126.9722, 75.921165]
    main1 = angles_and_dist[0:8]
    main2 = angles_and_dist[8:16]
    main3 = angles_and_dist[16:24]
    main4 = angles_and_dist[24:32]
    i=cnt_angles(am1)
    j=cnt_angles(am2)
    k=cnt_angles(am3)
    radical1 = angles_and_dist[32:34+i]
    radical2 = angles_and_dist[34+i:36+i+j]
    radical3 = angles_and_dist[36+i+j:38+i+j+k]
    amino1 = arr[am1]
    amino2 = arr[am2]
    amino3 = arr[am3]
    amino1 = rotRadicalInAmino(amino1, radical1, am1)
    amino2 = rotRadicalInAmino(amino2, radical2, am2)
    amino3 = rotRadicalInAmino(amino3, radical3, am3)
    pair1 = rotate_main(arr['A'], amino1, main1)
    #writeres('/tmp/AT.pdb', pair1)
    pair2 = rotate_main(pair1[1], amino2, main2)
    #writeres('/tmp/TA.pdb', pair2)
    pair3 = rotate_main(pair2[1], amino3, main3)
    #writeres('/tmp/AS.pdb', pair3)
    pair4 = rotate_main(pair3[1], arr['A'], main4)
    #writeres('/tmp/SA.pdb', pair4)
    writeres('/tmp/TAS.pdb', [pair2[0]]+[pair3[0]]+[pair4[0]])

if __name__ == "__main__":
    main()