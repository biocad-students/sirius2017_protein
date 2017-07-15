from Bio.SeqUtils import *
from geometry.rotation import *
from geometry.transform import *
from utils.io import *
from sampl_1 import rad, rotate_amino, rot1, culc_angle

def copy_to(amino1, amino2, arr_2=[]):
    arr=['H', 'N', 'CA', 'C', 'O', 'HA']+arr_2
    for atom_name in arr:
       if(amino2.has_id(atom_name)):
            amino1[atom_name].set_coord(amino2[atom_name].get_vector())

def rot_val(amino1, x1):
    amino = amino1.copy()
    rotated_amino = rot1(amino, x1, amino['N'].get_vector(), amino['CA'].get_vector(), amino['CB'].get_vector(), amino['CG1'].get_vector())
    copy_to(rotated_amino, amino1.copy())
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

def rot_trh(amino1, x1):
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
    amino = amino1.copy()
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

def main():
    structure = read('/home/ludmila/git/sirius2017_protein/sampling/aminos_out.pdb', "test")
    arr = {}
    for residue in structure:
        arr[seq1(residue.get_resname())] = residue
    amino = arr['P']
    writeres('/tmp/P1.pdb', [amino])
    writeres('/tmp/P2.pdb', [rot_pro(amino, rad(27), rad(-36))])


if __name__ == "__main__":
    main()