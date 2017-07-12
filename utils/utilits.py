from Bio.PDB import *
from Bio.SeqUtils import *
def zerolim(value):
    if(value < 5e-14):
        return 0
    else:
        return value

def letter(residue):
    """
        Идентификатор аминокислоты по 3 буквам
        Параметры:
            residue - аминокислота
    """
    return seq1(residue.get_resname())
