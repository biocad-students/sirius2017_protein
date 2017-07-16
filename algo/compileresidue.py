from Bio.PDB import *
import sys
sys.path.append("../")
from utils.io import *
from sampling.smartsampling import *

path = "../pdb/3-stor/"

def getresidue(name):
    return read(path+name+".pdb")

result = cleversamp('ADKAS')
writeres("out.pdb",result)
