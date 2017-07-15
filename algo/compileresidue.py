from Bio.PDB import *
import sys
sys.path.append("../")
from utils.io import *
from sampling.smartsampling import *

path = "../pdb/3-stor/"

def getresidue(name):
    return read(path+name+".pdb")


a = getresidue('ASG')
b = getresidue('SGA')
c = getresidue('AYL')
result = cleversamp([[a],[b],[c]])
writeres("out.pdb",result)
