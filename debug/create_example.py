import sys
sys.path.append("../")

from utils.io import *
from geometry.transform import *
from geometry.movements import *
from utils.final_check import final_check
from math import pi
from sampling.smartsampling import *

# SHIFT
# struct = read("../pdb/15477.pdb")
# struct = shift(struct,1.8)
# write("examples/ShiftExample.pdb",struct)

# ROTATION
# struct = read("../pdb/15477.pdb")
# struct = rot(struct,pi)
# write("examples/rotExample.pdb",struct)

# MERGE
# loop = read("../pdb/correct/cdr.pdb")
# protein = read("../pdb/correct/frs 22.pdb")
# loop = imposer(loop.child_list,protein[31],protein[32])
# writeres("examples/LoopAndFRSExample.pdb",loop)

# FINAL TEST
#final_check("../../../Desktop/test_rmsd/regions.txt","../../../Desktop/test_rmsd/source/","../../../Desktop/test_rmsd/result/","../../../Desktop/result.txt")
#final_check("../../../Desktop/regions.txt","../sirius_out/","../../../Desktop/out/","../../../Desktop/result.txt")

#SMART SAMPLING
# print("ssadads")
# a = read("../pdb/1.pdb")
# b = read("../pdb/2.pdb")
# c = read("../pdb/3.pdb")
# struct = smartsamp([c.child_list,b.child_list,a.child_list])
# writeres("out.pdb",struct)

# CLEVER SAMPLING
# finalresidue = cleversamp('ADKASG')
# writeres("final.pdb",finalresidue)

# GRAPHIC

struct = read("molec.pdb")
index = 0
for x in struct:
    x.id = (' ',index,' ')
    index+=1
for x in range(5):
    struct = rot(struct,pi/6,start = 1,isNewCalcReq = 0)
    writeres("rot"+str(x)+".pdb",struct)
