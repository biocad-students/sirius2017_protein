
import sys
sys.path.append("../")

from utils.io import *
from geometry.transform import *
from geometry.movements import *
from utils.final_check import final_check
from math import pi

# SHIFT
# struct = read("../pdb/15477.pdb")
# struct = shift(struct,1.8)
# write("examples/ShiftExample.pdb",struct)

# ROTATION
# struct = read("../pdb/15477.pdb")
# struct = rot(struct,pi)
# write("examples/rotExample.pdb",struct)

# MERGE
loop = read("../pdb/correct/cdr.pdb")
protein = read("../pdb/correct/frs 22.pdb")
loop = imposer(loop,protein[31],protein[32])
writeres("examples/LoopAndFRSExample.pdb",loop)


#final_check("../../../Desktop/sirius_out/regions.txt","../../../Desktop/our/","../../../Desktop/biocad/")
