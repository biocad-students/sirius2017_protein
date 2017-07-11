
import sys
sys.path.append("../")

from utils.io import *
from geometry.transform import *
from geometry.movements import *
from math import pi
#
# struct = read("../pdb/15477.pdb")
# struct = shift(struct,1.8)
# write("examples/ShiftExample.pdb",struct)
#
# struct = read("../pdb/15477.pdb")
# struct = rot(struct,pi)
# write("examples/rotExample.pdb",struct)

loop = read("../pdb/correct/cdr.pdb")
protein = read("../pdb/correct/frs 22.pdb")
loop = imposer(loop,protein[31],protein[32])
# write("out.pdb",loop)

chain = compileres(loop)
