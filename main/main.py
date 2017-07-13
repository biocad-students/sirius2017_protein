from Bio.PDB import *
import sys
sys.path.append("../")
from utils.io import *
from utils.calc import *
from utils.strstr import *
from utils.utilits import *
from utils.withdraw import *
from algo.CCD import CCD
from sampling.smartsampling import *

from Bio.PDB.Chain import Chain
regionPath = "../sirius_out/regions.txt"
structsPath = "../sirius_out/"

cdr3 = []
file1 = open(regionPath, "r")
_tmp = file1.read()
file1.close()
lines = _tmp.split('\n')
tmp = []
for x in lines:
    tmp.append(x.split(' '))
for x in range(0,10):
    cordstart = len(tmp[x][1]+tmp[x][2]+tmp[x][3]+tmp[x][4]+tmp[x][5])
    cordstop = cordstart+len(tmp[x][6])
    cdr3.append([tmp[x][0],cordstart,cordstop])

#for counter in range(0,len(cdr3)):
for counter in range(0,2):
    ourpdb = read(structsPath+str(cdr3[counter][0])+".pdb")
    ourres = ourpdb.child_list

    ourchain = Chain(0)
    firstRes = ourres[cdr3[counter][1]-1]
    lastRes =  ourres[cdr3[counter][2]+1]
    print('cords',firstRes,lastRes)
    for x in range(cdr3[counter][1],cdr3[counter][2]):
            ourchain.add(ourres[x])

    # 1 STEP
    letterList = [letter(x) for x in ourchain.child_list]
    preSub = loopSubstring(''.join(letterList),2,cdr3[counter][0])
#    print('tutut ',preSub,cdr3[counter][0])
    # 2 STEP
    sub = get_residues_by_pos(preSub[0])
    # 3 STEP
    merged = smartsamp(sub)
    #4 STEP
    combined = imposer(merged,firstRes,lastRes)
    writeres("beforeccd.pdb",combined)
    write("last.pdb",lastRes.copy())
    #5 STEP
    #print(combined)
    afterCCD = CCD(combined,lastRes)
    # 6 STEP
    firstPart = ourres[0:cdr3[counter][1]]
    secondPart = ourres[cdr3[counter][2]:]
    chainArray = firstPart+afterCCD+secondPart
#    print("POS:",cdr3[counter][1],cdr3[counter][2])
    writeres(str(counter)+".pdb",chainArray)
