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

def bres(an,structure):
    print(an)
    try:
        for x in structure:
            print(x)
            print(x.child_list)
    except:
        print("что-то пошло не так")
    print("")
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
    print(tmp[x][0],cordstart,cordstop, tmp[x])

#for counter in range(0,len(cdr3)):
for counter in range(1):
    print("File: ", cdr3[counter][0])
    ourpdb = read(structsPath+str(cdr3[counter][0])+".pdb")
    ourres = ourpdb.child_list

    ourchain = Chain(0)
    firstRes = ourres[cdr3[counter][1]-1]
    lastRes =  ourres[cdr3[counter][2]]
    for x in range(cdr3[counter][1]-1,cdr3[counter][2]+1):
            ourchain.add(ourres[x])

    # 1 STEP
    letterList = [letter(x) for x in ourchain.child_list]
    print("CDR3:", letterList)
    preSub = loopSubstring(''.join(letterList),2,cdr3[counter][0])
    # Собирает цепочки по буквам
    sub = get_residues_by_pos(preSub[0])
    for i in range(len(sub)):
        bres(str(i), sub[i])
    # Соединяет цепокич в одну
    merged = smartsamp(sub)
    bres("merged", merged)
    print("firstRes", firstRes)
    print("lastRes", lastRes)
    combined = imposer(merged,firstRes,lastRes)
    bres("combined", combined)
    writeres("combined.pdb",combined)
    writeres("combined_target.pdb",[lastRes])
    #5 CCD
    afterCCD = CCD(combined,lastRes,True)
    # 6 скл
    firstPart = ourres[0:cdr3[counter][1]]
    secondPart = ourres[cdr3[counter][2]:]
    chainArray = firstPart+afterCCD[1:-1]+secondPart
    writeres(str(counter)+".pdb",chainArray)
