from Bio.PDB import *
import sys
import os
sys.path.append("../")
# from utils.io import read,writeres
# from utils.calc import imposer
# from utils.strstr import loopSubstring
# from utils.utilits import letter
# from utils.withdraw import get_residues_by_pos
# from sampling.smartsampling import smartsamp
# from Bio.PDB.Chain import Chain
# from algo.CCD import CCD
# import timeit

from utils.io import *
from utils.calc import *
from utils.strstr import *
from utils.utilits import *
from utils.withdraw import *
from sampling.smartsampling import *
from Bio.PDB.Chain import Chain
from algo.CCD import *
import timeit


def debugI(name,structure):
    print(name)
    for x in structure:
        print(x)
    print("\n\n")

def smartWork():
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
    for counter in range(1):
        ourpdb = read(structsPath+str(cdr3[counter][0])+".pdb")
        ourres = ourpdb.child_list
        ourchain = Chain(0)
        firstRes = ourres[cdr3[counter][1]-1]
        lastRes =  ourres[cdr3[counter][2]]
        for x in range(cdr3[counter][1]-1,cdr3[counter][2]+1):
                ourchain.add(ourres[x])
        print("chain",ourchain.child_list)
        # 1 STEP
        letterList = [letter(x) for x in ourchain.child_list]
        preSub = loopSubstring(''.join(letterList),2,cdr3[counter][0])
        os.mkdir(str(cdr3[counter][0]))
        directory = str(cdr3[counter][0])+"/"
        # Собирает цепочки по буквам
        sub = get_residues_by_pos(preSub[0])
        debugI("sub",sub)
        # Соединяет цепокич в одну
        merged = smartsamp(sub)
        debugI("merged",merged)
        combined = imposer(merged,firstRes,lastRes)
        debugI("combined",combined)
        #5 CCD
        afterCCD = CCD(combined,lastRes)
        # 6 скл
        firstPart = ourres[0:cdr3[counter][1]]
        secondPart = ourres[cdr3[counter][2]:]
        chainArray = firstPart+afterCCD[1:-1]+secondPart
        writeres(directory+str(counter)+".pdb",chainArray)



def defaultWork():
    # Данные


    # Паша


    #Люда

    pass
# def timetest():
#     time = timeit.Timer(setup=smartWork).repeat(10)
#     summ = 0
#     for x in time:
#         summ+=x
#     print(summ/len(time))
# timetest()
smartWork()

# ОПТИМИЗАЦИИ
# 0.015892415899725166 - 1 результат
# 0.015419271100472542 - убрал pow
# 0.016497708200040505 - убрали библиотеки
# 0.01711203030008619 - красиво сделал
# 0.015177678400141304 - изменил все функции на * в модулях
