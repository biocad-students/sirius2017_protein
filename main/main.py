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
from sampling.sampl_1 import samples
import timeit

IssmartWork = True
IsDebugReq = False
COUNT = 10

def debugI(name,structure):
    if(IsDebugReq):
        file = open("debug.log","a")
        file.write(name)
        try:
            for x in structure:
                file.write(x)
        except:
                file.write(str(structure))
        file.write("\n\n")

def Work():
    print("Booting settings:\nIssmartWork {}\n IsDebugReq {}\n count {}".format(IssmartWork,IsDebugReq,COUNT))
    regionPath = "../sirius_out/regions.txt"
    structsPath = "../sirius_out/"

    folderwithresult = "result/"
    os.mkdir(folderwithresult)
    cdr3 = []
    file1 = open(regionPath, "r")
    _tmp = file1.read()
    file1.close()
    lines = _tmp.split('\n')
    tmp = []
    for x in lines:
        tmp.append(x.split(' '))

    for x in range(len(tmp)):
        cordstart = len(tmp[x][1]+tmp[x][2]+tmp[x][3]+tmp[x][4]+tmp[x][5])
        cordstop = cordstart+len(tmp[x][6])
        cdr3.append([tmp[x][0],cordstart,cordstop])
## MAIN PROGRAM ##
    for counter in range(9,len(tmp)):
        ourpdb = read(structsPath+str(cdr3[counter][0])+".pdb")
        ourres = ourpdb.child_list
        ourchain = Chain(0)
        firstRes = ourres[cdr3[counter][1]-1]
        lastRes =  ourres[cdr3[counter][2]]
        for x in range(cdr3[counter][1]-1,cdr3[counter][2]+1):
                ourchain.add(ourres[x])
        # 1 STEP
        letterList = [letter(x) for x in ourchain.child_list]
        os.mkdir(folderwithresult+str(cdr3[counter][0]))
        if(IssmartWork):
            _preSub = loopSubstring(''.join(letterList),COUNT,cdr3[counter][0])
            for fileenum in range(len(_preSub)):
                directory = str(cdr3[counter][0])+"/"
                # Собирает цепочки по буквам
                sub = get_residues_by_pos(_preSub[fileenum])
                for i in range(len(sub)):
                    debugI(str(i), sub[i])
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
                writeres(folderwithresult+directory+str(fileenum)+".pdb",chainArray)
        else:
                sa = samples(letterList,2)
                directory = str(cdr3[counter][0])+"/"
                for instance in range(len(sa)):
                    index = 0
                    print(type(sa),type(sa[0]))

                    combined = imposer(sa[instance],firstRes,lastRes)
                    debugI("combined",combined)
                    #5 CCD
                    afterCCD = CCD(combined,lastRes,True)
                    # 6 скл
                    firstPart = ourres[0:cdr3[counter][1]]
                    secondPart = ourres[cdr3[counter][2]:]
                    chainArray = firstPart+afterCCD[1:-1]+secondPart
                    writeres(folderwithresult+directory+str(instance)+".pdb",chainArray)
# def timetest():
#     time = timeit.Timer(setup=smartWork).repeat(10)
#     summ = 0
#     for x in time:
#         summ+=x
#     print(summ/len(time))
# timetest()
Work()

# ОПТИМИЗАЦИИ
# 0.015892415899725166 - 1 результат
# 0.015419271100472542 - убрал pow
# 0.016497708200040505 - убрали библиотеки
# 0.01711203030008619 - красиво сделал
# 0.015177678400141304 - изменил все функции на * в модулях
