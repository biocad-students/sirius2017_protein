from Bio.PDB import *
import os,sys,platform
import multiprocessing
sys.path.append("../")

from utils.io import *
from utils.calc import *
from utils.strstr import *
from utils.utilits import *
from utils.withdraw import *
from sampling.smartsampling import *
from Bio.PDB.Chain import Chain
from algo.CCD import *
from sampling.sampl_1 import samples

IssmartWork = True
IsDebugReq = False
THREADNUM = 1
COUNT = 10
regionPath = "../../../Desktop/sirius_out/regions.txt"
structsPath = "../../../Desktop/sirius_out/"
folderwithresult = "result/"


def debugI(name,structure):
	if(IsDebugReq):
		print(name)
		try:
			for x in structure:
				print(x)
		except:
				print(str(structure))

def preparing():
    print("Booting settings:\nIssmartWork {}\n IsDebugReq {}\n count {}".format(IssmartWork,IsDebugReq,COUNT))
    os.system("mkdir -p "+folderwithresult)
    try:
        file = open((folderwithresult+"info.log"),"w")
        file.write("Sys info:")
        file.write(str(platform.processor())+"\n")
        file.write(str(platform.machine())+"#"+"\n")
        file.write(str(platform.uname()))
        file.write("\nWe running:"+str(THREADNUM)+"\nNumber of files:"+str(COUNT))
        file.close()
    except:
        pass
    cdr3 = []
    with open(regionPath, "r") as f:
        string = f.readline().split()
        while string:
            cordstart = len(string[1]+string[2]+string[3]+string[4]+string[5])
            cordstop = cordstart+len(string[6])
            cdr3.append([string[0],cordstart,cordstop])
            string = f.readline().split()
    threads = []
    lengthP = (len(cdr3)+1)/THREADNUM
    calcpos = 0

    for n in range(THREADNUM):
        print("new thread: #",n," with range [",round(calcpos),",",round(calcpos+lengthP),"]\n")
        threads.append(multiprocessing.Process(target=Work,args = (cdr3,round(calcpos),round(calcpos+lengthP))))
        calcpos+=lengthP

    for thrd in threads:
        thrd.start()
    for thrd in threads:
        thrd.join()

def Work(cdr3,calcstart,calcstop):
    print("Booting thread #",os.getpid())
## MAIN PROGRAM ##
    for counter in range(calcstart,calcstop):
        ourpdb = read(structsPath+str(cdr3[counter][0])+".pdb")
        ourres = ourpdb.child_list
        ourchain = Chain(0)
        firstRes = ourres[cdr3[counter][1]-1]
        lastRes =  ourres[cdr3[counter][2]]
        for x in range(cdr3[counter][1]-1,cdr3[counter][2]+1):
            ourchain.add(ourres[x])
        print("Working with ",cdr3[counter])
        # 1 STEP
        letterList = [letter(x) for x in ourchain.child_list]
        os.mkdir(folderwithresult+str(cdr3[counter][0]))
        if(IssmartWork):
            _preSub = loopSubstring(''.join(letterList),COUNT,cdr3[counter][0],structsPath)
            for fileenum in range(len(_preSub)):
                directory = str(cdr3[counter][0])+"/"
                # Собирает цепочки по буквам
                sub = get_residues_by_pos(_preSub[fileenum],prefix = structsPath)
                for i in range(len(sub)):
                    debugI(str(i), sub[i])
                # Соединяет цепокич в одну
                merged = smartsamp(sub)
                debugI("merged",merged)
                combined = imposer(merged,firstRes,lastRes)
                debugI("combined",combined)
                #5 CCD
                afterCCD = CCD(combined,lastRes,feedback = False)
                # 6 скл
                firstPart = ourres[0:cdr3[counter][1]]
                secondPart = ourres[cdr3[counter][2]:]
                chainArray = firstPart+afterCCD[1:-1]+secondPart
                writeres(folderwithresult+directory+str(fileenum)+".pdb",chainArray)
        else:
            sa = samples(letterList,COUNT)
            directory = str(cdr3[counter][0])+"/"
            for instance in range(len(sa)):
                index = 0
                combined = imposer(sa[instance],firstRes,lastRes)
                debugI("combined",combined)
                #5 CCD
                afterCCD = CCD(combined,lastRes,feedback = False)
                # 6 скл
                firstPart = ourres[0:cdr3[counter][1]]
                secondPart = ourres[cdr3[counter][2]:]
                chainArray = firstPart+afterCCD[1:-1]+secondPart
                writeres(folderwithresult+directory+str(instance)+".pdb",chainArray)
preparing()
