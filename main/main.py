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
from params.Class import MyClass


TypeOfWork = 3
# 1 - 1 ветка
# 2 - 2 ветка
# 3 - 3 ветка
IsDebugReq = False
THREADNUM = 1
COUNT = 1
if len(sys.argv) > 1:
    regionPath = sys.argv[1] + "regions.txt"
    structsPath = sys.argv[1]
else:
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
    print("Booting settings:\n TypeOfWork {}\n IsDebugReq {}\n count {}".format(TypeOfWork,IsDebugReq,COUNT))
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
            loop = string[6]
            cdr3.append([string[0],cordstart,cordstop,loop])
            string = f.readline().split()
    threads = []
    lengthP = (len(cdr3))/THREADNUM
    calcpos = 0
    files = os.listdir(folderwithresult)
    mas = MyClass()
    for n in range(THREADNUM):
        print("new thread: #",n," with range [",round(calcpos),",",round(calcpos+lengthP),"]\n")
        threads.append(multiprocessing.Process(target=Work,args = (cdr3,round(calcpos),round(calcpos+lengthP),files,mas)))
        calcpos+=lengthP

    for thrd in threads:
        thrd.start()
    for thrd in threads:
        thrd.join()

def Work(cdr3,calcstart,calcstop,files,mas):
    #calcstart = 321
    #calcstop = 323
    print("Booting thread #",os.getpid())
## MAIN PROGRAM ##
    for counter in range(calcstart,calcstop):
        iscontinue = False
        for x in files:
            if(x == cdr3[counter][0]):
                print("Skipped ",cdr3[counter])
                iscontinue = True
        if(iscontinue):
            continue
        ourpdb = read(structsPath+str(cdr3[counter][0])+".pdb")
        ourres = ourpdb.child_list
        ourchain = Chain(0)
        firstRes = ourres[cdr3[counter][1]-1]
        #firstRes = ourres[cdr3[counter][1]]
        lastRes =  ourres[cdr3[counter][2]]
        for x in range(cdr3[counter][1]-1,cdr3[counter][2]+1):
            ourchain.add(ourres[x])
        print("Working with ",cdr3[counter])
        # 1 STEP
        letterList = [letter(x) for x in ourchain.child_list]
        os.system("mkdir -p "+folderwithresult+str(cdr3[counter][0]))
        directory = str(cdr3[counter][0])+"/"
        if(TypeOfWork == 2):
            _preSub = loopSubstring(''.join(letterList),COUNT,None,structsPath)
            for fileenum in range(len(_preSub)):
                # Собирает цепочки по буквам
                sub = get_residues_by_pos(_preSub[fileenum],"../files/",structsPath)
                for i in range(len(sub)):
                    debugI(str(i), sub[i])
                # Соединяет цепокич в одну
                merged = smartsamp(sub)
                debugI("merged",merged)
                combined = imposer(merged,firstRes,lastRes)
                debugI("combined",combined)
                #5 CCD
                #writeres("")
                afterCCD = CCD(combined,lastRes,feedback = False)
                # 6 скл
                if(afterCCD == None):
                    print("there is no peteylka :(. ",counter,instance)
                else:
                    firstPart = ourres[0:cdr3[counter][1]-1]
                    secondPart = ourres[cdr3[counter][2]:]
                    chainArray = firstPart+afterCCD[1:-1]+secondPart
                    writeres(folderwithresult+directory+str(fileenum)+".pdb",chainArray)
        elif(TypeOfWork == 1):
            sa = samples(letterList,COUNT)
            for instance in range(len(sa)):
                index = 0
                combined = imposer(sa[instance],firstRes,lastRes)
                debugI("combined",combined)
                #5 CCD
                afterCCD = CCD(combined,lastRes,feedback = False)
                # 6 скл
                if(afterCCD == None):
                    print("there is no peteylka :(. ",counter,instance)
                else:
                    firstPart = ourres[0:cdr3[counter][1]]
                    secondPart = ourres[cdr3[counter][2]:]
                    chainArray = firstPart+afterCCD[1:-1]+secondPart
                    writeres(folderwithresult+directory+str(instance)+".pdb",chainArray)
        elif(TypeOfWork == 3):
                listoflistofres = cleversamp(''.join(letterList),mas,COUNT)
                for iterator in range(len(listoflistofres)):
                    combined = imposer(listoflistofres[iterator],firstRes,lastRes)
                    debugI("combined",combined)
                    #5 CCD
                    afterCCD = CCD(combined,lastRes,feedback = False)
                    # 6 скл
                    if(afterCCD == None):
                        print("there is no peteylka :(. ",counter,iterator)
                    else:
                        firstPart = ourres[0:cdr3[counter][1]]
                        secondPart = ourres[cdr3[counter][2]:]
                        chainArray = firstPart+afterCCD[1:-1]+secondPart
                        writeres(folderwithresult+directory+str(iterator)+".pdb",chainArray)
        else:
            print("Unknown way ",TypeOfWork)
            break
preparing()
