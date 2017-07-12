from Bio.PDB import *
from Bio.PDB.Chain import Chain

def final_check(path1,path2,path3):
    """
        Читает 2 PDB файла и считает RMSD двух петлей
        path1 - путь к таблице петлей
        path2 - папка с нашими результатами
        path3 - папка с эталонами

    cdr3 хранит: имя файла | нач. индекс петельки | кон. индекс петельки
    """
    cdr3 = []
    file1 = open(path1, "r")
    _tmp = file1.read()
    file1.close()
    lines = _tmp.split('\n')
    tmp = []
    for x in lines:
        tmp.append(x.split(' '))
    cordstart = len(tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[5])
    cordstop = cordstart+len(tmp[6])+1
    cdr3.append([tmp[0][0],cordstart,cordstop])


    ourfile = open(path2+str(cdr3[0][0])+".pdb")
    ourtext = ourfile.read()
    ourfile.close()

    notourfile = open(path3+str(cdr3[0][0])+".pdb")
    notourtext = ourfile.read()
    notourfile.close()

    ourchain = Chain(0)
    notourchain = Chain(1)

    print(ourchain)
    print(notourchain)
