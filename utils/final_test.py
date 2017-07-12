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
    file = open(path1, "r")
    _tmp = 1
    while(_tmp):
        _tmp = file.readline()
        tmp = _tmp.split(' ')
        print(tmp)
        cordstart = len(tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[5])
        #cordstop = cordstart+len(tmp[6])+1
        #cdr3.append([tmp[0].split(' '),cordstart,cordstop])
        if(tmp[0] == '15648'):
            break
    ourfile = read(path2+cdr3[0][0]+".pdb")
    notourfile = read(path3+cdr3[0][0]+".pdb")

    ourchain = Chain(0)
    notourchain = Chain(1)
    ourresidues= ourfile.get_list()
    notourresidues= notourfile.get_list()
    for x in range(cdr3[0][1],cdr3[0][2]):
        ourchain.add(ourresidues[x])
    for x in range(cdr3[0][1],cdr3[0][2]):
        notourchain.add(notourresidues[x])
    print(ourchain)
    print(notourchain)
