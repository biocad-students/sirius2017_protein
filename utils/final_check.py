from Bio.PDB import *
from Bio.PDB.Chain import Chain
from utils.calc import rmsd
from utils.io import *

def final_check(path1,path2,path3,path4):
    """
        Читает 2 PDB файла и считает RMSD двух петлей
        path1 - путь к таблице петлей
        path2 - папка с нашими результатами
        path3 - папка с эталонами
        path4 - файл с результатами сравнения
        extrashift - сдвиг нумерации

    cdr3 хранит: имя файла | нач. индекс петельки | кон. индекс петельки
    """
    cdr3 = []
    file1 = open(path1, "r")
    _tmp = file1.read()
    file1.close()
    lines = _tmp.split('\n')
    tmp = []
    result = open(path4,"w")
    for x in lines:
        tmp.append(x.split(' '))
    for x in range(len(tmp)):
        cordstart = len(tmp[x][1]+tmp[x][2]+tmp[x][3]+tmp[x][4]+tmp[x][5])
        cordstop = cordstart+len(tmp[x][6])
        cdr3.append([tmp[x][0],cordstart,cordstop])
    for counter in range(0,len(cdr3)):
        try:
            ourpdb = read(path2+str(cdr3[counter][0])+".pdb")
            ourres = ourpdb.child_list
            notourpdb = read(path3+str(cdr3[counter][0])+".pdb")
            notourres = notourpdb.child_list
            if(ourres[0].resname == 'ACE'):
                ourshift = 1
            else:
                ourshift = 0
            if(notourres[0].resname == 'ACE'):
                notourshift = 1
            else:
                notourshift = 0
            ourchain = Chain(0)
            notourchain = Chain(1)
            for x in range(cdr3[counter][1],cdr3[counter][2]):
                ourchain.add(ourres[x+ourshift])
            for x in range(cdr3[counter][1],cdr3[counter][2]):
                notourchain.add(notourres[x+notourshift])

            rms = rmsd(ourchain,notourchain)

            result.write(str(cdr3[counter][0])+" "+str(rms)+"\n")
            print("Done: ",)
        except:
            print("Error: reading ",str(cdr3[counter][0]))
    result.close()
