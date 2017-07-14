from Bio.PDB import *
from Bio.PDB.Chain import Chain

def getletter(structure):
    """
        Достаёт букву из структуры
        Параметры:
            structure - структура
    """
    return structure[0].child_list[0].id

def read(path,name = "test"):
    """
        Читает цепочку аминокислот
        Параметры:
            path - путь к файлу
            name - имя структуры
    """
    parser = PDBParser()
    structure = parser.get_structure(name,path)
    return structure[0][getletter(structure)]

def write(path,record):
    """
        Записывает цепочку в файл
        Параметры:
            path - путь к файлу
            record - структура
    """
    io = PDBIO()
    io.set_structure(record)
    io.save(path)

def addDash(s,index):
    """
        Вставляет тире в строке по индексу
        Параметры:
            s - строка
            index - индекс
    """
    return s[:index] + '-' + s[index:]

def tta(tup,extra = 0):
    a = [' ']
    for x in tup:
        a.append(x)
    a.append(' ')
    return a

def compileres(residues,id = 'H'):
    """
        Собирает список из residue в chain
        Параметры:
            residues - список residue
            id - идентификатор цепочки
    """
    chain = Chain("H")
    startid = 0
    for x in residues:
        x.detach_parent()
        x.id = (' ',startid,' ')
        startid+=1
        chain.add(x)
    return chain

def writeres(path,residues,id = 'H'):
    write(path,compileres(residues,id))
