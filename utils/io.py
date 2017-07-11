from Bio.PDB import *
from Bio.PDB.Chain import Chain
def read(path,name = "test"):
    """
        Читает цепочку аминокислот
        Параметры:
            path - путь к файлу
            name - имя структуры
    """
    parser = PDBParser()
    structure = parser.get_structure(name,path)
    return structure[0]['F']

def write(path,record):
    """
        Записывает цепочку в файл
        Параметры:
            path - путь к файлу
            name - имя структуры
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

def compileres(residues):
    """
        Собирает список из residue в chain
        Параметры:
            residues - список residue
    """
    #sb = StructureBuilder()
    chain = Chain(0)
    for x in residues:
        chain.add(x)
    return chain
