from Bio.PDB import *
def read(path,name = "test"):
    """
        Читает цепочку аминокислот
        Параметры:
            path - путь к файлу
            name - имя структуры
    """
    parser = PDBParser()
    structure = parser.get_structure(name,path)
    return structure[0]['H']

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
