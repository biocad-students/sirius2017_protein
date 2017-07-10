def movechain(structure,coord):
    """
    Сдвигает все атомы на cord
    Параметры:
        structure - цепочка
        coord - вектор сдвига
    """
    for spx in range(j,len(structure)):
        for i in structure[spx]:
                vect = structure[spx][i.get_name()].get_vector()+coord
                structure[spx][i.get_name()].set_coord(vect)
    return structure

def combine(structure1,structure2):
    
