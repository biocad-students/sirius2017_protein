def movechain(structure,coord,extra = 0):
    """
        Сдвигает все атомы на cord
        Параметры:
            structure - цепочка
            coord - вектор сдвига
            extra1 - сдвиг элементов
    """
    for spx in range(extra,len(structure)+extra):
        for i in structure[spx]:
                vect = structure[spx][i.get_name()].get_vector()+coord
                structure[spx][i.get_name()].set_coord(vect)
    return structure
