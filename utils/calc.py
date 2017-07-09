def distance(cordA,cordB):
    """
         Возвращает расстояние между точками в пространстве
         Параметры:
            cordA, cordB - координаты точек
    """
    return sqrt(pow(cordA[0]-cordB[0],2)+pow(cordA[1]-cordB[1],2)+pow(cordA[2]-cordB[2],2))
