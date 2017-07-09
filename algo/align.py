def calculate(string1,string2):
    """
        Алгоритм Вагнера—Фишера
        Параметры:
            string1,string2 - строки для сравнения
    """
    m = len(string1)+2
    n = len(string2)+2
    matrix = np.zeros((n,m), dtype=np.int)
    for i,x in enumerate(string1):
        matrix[0,i+2] = ord(x)
        matrix[1,i+1] = i
    matrix[1,i+2] = i+1
    for i,x in enumerate(string2):
        matrix[i+2,0] = ord(x)
        matrix[i+1,1] = i
    matrix[i+2,1] = i+1
    for x in range(1,n):
        for y in range(2,m):
            above = matrix[x-1,y]
            side = matrix[x,y-1]
            aaa = matrix[x-1,y-1]
            if(side < aaa or above<aaa):
                outValue = min(above+1,side+1)
            else:
                if(matrix[0,y] != matrix[x,0]):
                    outValue = aaa+1
                else:
                    outValue = aaa
            matrix[x,y] = outValue
    print(matrix)
    string1 = ""
    string2 = ""
    for i in range(2,m):
        string1+=chr(matrix[0,i])
    for i in range(2,n):
        string2+=chr(matrix[i,0])
    x = n-2
    while(x>2):
        y = m-2
        while(y>2):
            above= matrix[x-1,y]
            side = matrix[x,y-1]
            aaa  = matrix[x-1,y-1]
            curr = matrix[x,y]
            if(above+1==curr):
                string1 = addDash(string1,x)
                x-=1
            elif(side+1 == curr):
                string2 = addDash(string2,y-2)
                y-=1
            else:
                x-=1
                y-=1
    print(string1)
    print(string2)
    return string1,string2,matrix[n-2,m-2]
