#требуется обьявить обьект с идентификатором (3 буквы)
class MyClass:
    def __init__(self):
        f = open('kmersnd.txt'#вставьте сюда путь до этого файла)
        arr = {}
        mas = []
        masnext = []
        mas2 = []
        ind = 'TAS'
        arr = {}
        for line in f:
            st = line
            mas = st.split(' ')
            if ind == mas[0]:
                n = 1
                while n < len(mas):
                    try:
                        mas2.append(float(mas[n]))
                    except:
                        print("OUT:",mas[n])
                    n += 1
                masnext.append(mas2)
                mas2 = []
            else:
                arr[ind] = masnext
                masnext = []
                ind = mas[0]
                n = 1
                while n < len(mas):
                    mas2.append(float(mas[n]))
                    n += 1
                masnext.append(mas2)
                mas2 = []
        arr[ind] = masnext
        try:
            value = arr[num]
        except KeyError:
            self.pos = []
        else:
            self.pos = arr[num]
