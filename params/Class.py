#требуется обьявить обьект с идентификатором
f = open('../files/kmers.txt')
class MyClass:
    def __init__(self):
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
        self.pos = arr

    def getValue(self,key):
        return self.pos[key]
