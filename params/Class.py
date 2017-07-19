#требуется обьявить обьект с идентификатором
f = open('../sampling/kmersnd.txt')
class MyClass:
    def __init__(self, num):
        mas = []
        masnext = []
        mas2 = []
        ind = 'TAS'
        for line in f:
            st = line
            mas = st.split()
            if ind == mas[0]:
                n = 1
                while n < len(mas):
                    mas2.append(float(mas[n]))
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
        self.pos = arr[num]
