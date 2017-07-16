def histo(way, param):
    f = open(way)
    a = open('results_hist.txt', 'w')
    fl = []
    maxy = 0
    for line in f:
        st = line
        mas = st.split()
        fl.append(float(mas[1]))
        if float(mas[1]) > maxy:
            maxy = float(mas[1])
    n = 0
    n = float(n)
    while n < maxy * 10 + 0.1:
        k = 0
        ch = 0
        while k < len(fl):
            if float(fl[k]) >= n/10 and float(fl[k]) < n/10 + param:
                ch += 1
            k += 1
        n += param * 10
        a.write (str(n/10) + ' ' + str(ch) + '\n')
    f.close()
    a.close()
