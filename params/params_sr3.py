#Случайные углы
from os import listdir
from Bio.SeqUtils import *
from Bio.PDB import *
from math import *
import numpy as np
import random
kolvo = []
def varrand(mas, b):
    def rand(am1, am2, j):
        Pam = (am1 + ' ' + am2)
        f = open('text.txt')
        q = 0
        arr = []
        for line in f:
            string = line
            if line.find(Pam) != -1:
                arr.append(string)
                q += 1
        p = 1
        while p <= int(j):
            r = random.randint(0, q - 1)
            g = arr[r].split()
            an = []
            an.append(g[2])
            an.append(g[4])
            an.append(g[3])
            kolvo.append(an)
            p += 1
        f.close()
    z = 1
    x = 0
    itog = [[[1, 2, 3]] * (len(mas) - 1) for o in range(int(b))]
    while z < len(mas):
        rand(mas[z - 1], mas[z], int(b))
        while x < int(b):
            itog [x] [z - 1] = kolvo[x + (z - 1) * int(b)]
            x += 1
        z += 1
        x = 0
    return (itog)
varrand ('QQQ', 7)
