ffrom Bio.SeqUtils import *
from Bio.PDB import *
from math import *
import numpy as np
def atom_distance(arr):
    nar = []
    mas = {}
    mas['CC'] = 2.37
    mas['CH'] = 1.93
    mas['HC'] = 1.93
    mas['CN'] = 2.29
    mas['NC'] = 2.29
    mas['CO'] = 2.60
    mas['OC'] = 2.60
    mas['NH'] = 1.75
    mas['HN'] = 1.75
    mas['OH'] = 1.51
    mas['HO'] = 1.51
    mas['NO'] = 2.18
    mas['ON'] = 2.18
    mas['HH'] = 1.27
    mas['OO'] = 2.47
    mas['NN'] = 2.57
    n = 2
    while n < len(arr):
        amin = arr[n]
        if str(seq1(str(amin.get_resname()))) != 'P':
            H = amin['H']
            k = n + 1
            while k < len(arr):
                amin2 = arr[k]
                N1 = amin2['N']
                if H - N1 < mas['HN']:
                    nar.append(n)
                    nar.append(k)
                CA1 = amin2['CA']
                if H - CA1 < mas['HC']:
                    nar.append(n)
                    nar.append(k)
                if str(seq1(str(amin2.get_resname()))) == 'G':
                    HA8 = amin2['HA1']
                    if H - HA8 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                    HA9 = amin2['HA2']
                    if H - HA9 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                else:
                    HA1 = amin2['HA']
                    if H - HA1 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                C1 = amin2['C'] 
                if H - C1 < mas['HC']:
                    nar.append(n)
                    nar.append(k)
                O1 = amin2['O']  
                if H - O1 < mas['HO']:
                    nar.append(n)
                    nar.append(k)
                if str(seq1(str(amin2.get_resname()))) != 'P':   
                    H1 = amin2['H']
                    if H - H1 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                k += 1
        N = amin['N']
        k = n + 1
        while k < len(arr):
            amin2 = arr[k]
            N1 = amin2['N']
            if N - N1 < mas['NN']:
                nar.append(n)
                nar.append(k)
            CA1 = amin2['CA']
            if N - CA1 < mas['NC']:
                nar.append(n)
                nar.append(k)
            if str(seq1(str(amin2.get_resname()))) == 'G':
                    HA8 = amin2['HA1']
                    if N - HA8 < mas['NH']:
                        nar.append(n)
                        nar.append(k)
                    HA9 = amin2['HA2']
                    if N - HA9 < mas['NH']:
                        nar.append(n)
                        nar.append(k)
            else:
                    HA1 = amin2['HA']
                    if N - HA1 < mas['NH']:
                        nar.append(n)
                        nar.append(k)
            C1 = amin2['C']
            if N - C1 < mas['NC']:
                nar.append(n)
                nar.append(k)
            O1 = amin2['O']
            if N - O1 < mas['NO']:
                nar.append(n)
                nar.append(k)
            if str(seq1(str(amin2.get_resname()))) != 'P':   
                H1 = amin2['H']
                if N - H1 < mas['NH']:
                    nar.append(n)
                    nar.append(k)
            k += 1
        CA = amin['CA']
        k = n + 1
        while k < len(arr):
            amin2 = arr[k]
            N1 = amin2['N']
            if CA - N1 < mas['CN']:
                nar.append(n)
                nar.append(k)
            CA1 = amin2['CA']
            if CA - CA1 < mas['CC']:
                nar.append(n)
                nar.append(k)
            if str(seq1(str(amin2.get_resname()))) == 'G':
                    HA8 = amin2['HA1']
                    if CA - HA8 < mas['CH']:
                        nar.append(n)
                        nar.append(k)
                    HA9 = amin2['HA2']
                    if CA - HA9 < mas['CH']:
                        nar.append(n)
                        nar.append(k)
            else:
                    HA1 = amin2['HA']
                    if CA - HA1 < mas['CH']:
                        nar.append(n)
                        nar.append(k)
            C1 = amin2['C'] 
            if CA - C1 < mas['CC']:
                nar.append(n)
                nar.append(k)
            O1 = amin2['O']  
            if CA - O1 < mas['CO']:
                nar.append(n)
                nar.append(k)
            if str(seq1(str(amin2.get_resname()))) != 'P':   
                H1 = amin2['H']
                if CA - H1 < mas['CH']:
                    nar.append(n)
                    nar.append(k)
            k += 1
        if str(seq1(str(amin.get_resname()))) == 'G':
            HA3 = amin['HA1']
            k = n + 1
            while k < len(arr):
                amin2 = arr[k]
                N1 = amin2['N']
                if HA3 - N1 < mas['HN']:
                    nar.append(n)
                    nar.append(k)
                CA1 = amin2['CA']
                if HA3 - CA1 < mas['HC']:
                    nar.append(n)
                    nar.append(k)
                if str(seq1(str(amin2.get_resname()))) == 'G':
                    HA8 = amin2['HA1']
                    if HA3 - HA8 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                    HA9 = amin2['HA2']
                    if HA3 - HA9 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                else:
                    HA1 = amin2['HA']
                    if HA3 - HA1 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                C1 = amin2['C']  
                if HA3 - C1 < mas['HC']:
                    nar.append(n)
                    nar.append(k)
                O1 = amin2['O'] 
                if HA3 - O1 < mas['HO']:
                    nar.append(n)
                    nar.append(k)
                if str(seq1(str(amin2.get_resname()))) != 'P':   
                    H1 = amin2['H']
                    if HA3 - H1 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                k += 1
            HA4 = amin['HA2']
            k = n + 1
            while k < len(arr):
                amin2 = arr[k]
                N1 = amin2['N']
                if HA4 - N1 < mas['HN']:
                    nar.append(n)
                    nar.append(k)
                CA1 = amin2['CA']
                if HA4 - CA1 < mas['HC']:
                    nar.append(n)
                    nar.append(k)
                if str(seq1(str(amin2.get_resname()))) == 'G':
                    HA8 = amin2['HA1']
                    if HA4 - HA8 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                    HA9 = amin2['HA2']
                    if HA4 - HA9 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                else:
                    HA1 = amin2['HA']
                    if HA4- HA1 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                C1 = amin2['C']  
                if HA4- C1 < mas['HC']:
                    nar.append(n)
                    nar.append(k)
                O1 = amin2['O'] 
                if HA4- O1 < mas['HO']:
                    nar.append(n)
                    nar.append(k)
                if str(seq1(str(amin2.get_resname()))) != 'P':   
                    H1 = amin2['H']
                    if HA4- H1 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                k += 1
        else:
            HA = amin['HA']
            k = n + 1
            while k < len(arr):
                amin2 = arr[k]
                N1 = amin2['N']
                if HA - N1 < mas['HN']:
                    nar.append(n)
                    nar.append(k)
                CA1 = amin2['CA']
                if HA - CA1 < mas['HC']:
                    nar.append(n)
                    nar.append(k)
                if str(seq1(str(amin2.get_resname()))) == 'G':
                    HA8 = amin2['HA1']
                    if HA- HA8 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                    HA9 = amin2['HA2']
                    if HA- HA9 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                else:
                    HA1 = amin2['HA']
                    if HA- HA1 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                C1 = amin2['C']  
                if HA - C1 < mas['HC']:
                    nar.append(n)
                    nar.append(k)
                O1 = amin2['O'] 
                if HA - O1 < mas['HO']:
                    nar.append(n)
                    nar.append(k)
                if str(seq1(str(amin2.get_resname()))) != 'P':   
                    H1 = amin2['H']
                    if HA - H1 < mas['HH']:
                        nar.append(n)
                        nar.append(k)
                k += 1
        C = amin['C']
        k = n + 1
        while k < len(arr):
            amin2 = arr[k]
            N1 = amin2['N']
            if C - N1 < mas['CN']:
                if k - n != 1:
                    print(C - N1)
                    nar.append(n)
                    nar.append(k)
            CA1 = amin2['CA']
            if C - CA1 < mas['CC']:
                nar.append(n)
                nar.append(k)
            if str(seq1(str(amin2.get_resname()))) == 'G':
                    HA8 = amin2['HA1']
                    if C- HA8 < mas['CH']:
                        nar.append(n)
                        nar.append(k)
                    HA9 = amin2['HA2']
                    if C - HA9 < mas['CH']:
                        nar.append(n)
                        nar.append(k)
            else:
                    HA1 = amin2['HA']
                    if C - HA1 < mas['CH']:
                        nar.append(n)
                        nar.append(k)
            C1 = amin2['C'] 
            if C - C1 < mas['CC']:
                nar.append(n)
                nar.append(k)
            O1 = amin2['O']  
            if C - O1 < mas['CO']:
                nar.append(n)
                nar.append(k)
            if str(seq1(str(amin2.get_resname()))) != 'P':   
                H1 = amin2['H']
                if C - H1 < mas['CH']:
                    nar.append(n)
                    nar.append(k)
            k += 1
        O = amin['O']
        k = n + 1
        while k < len(arr):
            amin2 = arr[k]
            N1 = amin2['N']
            if O - N1 < mas['ON']:
                nar.append(n)
                nar.append(k)
            CA1 = amin2['CA']
            if O - CA1 < mas['OC']:
                nar.append(n)
                nar.append(k)
            if str(seq1(str(amin2.get_resname())))== 'G':
                    HA8 = amin2['HA1']
                    if O - HA8 < mas['OH']:
                        nar.append(n)
                        nar.append(k)
                    HA9 = amin2['HA2']
                    if O - HA9 < mas['OH']:
                        nar.append(n)
                        nar.append(k)
            else:
                    HA1 = amin2['HA']
                    if O - HA1 < mas['OH']:
                        nar.append(n)
                        nar.append(k)
            C1 = amin2['C']   
            if O - C1 < mas['OC']:
                nar.append(n)
                nar.append(k)
            O1 = amin2['O'] 
            if O - O1 < mas['OO']:
                nar.append(n)
                nar.append(k)
            if str(seq1(str(amin2.get_resname()))) != 'P':   
                H1 = amin2['H']
                if O - H1 < mas['OH']:
                    nar.append(n)
                    nar.append(k)
            k += 1
            
        n += 1
    l = 0
    itog = []
    try:
        value = nar[0]
    except IndexError:
        return(itog)
    else:
        while l < len(nar):
            s = 0
            ch = 0
            while s < l:
                if nar[l] == nar[s]:
                       ch += 1
                s += 1
            if ch == 0:
                itog.append(nar[l])
            l += 1
        return(itog)
