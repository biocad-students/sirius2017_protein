import numpy as np
from Bio import PDB as pdb
from geometry.quaternion import Quaternion as Quat
from math import sqrt, atan2, sin, cos
from utils.calc import normalize

lst = list
summa = sum

def get_align(loop, endpoint, char): #Отстояние атомa char конца петли от его конечной позиции
	return  np.linalg.norm(arr(loop[-1][char].get_vector())-arr(endpoint[char].get_coord()))

def norm(vec): #Нормировка вектора по длине 1
	if summa(lst(vec)) > 0:
		return arr(vec)/np.linalg.norm(arr(vec))
	return arr(vec)

def arr(vec):  #Конвертация из Vector Biopython в array[3] numpy
	return np.array(lst(vec))

def get_intersection(dot1, dot2, vec): #Нахождение точки пересечения прямой (задаётся dot1 и dot2) и перпендикуляра из vec
		return np.dot(vec - dot1, norm(dot2-dot1)) * norm(dot2-dot1) + dot1

def calc_O(dot1, dot2, loop, char):
	return get_intersection(arr(dot1),arr(dot2), arr(loop[-1][char].get_vector()))

def calc_multiplier(loop, O_point, char):
	return np.linalg.norm(arr(loop[-1][char].get_vector()) - O_point)

def calc_c(loop, endpoint, dot1, dot2, multiplier, O_point, char, diff):
	c_diff = np.cross(diff, norm(dot2-dot1))	#si
	c_diff = np.dot(c_diff, arr(endpoint[char].get_vector()) - O_point)	#si * OiFi

	c = multiplier * c_diff #ci = -2 * (si * OiFi)
	return c

def calc_b(loop, endpoint, dot1, dot2, multiplier, O_point, char, diff):	# ri
	b_diff = np.dot(diff, arr(endpoint[char].get_vector()) - O_point)	#ri * OiFi
	b =  multiplier * b_diff #bi = -2 * (ri * OiFi)
	return b

def do_rotation(loop, endpoint, res_index, char1, char2):
	dot1 = loop[res_index][char1].get_vector()
	dot2= loop[res_index][char2].get_vector()

	b, c = 0., 0.
	for char in ['N','CA','C']:
		O_point = calc_O(dot1, dot2, loop, char)
		difference = norm(loop[-1][char].get_vector() - O_point)
		multiplier = calc_multiplier(loop, O_point, char)
		c += 2 * calc_c(loop, endpoint, dot1, dot2, multiplier, O_point, char, difference)
		b += 2 * calc_b(loop, endpoint, dot1, dot2, multiplier, O_point, char, difference)
	sqt = sqrt(b**2+c**2)

	if char1 == 'N':
		residues = [[atom for atom in loop[res_index] if not atom.get_name() == 'H']] #φ-вращение
	else:
		residues = [[loop[res_index]['O']]]	#φ-вращение
	residues = residues + loop[res_index+1:]	#ψ-вращение

	angle = atan2(-c/sqt,b/sqt)
	sind = sin(angle/2)
	cosd = cos(angle/2)

	dot1, dot2 = dot1.get_array(), dot2.get_array()
	q = Quat.from_axistrig(sind, cosd, dot2-dot1)
	for res in residues:
		for atom in res:
			coord = atom.get_vector().get_array()
			atom.set_coord(q * (coord - dot1)+dot1)
			#atom.set_coord(matrix_axan_rotation(norm(dot2-dot1), coord, -c/sqt,b/sqt, Eab#))
			#atom.set_coord(cpp_rot(dot1, dot2, coord, -c/sqt, b/sqt, Eab))
			#atom.set_coord(calcnewcord(dot1, dot2, coord, -c/sqt,b/sqt,Eab))


def CCD(loop, endpoint, feedback=False, save_every_loop=False):
	if feedback: #Отладочная информация
		N_dist = get_align(loop, endpoint, 'N')
		CA_dist = get_align(loop, endpoint, 'CA')
		C_dist = get_align(loop, endpoint, 'C')
		print(0, summa([N_dist**2, CA_dist**2, C_dist**2]))
	MAX_ITERATIONS = 5002
	MAX_DISTANCE = 0.08	#Максимальное расстояние между атомами
	for i in range (MAX_ITERATIONS):
		for i1, residue in enumerate(loop):
			if i>0:	#Двигаем ото всех кроме первого residue
				do_rotation(loop, endpoint, i1, 'N', 'CA') #Сначала φ-вращение
				do_rotation(loop, endpoint, i1, 'CA', 'C') #Затем ψ-вращение

		N_dist = get_align(loop, endpoint, 'N')
		CA_dist = get_align(loop, endpoint, 'CA')
		C_dist = get_align(loop, endpoint, 'C')
		if feedback: #Отладочная информация
			print(i, summa([N_dist**2, CA_dist**2, C_dist**2]))
		if N_dist < MAX_DISTANCE and C_dist < MAX_DISTANCE and CA_dist < MAX_DISTANCE:
			return loop
	return None
