import numpy as np
from Bio import PDB as pdb
from geometry.rotation import rotate_vector, vec_mult_vec, calcnewcord
from math import sqrt

def get_align(loop, endpoint, char): #Отстояние атомa char конца петли от его конечной позиции
	return  np.linalg.norm(arr(loop[-1][char].get_vector())-arr(endpoint[char].get_coord()))


def norm(vec): #Нормировка вектора по длине 1
	if sum(list(vec)) > 0:
		return arr(vec)/np.linalg.norm(arr(vec))
	return arr(vec)

def arr(vec):  #Конвертация из Vector Biopython в array[3] numpy
	return np.array(list(vec))

def get_intersection(dot1, dot2, vec): #Нахождение точки пересечения прямой (задаётся dot1 и dot2) и перпендикуляра из vec
		return np.dot(vec - dot1, norm(dot2-dot1)) * norm(dot2-dot1) + dot1

def do_rotation(loop, endpoint, res_index, char1, char2):
	print("res",loop[res_index])
	dot1 = loop[res_index][char1].get_vector()
	dot2= loop[res_index][char2].get_vector()
	O_points = np.array([get_intersection(arr(dot1),arr(dot2), arr(loop[-1][x].get_vector())) for x in ['N', 'CA', 'C']])	#Координаты точек пересечения осей и перпедикуляров из атомов последнего residue на них (Oi)
	multipliers = np.array([np.linalg.norm(arr(loop[-1][x].get_vector()) - O_points[i2]) for i2,x in enumerate(['N','CA','C'])]) #|| OiM'i ||

	c_diff = [vec_mult_vec(norm(loop[-1][x].get_vector() - O_points[i2]), norm(dot2-dot1)) for i2, x in enumerate(['N','CA','C'])]	#si

	c_diff = [np.dot(c_diff[i2], arr(endpoint[x].get_vector()) - O_points[i2]) for i2, x in enumerate(['N','CA','C'])]	#si * OiFi

	c = 2 * sum([multipliers[i] * c_diff[i] for i in range(3)]) #ci = -2 * (si * OiFi)

	b_diff = [norm(loop[-1][x].get_vector() - O_points[i2]) for i2,x in enumerate(['N','CA','C'])]	# ri
	b_diff = [np.dot(b_diff[i2], arr(endpoint[x].get_vector()) - O_points[i2]) for i2, x in enumerate(['N','CA','C'])]	#ri * OiFi
	b = 2 * sum([multipliers[i] * b_diff[i] for i in range(3)]) #bi = -2 * (ri * OiFi)

	sqt = sqrt(b**2+c**2)

	if char1 == 'N':
		residues = [[atom for atom in loop[res_index] if not atom.get_name() == 'H']] #φ-вращение
	else:
		residues = [[loop[res_index]['O']]]	#φ-вращение
	residues = residues + loop[res_index+1:]	#ψ-вращение

	for res in residues:
		for atom in res:
			coord = atom.get_vector()
			atom.set_coord(calcnewcord(dot1, dot2, coord, -c/sqt,b/sqt))


def CCD(loop, endpoint, feedback=False):
	if feedback: #Отладочная информация
		N_dist = get_align(loop, endpoint, 'N')
		CA_dist = get_align(loop, endpoint, 'CA')
		C_dist = get_align(loop, endpoint, 'C')
		print(0, sum([N_dist**2, CA_dist**2, C_dist**2]))
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
			print(i, sum([N_dist**2, CA_dist**2, C_dist**2]))
		if N_dist < MAX_DISTANCE and C_dist < MAX_DISTANCE and CA_dist < MAX_DISTANCE:
			return loop
	return None
