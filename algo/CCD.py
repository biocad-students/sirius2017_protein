import numpy as np
from Bio import PDB as pdb
from geometry.rotation import rotate_vector, vec_mult_vec, calcnewcord 
from math import acos, sqrt, atan2, sin, cos


def norm(vec): #Нормировка вектора по длине 1
	return arr(vec)/np.linalg.norm(arr(vec))

def arr(vec):  #Конвертация из Vector Biopython в array[3] numpy
	return np.array(list(vec))

def get_intersection(dot1, dot2, vec): #Нахождение точки пересечения прямой (задаётся dot1 и dot2) и перпендикуляра из vec
		return np.dot(vec - dot1, norm(dot2-dot1)) * norm(dot2-dot1)

def CCD(loop, endpoint, feedback=False):
	if feedback: #Отладочная информация
		N_dist = np.linalg.norm(np.array(list(loop[-1]['N'].get_vector()))-np.array(list(endpoint['N'].get_coord())))
		CA_dist = np.linalg.norm(np.array(list(loop[-1]['CA'].get_vector()))-np.array(list(endpoint['CA'].get_coord())))
		C_dist = np.linalg.norm(np.array(list(loop[-1]['C'].get_vector()))-np.array(list(endpoint['C'].get_coord())))
		print(0, sum([N_dist**2, CA_dist**2, C_dist**2]))	
	MAX_ITERATIONS = 5002	
	MAX_DISTANCE = 0.08	#Максимальное расстояние между атомами
	for i in range (MAX_ITERATIONS):
		for i1, residue in enumerate(loop):
			if i1 != 0:	#Двигаем ото всех кроме первого residue
				N_vec = residue['N'].get_vector()
				CA_vec= residue['CA'].get_vector()
				C_vec = residue['C'].get_vector()

				O_points = np.array([get_intersection(arr(N_vec),arr(CA_vec), arr(loop[-1][x].get_vector())) for x in ['N', 'CA', 'C']])	#Координаты точек пересечения осей и перпедикуляров из атомов последнего residue на них (Oi)
				multipliers = np.array([np.linalg.norm(arr(loop[-1][x].get_vector()) - O_points[i2]) for i2,x in enumerate(['N','CA','C'])]) #|| OiM'i ||

				c_diff = [vec_mult_vec(norm(loop[-1][x].get_vector() - O_points[i2]), norm(CA_vec-N_vec)) for i2, x in enumerate(['N','CA','C'])]	#si

				c_diff = [np.dot(c_diff[i2], arr(endpoint[x].get_vector()) - O_points[i2]) for i2, x in enumerate(['N','CA','C'])]	#si * OiFi
								
				c = -2 * sum([multipliers[i] * c_diff[i] for i in range(3)]) #ci = -2 * (si * OiFi)

				b_diff = [norm(loop[-1][x].get_vector() - O_points[i2]) for i2,x in enumerate(['N','CA','C'])]	# ri
				b_diff = [np.dot(b_diff[i2], arr(endpoint[x].get_vector()) - O_points[i2]) for i2, x in enumerate(['N','CA','C'])]	#ri * OiFi
				b = -2 * sum([multipliers[i] * b_diff[i] for i in range(3)]) #bi = -2 * (ri * OiFi)

				sqt = sqrt(b**2+c**2)

				residues = [[atom for atom in residue if not atom.get_name() == 'H']] + loop[i1+1:] #Сначала φ-вращение. Все атомы последующих residue, все атомы текущего residue кроме H			
				
				for res in residues:
					for atom in res:
						
						coord = atom.get_vector()
						atom.set_coord(calcnewcord(N_vec, CA_vec, coord, c/sqt,-b/sqt)) 

				N_vec = residue['N'].get_vector()
				CA_vec= residue['CA'].get_vector()
				C_vec = residue['C'].get_vector()

				O_points = np.array([get_intersection(arr(CA_vec),arr(C_vec), arr(loop[-1][x].get_vector())) for x in ['N', 'CA', 'C']])	#Координаты точек пересечения осей и перпедикуляров из атомов последнего residue на них (Oi)

				multipliers = np.array([np.linalg.norm(arr(loop[-1][x].get_vector()) - O_points[i2]) for i2,x in enumerate(['N','CA','C'])])	#|| OiM'i ||

				c_diff = [vec_mult_vec(norm(loop[-1][x].get_vector() - O_points[i2]), norm(C_vec-CA_vec)) for i2, x in enumerate(['N','CA','C'])]	#si

				c_diff = [np.dot(c_diff[i2], arr(endpoint[x].get_vector()) - O_points[i2]) for i2, x in enumerate(['N','CA','C'])]	#si * OiFi
								
				c = -2 * sum([multipliers[i] * c_diff[i] for i in range(3)])	#ci = -2 * (si * OiFi)

				b_diff = [norm(loop[-1][x].get_vector() - O_points[i2]) for i2,x in enumerate(['N','CA','C'])]	#ri
				b_diff = [np.dot(b_diff[i2], arr(endpoint[x].get_vector()) - O_points[i2]) for i2, x in enumerate(['N','CA','C'])]	#ri * OiFi
				b = -2 * sum([multipliers[i] * b_diff[i] for i in range(3)])	#bi = -2 * (ri * OiFi)

				sqt = sqrt(b**2+c**2)
				residues = [[residue['O']]] + loop[i1+1:] #Затем ψ-вращение. Все атомы последующих residue и O текущего residue
				for res in residues:
					for atom in res:
						coord = atom.get_vector()
						atom.set_coord(calcnewcord(CA_vec, C_vec, coord, c/sqt,-b/sqt))				
		N_dist = np.linalg.norm(arr(loop[-1]['N'].get_vector())-arr(endpoint['N'].get_coord()))
		CA_dist = np.linalg.norm(arr(loop[-1]['CA'].get_vector())-arr(endpoint['CA'].get_coord()))
		C_dist = np.linalg.norm(arr(loop[-1]['C'].get_vector())-arr(endpoint['C'].get_coord()))
		if feedback: #Отладочная информация	
			print(i, sum([N_dist**2, CA_dist**2, C_dist**2]))
		if N_dist < MAX_DISTANCE and C_dist < MAX_DISTANCE and CA_dist < MAX_DISTANCE:
			return loop
	return None

	
