import numpy as np
from Bio import PDB as pdb
from geometry.rotation import rotate_vector, vec_mult_vec #(point1, point2, vec, angle)
from math import acos, sqrt, atan2

def arr(vec):
	return np.array(list(vec))

def get_intersection(dot1, dot2, vec):
		return np.dot(vec - dot1, (dot2 - dot1)) * (dot2-dot1) / np.linalg.norm(dot2 - dot1)

def get_char(b,c):
	if acos(b/sqrt(b**2+c**2)) == acos(b/sqrt(b**2+c**2)):
		return 1
	return -1

def do_CCD(loop, endpoint, feedback=False):	
	MAX_ITERATIONS = 5002
	MAX_DISTANCE = 0.08
	for i in range (MAX_ITERATIONS):
		for i1, residue in enumerate(loop):
			if i1 != 0:

				N_vec = residue['N'].get_vector()
				CA_vec= residue['CA'].get_vector()
				C_vec = residue['C'].get_vector()

				N_angle = pdb.calc_dihedral(endpoint['N'].get_vector(), N_vec, CA_vec, loop[-1:][0]['N'].get_vector())
				CA_angle = pdb.calc_dihedral(endpoint['CA'].get_vector(), N_vec, CA_vec, loop[-1:][0]['CA'].get_vector())
				C_angle = pdb.calc_dihedral(endpoint['C'].get_vector(), N_vec, CA_vec, loop[-1:][0]['C'].get_vector())
				
				O_points = np.array([get_intersection(np.array(list(N_vec)),np.array(list(CA_vec)), np.array(list(loop[-1:][0][x].get_vector()))) for x in ['N', 'CA', 'C']])
				multipliers = np.array([np.linalg.norm(np.array(list(loop[-1:][0][x].get_vector())) - O_points[i2]) for i2,x in enumerate(['N','CA','C'])]) 
				b_diffs = np.array([np.array(endpoint[x].get_coord()) - O_points[i2] for i2,x in enumerate(['N','CA','C'])])
				b_diffs = np.array([vec_mult_vec((np.array(list(CA_vec)) - O_points[x]) / np.linalg.norm(np.array(list(CA_vec)) - O_points[x]), b_diffs[x] / np.linalg.norm(b_diffs[x])) for x in range(3)])
				
				b = -2 * sum([multipliers[x] * np.dot(O_points[x], b_diffs[x]) for x in range(3)])
				c = -2 * sum([multipliers[i] * np.dot((np.array(endpoint[x].get_coord()) - O_points[i]) / np.linalg.norm(np.array(endpoint[x].get_coord()) - O_points[i]), O_points[i])for i, x in enumerate(['N','CA','C'])])
				
				angle = -atan2(c/sqrt(b**2+c**2),b/sqrt(b**2+c**2))
				
				
				distangle = {N_angle:N_dist, CA_angle:CA_dist, C_angle:C_dist}
				angles = {abs(i):i for i in [N_angle, CA_angle, C_angle] if i !=0 and distangle[i]>MAX_DISTANCE}
				angles = {abs(i):i for i in [N_angle, CA_angle, C_angle] if i !=0 and distangle[i]>MAX_DISTANCE/1.2}
				if angles.keys():
					angle = angles[min(angles.keys())]
				else:
					angle = 0
				residues = [[atom for atom in residue if not atom.get_name() == 'H']] + loop[i1+1:] #Сначала φ-вращение
				for res in residues:
					for atom in res:
						
						coord = atom.get_vector()
						
						atom.set_coord(rotate_vector(arr(N_vec), arr(CA_vec), arr(coord), -angle))

				N_angle = pdb.calc_dihedral(endpoint['N'].get_vector(), CA_vec, C_vec, loop[-1:][0]['N'].get_vector())
				CA_angle = pdb.calc_dihedral(endpoint['CA'].get_vector(), CA_vec, C_vec, loop[-1:][0]['CA'].get_vector())
				C_angle = pdb.calc_dihedral(endpoint['C'].get_vector(), CA_vec, C_vec, loop[-1:][0]['C'].get_vector())
				

				O_points = np.array([get_intersection(np.array(list(N_vec)),np.array(list(CA_vec)), np.array(list(loop[-1:][0][x].get_vector()))) for x in ['N', 'CA', 'C']])

				multipliers = np.array([np.linalg.norm(np.array(list(loop[-1:][0][x].get_vector())) - O_points[i2]) for i2,x in enumerate(['N','CA','C'])]) 



				b_diffs = np.array([np.array(endpoint[x].get_coord()) - O_points[i2] for i2,x in enumerate(['N','CA','C'])])
				

				b_diffs = np.array([vec_mult_vec((np.array(list(CA_vec)) - O_points[x]) / np.linalg.norm(np.array(list(CA_vec)) - O_points[x]), b_diffs[x] / np.linalg.norm(b_diffs[x])) for x in range(3)])
				

				b = -2 * sum([multipliers[x] * np.dot(O_points[x], b_diffs[x]) for x in range(3)])
				c = -2 * sum([multipliers[i] * np.dot((np.array(endpoint[x].get_coord()) - O_points[i]) / np.linalg.norm(np.array(endpoint[x].get_coord()) - O_points[i]), O_points[i])for i, x in enumerate(['N','CA','C'])])
				
				angle = -atan2(c/sqrt(b**2+c**2),b/sqrt(b**2+c**2))
				
				distangle = {N_angle:N_dist, CA_angle:CA_dist, C_angle:C_dist}
				angles = {abs(i):i for i in [N_angle, CA_angle, C_angle] if i !=0 and distangle[i]>MAX_DISTANCE/1.2}
				if angles.keys():
					angle = angles[min(angles.keys())]
				else:
					angle = 0
				residues = [[residue['O']]] + loop[i1+1:] 
				
				for res in residues:
					
					for atom in res:
						coord = atom.get_vector()
						atom.set_coord(rotate_vector(arr(CA_vec), arr(C_vec), arr(coord), -angle))
						
			N_dist = np.linalg.norm(np.array(list(loop[-1:][0]['N'].get_vector()))-np.array(list(endpoint['N'].get_coord())))
			CA_dist = np.linalg.norm(np.array(list(loop[-1:][0]['CA'].get_vector()))-np.array(list(endpoint['CA'].get_coord())))
			C_dist = np.linalg.norm(np.array(list(loop[-1:][0]['C'].get_vector()))-np.array(list(endpoint['C'].get_coord())))
		 if feedback:
			print(i, sum([N_dist**2, CA_dist**2, C_dist**2]))
		if N_dist < MAX_DISTANCE and C_dist < MAX_DISTANCE and CA_dist < MAX_DISTANCE:
			return loop
	return None

	
