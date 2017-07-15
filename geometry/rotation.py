# *-* encoding: utf-8 *-*

import numpy as np
from math import sin, cos
from utils.calc import *
from numpy import *

def quat_invert(quat):
	'''
	Инверсия кватерниона 
	quat --> quat
	q^(-1) = [-v, w] / [xx + yy + zz + ww]
	'''
	res = quat * -1
	res[3] = quat[3]
	return res / np.linalg.norm(res)

def quat_mul_quat(quat1, quat2):
	'''
	Произведение кватернионов
	quat, quat --> quat
	qq' = [ vv' + wv' + w'v, ww' – v•v' ]
	'''
	res = cross(np.resize(quat1,3), np.resize(quat2,3)) + np.resize(quat1,3) * quat2[3] + np.resize(quat2,3) * quat1[3]
	res = np.resize(res,4)	
	res[3] = quat1[3] * quat2[3] - np.dot(quat1, quat2)
	return res

def vec_mult_vec(vec1, vec2):
	'''
	Векторное произведение векторов
	vec, vec --> vec
	vv' = [yz' - zy', zx' - xz', xy' - yx']
	'''
	vec = np.zeros(3)
	vec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1]
	vec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2]
	vec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0]
	return vec


##def get_rot_matrix(axis, c, s):
#	"""
#        Генерация матрицы поворота вокруг заданной оси на угол с заданными синусом и косинусом
#        Параметры:
#            axis - вектор оси
#            c -  косинус угла поворота
#            s - синус угла поворота
#    	"""
#    	#axis = axis/math.sqrt(np.dot(axis, axis))
#	pass
	

#def euler_rod_rotation(_veca,_vecb,_vecc, sinangle, cosangle):
#    """
#        Вращение вектора на угол с заданным синусом и косинусом вокруг оси, заданной двумя точками
#        Параметры:
#            _veca,vecb - векторы оси вектора PDBшные
#            _vecc -  начальный вектор
#            sinangle - синус угла поворота
#            cosangle - косинус угла поворота
#    """
#    veca = _veca.get_array()
#    vecb = _vecb.get_array()
#    vecc = _vecc.get_array()



#def rod_calcnewcord(_veca,_vecb,_vecc, sinangle, cosangle):
#    """
#        Вращение вектора на угол с заданным синусом и косинусом вокруг оси, заданной двумя точками
#        Параметры:
#            _veca,vecb - векторы оси вектора PDBшные
#            _vecc -  начальный вектор
#            sinangle - синус угла поворота
#            cosangle - косинус угла поворота
#    """
#	veca = _veca.get_array()
#    	vecb = _vecb.get_array()
#    	vecc = _vecc.get_array()
#	rot_vec = vecb - veca
#	roted_vec = vecc - veca
	



def matrix_axan_rotation(vec,_vecc, s, c):
	"""
        Вращение вектора на угол с заданным синусом и косинусом вокруг оси, заданной двумя точками
        Параметры:
            _veca,vecb - векторы оси вектора PDBшные
            _vecc -  начальный вектор
            s - синус угла поворота
            c - косинус угла поворота
	"""


	vecc = _vecc.get_array()
	
	x, y, z = vec[0], vec[1], vec[2]
	vecca = vecc - vec
	rvec = np.array([[vecca[0]], [vecca[1]], [vecca[2]]])
	#matrix = np.array([ 
	#[cosangle + (1 - cosangle) * x**2, (1-cosangle)*x*y-(sinangle*z), (1-cosangle)*x*z+sinangle*y],
	#[(1 - cosangle)*y*x + sinangle*z, cosangle + (1 - cosangle) * y **2, (1-cosangle)*y*z - sinangle*x],
	#[(1 - cosangle)*z*x - sinangle*y, (1-cosangle)*y*z + sinangle*x, cosangle + (1-cosangle)*z**2]
	#])
	c1 = 1 - c
	matrix = np.array([
	[x**2+c, x*y*c1 - z*s, x*z*c1+y*s],
	[x*y*c1 + z*s, y*y*c1+c, y*z*c1-x*s],
	[x*z*c1 - y*s, x*z*c1 + x*s, z*z*c1 + c]
	])
	print(distance(vec, vecc), distance(vec, np.transpose(matrix @ rvec)[0]))
	return np.transpose(matrix @ rvec)[0] + vecc
	
	Eab = normalize(vecb - veca)
	AC = vecc - veca
	O = veca + Eab * dot(AC,Eab)
	Eoc = normalize(vecc - O)
	S = cross(Eab,Eoc)
	OCd = distance(O,vecc)
	M = Eoc * OCd * cosangle + S * OCd * sinangle + O
	return M
cs = cross
def calcnewcord(_veca,_vecb,_vecc, sinangle, cosangle, Eab):
    """
        Вращение вектора на угол с заданным синусом и косинусом вокруг оси, заданной двумя точками
        Параметры:
            _veca,vecb - векторы оси вектора PDBшные
            _vecc -  начальный вектор
            sinangle - синус угла поворота
            cosangle - косинус угла поворота
    """
    veca = _veca.get_array()
    vecb = _vecb.get_array()
    vecc = _vecc.get_array()

   # Eab = normalize(vecb - veca)
    AC = vecc - veca
    O = veca + Eab * dot(AC,Eab)
    Eoc = normalize(vecc - O)
    S = cs(Eab,Eoc)
    OCd = distance(O,vecc)
    M = Eoc * OCd * cosangle + S * OCd * sinangle + O
    return M
   

def rotate_vector(point1, point2, vec, angle):
	'''
	Вращение вектора на заданный угол вокруг оси, заданной двумя точками
	point1, point2, vec, angle --> vec
	V' = qvq^(–1)
	'''
	new_vec = np.resize(vec - point1,4)
	new_vec[3] = 0
	quat = make_quaternion(point1, point2, angle)
	t = quat_mul_quat(quat, new_vec)
	t = quat_mul_quat(t, quat_invert(quat))
	return np.resize(t,3) + point1

def make_quaternion(point1, point2, angle):
	'''
	Создание кватерниона, заданного осью вращения и углом
	point1, point2, angle --> quat
	q = [cos(a/2), [x sin(a/2), y sin(a/2), z sin(a/2)]]
	'''
	startpoint = point2 - point1
	rotate_vector = startpoint / np.linalg.norm(startpoint)
	quat = np.zeros(4)
	quat[0] = rotate_vector[0] * sin(angle / 2 )
	quat[1] = rotate_vector[1] * sin(angle / 2 )
	quat[2] = rotate_vector[2] * sin(angle / 2 )
	quat[3] = cos(angle / 2 )
	return quat

