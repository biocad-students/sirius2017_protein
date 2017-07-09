# *-* encoding: utf-8 *-*

import numpy as np
from math import sin, cos

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
	res = vec_mult_vec(np.resize(quat1,3), np.resize(quat2,3)) + np.resize(quat1,3) * quat2[3] + np.resize(quat2,3) * quat1[3]
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

