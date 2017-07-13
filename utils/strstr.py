import numpy as np


def get_residue_number(char_number, loops):	#Номер  pdb файла по его индексу в regions.txt
	splt = loops[:char_number+1].split()
	return len(splt)-1, len(splt[-1])-1

def dict_ziper(dictionary):
	keys = sorted(dictionary.keys())
	result = []
	for key in keys:
		result.append(dictionary[key])
	return result

def maxSubstring(string1, string2, lenner): #Максимальная общая подстрока string1 и string2
	maxpos = (0, 0, 0)
	matrix = np.zeros((len(string1)+1,len(string2)+1), dtype=np.int)
	for i in range(len(string1)+1):
		for j in range(len(string2)+1):
			if min(i,j) > 0 and string1[i-1] == string2[j-1]:
				matrix[i,j] = matrix[i-1,j-1] + 1
				if matrix[i,j] - lenner > maxpos[2]: #Уменьшаем потенциальнуя строку на lenner для вариативности
					maxpos = (i,j,matrix[i,j]- lenner)
	return maxpos[0] - maxpos[2], maxpos[1] - maxpos[2], maxpos[2]

def loopSubstring(to_parse, var_number, hidden_num=None):
	'''
		loopSubstring(to_parse, var_number, hidden_num=None)
		Возвращает var_number вариантов разбиения to_parce на подстроки из БД петель
		Возвращает список списков кортежей:
		loopSubstring() = [[(protId, start, len), (protId, start,len)...],[(protId, start, len), (protId, start,len)...],...]
		Где protId -- id белка, start -- первый residue в последовательности, len -- длина последовательности
	'''
	f = open('../sirius_out/loops.txt', 'r')
	loops = f.read()
	f.close()
	f = open('../sirius_out/ids.txt', 'r')
	ids = f.read().split()
	f.close()
	if hidden_num:		#Прячем "скрытый номер"
		index = ids.index(hidden_num)
		loops = loops.split()
		del(loops[index])
		del(ids[index])
		loops = ' '.join(loops)
	results=[]
	for i in range(var_number):
		parsing = to_parse
		results.append({})
		i2 = 0
		while not all([s == '' for s in parsing.split('%')]):	#Жадно разбиваем на подстроки
			out = maxSubstring(loops,parsing, min(max(i-i2,0), max(max([len(s) for s in parsing.split('%')]) - 2, 0)))
			nums = get_residue_number(out[0], loops)
			results[i][out[1]] = (ids[nums[0]], nums[1], out[2])
			parsing = parsing[:out[1]] + '%' * out[2] + parsing[out[1]+out[2]:]
			i2+=1
		results[i] = dict_ziper(results[i])
	return results
#Возвращает список списков (вариантов разбиения) кортежей (каждый из которых имеет)
