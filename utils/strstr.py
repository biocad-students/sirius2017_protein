import numpy as np


def get_residue_number(char_number, loops):
	splt = loops[:char_number+1].split()
	return len(splt)-1, len(splt[-1])-1

def dict_ziper(dictionary):
	keys = sorted(dictionary.keys())
	result = []
	for key in keys:
		result.append(dictionary[key])
	return result

def maxSubstring(string1, string2, lenner):
	maxpos = (0, 0, 0)
	matrix = np.zeros((len(string1)+1,len(string2)+1), dtype=np.int)
	for i in range(len(string1)+1):
		for j in range(len(string2)+1):
			if min(i,j) > 0 and string1[i-1] == string2[j-1]:
				matrix[i,j] = matrix[i-1,j-1] + 1
				if matrix[i,j] - lenner > maxpos[2]:
					maxpos = (i,j,matrix[i,j]- lenner)
	return maxpos[0] - maxpos[2], maxpos[1] - maxpos[2], maxpos[2]

def loopSubstring(to_parse, var_number, hidden_num=None):
	to_parse = 'SDVVVVTSIMRPAYYYGM'
	f = open('../sirius_out/loops.txt', 'r')
	loops = f.read()
	f.close()
	f = open('../sirius_out/ids.txt', 'r')
	ids = f.read().split()
	f.close()
	if hidden_num:
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
		while not all([s == '' for s in parsing.split('%')]):
			out = maxSubstring(loops,parsing, min(max(i-i2,0), max(max([len(s) for s in parsing.split('%')]) - 2, 0)))
			nums = get_residue_number(out[0], loops)
			results[i][out[1]] = (ids[nums[0]], nums[1], out[2])
			parsing = parsing[:out[1]] + '%' * out[2] + parsing[out[1]+out[2]:]
			i2+=1
		results[i] = dict_ziper(results[i])
	return results

