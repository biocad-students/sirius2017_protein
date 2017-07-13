import numpy as np

def get_residue_number(char_number, helis):
	splt = helis[:char_number+1].split()
	return len(splt)-1, len(splt[-1])-1

def dict_ziper(dictionary):
	keys = sorted(dictionary.keys())
	result = []
	for key in keys:
		result.append(dictionary[key])
	return result

def maxSubstring(string1, string2):
	maxpos = (0, 0, 0)
	matrix = np.zeros((len(string1)+1,len(string2)+1), dtype=np.int)
	for i in range(len(string1)+1):
		for j in range(len(string2)+1):
			if min(i,j) > 0 and string1[i-1] == string2[j-1]:
				matrix[i,j] = matrix[i-1,j-1] + 1
				if matrix[i,j] > maxpos[2]:
					maxpos = (i,j,matrix[i,j])
	substring = ''
	for i in range(maxpos[2]):
		substring = string1[maxpos[0]-1-i]+ substring
	return maxpos[0] - maxpos[2], maxpos[1] - maxpos[2], maxpos[2], substring

def heliSubstring(to_parse, var_number, hidden_num=None):
	to_parse = 'SDVVVVTSIMRPAYYYGM'
	f = open('../sirius_out/helis.txt', 'r')
	helis = f.read()
	f.close()
	f = open('../sirius_out/ids.txt', 'r')
	ids = f.read().split()
	f.close()
	if hidden_num:
		index = ids.index(hidden_num)
		helis = helis.split()
		del(helis[index])
		del(ids[index])
		helis = ' '.join(helis)
	results={}
	while not all([s == '' for s in to_parse.split('%')]):
		out = maxSubstring(helis,to_parse)
		nums = get_residue_number(out[0], helis)
		results[out[1]] = (ids[nums[0]], nums[1], out[2])
		to_parse = to_parse[:out[1]] + '%' * out[2] + to_parse[out[1]+out[2]:]
		#print(to_parse, out)
	return dict_ziper(results)
