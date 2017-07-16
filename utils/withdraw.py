from Bio import PDB as pdb

def regions_to_dict(regions): #Превращает строку файла regions в словарь id:[имена residue]
	result = {}
	for string in regions.split('\n'):
		if string:
			result[string.split()[0]] = string.split()[1:]
	return result

def get_residues_by_pos(resList, prefix):
	'''
		Принимает список кортежей в формате вывода strstr.loopSubstring, возвращает список списков соответствующих residue
	'''
	parser = pdb.PDBParser()
	residues = []
	regions = regions_to_dict(open(prefix + 'regions.txt', 'r').read())
	for res_tuple in resList:
		structure = parser.get_structure('strk', prefix+res_tuple[0]+'.pdb')
		align = sum([ len(s) for s in regions[res_tuple[0]][:-2]]) + res_tuple[1]
		residues.append(list(structure.get_residues())[align:align+res_tuple[2]+1])
	return residues
