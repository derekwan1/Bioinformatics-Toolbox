def candidate_search(filename, candidate):
	proteins = []
	lines = reader(filename)
	grouped_proteins = group_by_protein(lines)
	for protein in grouped_proteins:
		protein_name = protein[0]
		if search_in_DNA(protein[1], candidate):
			proteins.append(protein)
	return proteins

def hypotheticals_with_C(proteins):
	hypotheticals_that_have_C = []
	for elem in proteins:
		if 'hypothetical' in elem[0]:
			hypotheticals_that_have_C.append(elem)
	return hypotheticals_that_have_C

def search_in_DNA(DNA, candidate):
	concatenated_DNA = ''
	for line in DNA:
		concatenated_DNA += line
	return True if candidate in concatenated_DNA else False

def reader(filename):
    f = open(filename, 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()
    return lines

def get_protein_indices(lines):
	index = 0
	indices = []
	for line in lines:
		if '>' in line:
			indices.append(index)
		index += 1
	return indices

def group_by_protein(lines):
	index = 0 
	correction_factor = 0
	indices = get_protein_indices(lines)
	all_proteins = []
	while index <= len(indices) - 1:
		curr_protein = []
		protein_name = lines[0][1:]
		curr_protein.append(protein_name)
		if len(indices) == index + 1:
			curr_protein.append(lines[1:])
		else:
			next_protein_index = indices[index+1] - correction_factor
			curr_protein.append(lines[1:next_protein_index])
			correction_factor += next_protein_index
			lines[0:next_protein_index] = []
		index+=1
		all_proteins.append(curr_protein)
	return all_proteins

def concat(lst):
	string = ''
	for elem in lst:
		string += elem
	return string

def compare(a, b):
	z = 0
	while z<=len(a) - 1:
		elem1 = a[z]
		name1 = elem1[0]
		sequence1 = concat(elem1[1])
		k = 0
		while k<=len(b) - 1:
			elem2 = b[k]
			name2 = elem2[0]
			sequence2 = concat(elem2[1])
			print(sequence2)
			break
			if percent_similarity(sequence1, sequence2) >= 0.99:
				elem1.append(name2)
				hide = b.pop(k)
				break
			elif k == len(b) - 1:
				hide = a.pop(z)
				z -= 1
				break
			k += 1
		z += 1
	return a

def percent_similarity(string1, string2):
	k, minimum, count = 0, min(len(string1), len(string2)), 0
	while k<=minimum -1:
		if string1[k] == string2[k]:
			count += 1
		k += 1
	return count / minimum

