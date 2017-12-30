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

def candidate_search(filename, candidate):
	proteins = []
	lines = reader(filename)
	grouped_proteins = group_by_protein(lines)
	for protein in grouped_proteins:
		protein_name = protein[0]
		if search_in_DNA(protein[1], candidate):
			proteins.append(protein_name)
	return proteins

def search_in_DNA(DNA, candidate):
	concatenated_DNA = ''
	for line in DNA:
		concatenated_DNA += line
	return True if candidate in concatenated_DNA else False

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

def index_candidate(DNA, candidate):
	concatenated_DNA = ''
	for line in DNA:
		concatenated_DNA += line
	index = 0
	for nucleotide in concatenated_DNA:
		if concatenated_DNA[index:index+len(candidate)] == candidate:
			return index
		index+=1

def is_same_location(filenames, genes, candidate):
	all_bools = []
	for gene in genes:
		indices = []
		for filename in filenames:
			lines = reader(filename)
			grouped_proteins = group_by_protein(lines)
			for protein in grouped_proteins:
				if gene in protein[0]:
					if search_in_DNA(protein[1], candidate):
						location = index_candidate(protein[1], candidate)
						indices.append(location)
						break
			'HANDLE CASES IN WHICH YOU HAVE MULTIPLE PROTEINS OF THE SAME NAME'
		default_bool = True
		print(indices, gene)
		for index in indices[1:]:
			if index != indices[0]:
				default_bool = False
		all_bools.append(default_bool)
	return all_bools