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

def common_proteins(all_proteins):
	all_proteins = get_real_name(all_proteins)
	intersections = set(all_proteins[0])
	for protein in all_proteins[1:]:
		intersections = intersections.intersection(set(protein))
	return intersections

def unique_proteins(all_proteins):
	all_proteins = get_real_name(all_proteins)
	unions = set(all_proteins[0])
	for protein in all_proteins[1:]:
		unions = unions.union(set(protein))
	return unions

def get_real_name(all_proteins):
	k1 = 0
	while k1 <= len(all_proteins) - 1:
		k2 = 0
		protein_group = all_proteins[k1]
		while k2 <= len(protein_group) - 1:
			name_index = index_of_actual_name(protein_group[k2])
			protein_group[k2] = protein_group[k2][name_index:]
			k2 += 1
		k1+=1
	return all_proteins

def index_of_actual_name(string):
	index = 0
	while index <= len(string) - 1:
		if string[index] == ' ':
			return index + 1
		else:
			index += 1
	
def hypothetical_protein_names(lines):
	k, hypothetical_proteins = 0, []
	grouped_proteins = group_by_protein(lines)
	while k <= len(grouped_proteins) - 1:
		if 'hypothetical protein' in grouped_proteins[k][0]:
			hypo = grouped_proteins[k]
			b4 = grouped_proteins[k-1][0] if k!= 0 else grouped_proteins[-1][0]
			one_after = grouped_proteins[k+1][0] if k != len(grouped_proteins)-1 else grouped_proteins[0][0]
			if 'hypothetical protein' in one_after and 'hypothetical protein' in b4:
				before, after = 'None', 'None'
				hypothetical_proteins.append([before, hypo, after])
			elif 'hypothetical protein' in b4 and 'hypothetical protein' not in one_after:
				before = 'None'
				full_name_after = grouped_proteins[k+1][0]
				index_after = index_of_actual_name(full_name_after)
				after = full_name_after[index_after:]
				hypothetical_proteins.append([before, hypo, after])
			elif 'hypothetical protein' in one_after and 'hypothetical protein' not in b4:
				after = 'None'
				full_name_before = grouped_proteins[k-1][0]
				index_before = index_of_actual_name(full_name_before)
				before = full_name_before[index_before:]
				hypothetical_proteins.append([before, hypo, after])
			else:
				full_name1 = grouped_proteins[k-1][0]
				full_name2 = grouped_proteins[k+1][0]
				index1 = index_of_actual_name(full_name1)
				index2 = index_of_actual_name(full_name2)
				before = full_name1[index1:]
				after = full_name2[index2:]
				hypothetical_proteins.append([before, hypo, after])
		k+=1
	return hypothetical_proteins

def similarity(seq1, seq2):
	length1 = len_seq(seq1)
	length2 = len_seq(seq2)
	if min(length1, length2) / max(length1, length2) < 0.99:
		return False
	else:
		denom = max(length1, length2)
		k, numerator = 0, 0
		seq1, seq2 = concat(seq1), concat(seq2)
		while k <= min(length1, length2) - 1:
			if seq1[k] == seq2[k]:
				numerator += 1
			k+=1
		percent_similar = numerator/denom
		return True if percent_similar >= 0.99 else False

def len_seq(sequence):
	count = 0
	for line in sequence:
		for nucleotide in line:
			if nucleotide == 'A' or nucleotide == 'a' or nucleotide =='G' or nucleotide == 'g' or nucleotide == 'C' or nucleotide == 'c' or nucleotide == 'T' or nucleotide == 't':
				count+=1
	return count

def concat(seq):
	concatenated_DNA = ''
	for line in seq:
		concatenated_DNA += line
	return concatenated_DNA

def candidate_search_hypothetical(data, candidate):
	new_candidates = []
	for triplet in data:
		if search_in_DNA(triplet[1][1], candidate):
			new_candidates.append(triplet)
	return new_candidates

def protein_matcher(filtered_data):
	dictionary = {}
	for triplet in filtered_data:
		if triplet[0] + triplet[2] not in dictionary.keys():
			dictionary[triplet[0] + triplet[2]] = [triplet[1]]
		else:
			try:
				dictionary[triplet[0] + triplet[2]] = dictionary[triplet[0] + triplet[2]].append(triplet[1])
			except AttributeError:
				if triplet[0] + triplet[2] == '':
					print(triplet[1])
				dictionary[triplet[0] + triplet[2]] = [triplet[1]]
				print('hi', dictionary[triplet[0] + triplet[2]])
	return dictionary

def id_unmatched_proteins(filtered_data, matches):
	unmatched = []
	for triplet in filtered_data:
		matched = False
		for pair in matches:
			if triplet[1][0] in pair:
				matched = True
				break
		if matched == False:
			unmatched.append(triplet[1][0])
	return unmatched

def remove_duplicates(filtered_data):
	no_duplicates_data = []
	for pair1 in filtered_data:
		if pair1 not in no_duplicates_data and [pair1[1], pair1[0]] not in no_duplicates_data:
			no_duplicates_data.append(pair1)
	return no_duplicates_data

def combine_all_matches(no_duplicates_data):
	reference_copy = no_duplicates_data[:]
	all_unions = []
	for pair in no_duplicates_data:
		cont = False
		for elem in no_duplicates_data:
			for sets in all_unions:
				if elem[0] in sets or elem[1] in sets:
					cont = True
		if cont == True:
			continue
		union_matches = set(pair)
		for pairCopy in reference_copy:
			for elem in union_matches:
				if elem in pairCopy:
					union_matches = union_matches.union(set(pairCopy))
					print(all_unions)
					break
		print(all_unions)
		already_used = False
		for each_set in all_unions:
			if each_set == union_matches:
				already_used = True
				break
		if already_used == False:
			all_unions.append(union_matches)
		print(all_unions)
	return all_unions

def intersection_of_hypotheticals(all_unions, num_species):
	intersection = []
	for each_set in all_unions:
		if len(each_set) == num_species:
			intersection.append(each_set)
	return intersection

def union_of_hypotheticals(all_unions):
	new_union = {}
	for each_set in all_unions:
		new_union = each_set.union(new_union)
	return new_union

def names_for_ecoli(filename, candidate):
	c1 = reader(filename)
	real_names = []
	for elem in c1:
		k = 0
		while k<=len(elem)-1:
			if elem[k:k+8] == 'protein=':
				y = k + 8
				while y<=len(elem) - 1:
					if elem[y] == ']':
						real_names.append(elem[k+8:y])
						break
					y += 1
			k += 1
	f = open('EColiProteinNames.txt', 'w')
	for elem in real_names:
		f.write(elem + '\n')
	f.close()
	return None
