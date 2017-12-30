def reader(filename):
    f = open(filename, 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()
    return lines

def remove_quotes(lines, manual_break):
	k = 0
	new_lines = []
	while k<=len(lines) - 1:
		curr_protein = lines[k]
		if curr_protein != manual_break:
			curr_protein = eval(curr_protein)
		new_lines.append(curr_protein)
		k += 1
	return new_lines

def get_real_name(new_lines, characters_before_number, manual_break):
	k1 = 0
	updated_lines = []
	while k1 <= len(new_lines) - 1:
		curr_protein = new_lines[k1]
		name = curr_protein[0]
		if name != manual_break[0]:
			index = index_of_actual_name(name)
			number = name[characters_before_number:(index-1)]
			number = int(number)
			name = name[index:]
			updated_lines.append([name, number, curr_protein[1:]])
		else:
			updated_lines.append(manual_break)
		k1+=1
	return updated_lines

def index_of_actual_name(string):
	index = 0
	while index <= len(string) - 1:
		if string[index] == ' ':
			return index + 1
		else:
			index += 1

def match_to_common_protein(updated_lines, common_protein_filename, manual_break):
	lines = reader(common_protein_filename)
	common_proteins_with_Chi = []
	other_proteins = []
	before_break = True
	for protein in updated_lines:
		if protein == manual_break:
			before_break = False
		else:
			if before_break == True:
				protein_name = protein[0]
				if protein_name in lines:
					common_proteins_with_Chi.append(protein)
				else:
					other_proteins.append(protein)
			else:
				other_proteins.append(protein)
	return [common_proteins_with_Chi, 'BREAK', other_proteins]

def organize_other_proteins(other_proteins):
	organized = []
	numbers = []
	for protein in other_proteins:
		numbers.append(protein[1])
	numbers = sorted(numbers)
	k = 0
	while k <= len(other_proteins) - 1:
		if other_proteins[k][1] == numbers[0]:
			organized.append(other_proteins[k])
			numbers = numbers[1:]
			hidden_output = other_proteins.pop(k)
			k = -1
		if len(numbers) == 0:
			k = len(other_proteins)
		k += 1
	return organized

def sparse_regions(organized, length_of_candidate):
	k, count = 0, 0
	regions = []
	while k <= len(organized) - 1:
		curr_protein_index = organized[k][1]
		if k == len(organized) - 1:
			count += sum(organized[k][2]) + (length_of_candidate*(len(organized[k][2])-1))
			regions.append([curr_protein_index, count])
			break
		next_protein_index = organized[k+1][1]
		count += sum(organized[k][2]) + (length_of_candidate*(len(organized[k][2])-1))
		if next_protein_index != curr_protein_index + 1:
			regions.append([curr_protein_index, count])
			count = 0
		k += 1
	return regions

def linearize(chi_and_common, sparse, starting_protein, length_of_candidate):
	k = 0
	while k<=len(chi_and_common)-1:
		name = chi_and_common[k][0] 
		if name == starting_protein:
			chi_and_common.extend(chi_and_common[:k])
			chi_and_common[:k] = []
			break
		k += 1
	starting_number = chi_and_common[0][1]
	k = 0
	while k <= len(sparse) - 1:
		number = sparse[k][0]
		if number > starting_number:
			sparse.extend(sparse[:k])
			sparse[:k] = []
			break
		k += 1
	return [chi_and_common, 'BREAK', sparse]

def cumulative_counts(chi_and_common, sparse, length_of_candidate, max_gene_number):
	combined_with_counts = []
	while chi_and_common and sparse:
		index1 = chi_and_common[0][1]
		index2 = sparse[0][0]
		if index2 == max_gene_number:
			curr_protein = chi_and_common[0]
			sequence = curr_protein[2]
			len_sequence = ((len(sequence)-1)*length_of_candidate) + sum(sequence)
			if len(combined_with_counts) != 0:
				cumulative_count = len_sequence + combined_with_counts[-1][-1]
			else:
				cumulative_count = len_sequence
			chi_and_common = chi_and_common[1:]
			curr_protein.append(cumulative_count)
			combined_with_counts.append(curr_protein)
			# Need to put this code in here because there's no other way the function can handle this base case
			curr_protein = sparse[0]
			count = curr_protein[1]
			cumulative_count = combined_with_counts[-1][-1] + count
			sparse = sparse[1:]
			curr_protein.append(cumulative_count)
			combined_with_counts.append(curr_protein)
			continue
		elif index1 == max_gene_number:
			curr_protein = sparse[0]
			count = curr_protein[1]
			cumulative_count = combined_with_counts[-1][-1] + count
			sparse = sparse[1:]
			curr_protein.append(cumulative_count)
			combined_with_counts.append(curr_protein)
			curr_protein = chi_and_common[0]
			sequence = curr_protein[2]
			len_sequence = ((len(sequence)-1)*length_of_candidate) + sum(sequence)
			if len(combined_with_counts) != 0:
				cumulative_count = len_sequence + combined_with_counts[-1][-1]
			else:
				cumulative_count = len_sequence
			chi_and_common = chi_and_common[1:]
			curr_protein.append(cumulative_count)
			combined_with_counts.append(curr_protein)
			continue
		elif index2 < index1:
			curr_protein = sparse[0]
			count = curr_protein[1]
			cumulative_count = combined_with_counts[-1][-1] + count
			sparse = sparse[1:]
		else:
			curr_protein = chi_and_common[0]
			sequence = curr_protein[2]
			len_sequence = ((len(sequence)-1)*length_of_candidate) + sum(sequence)
			if len(combined_with_counts) != 0:
				cumulative_count = len_sequence + combined_with_counts[-1][-1]
			else:
				cumulative_count = len_sequence
			chi_and_common = chi_and_common[1:]
		curr_protein.append(cumulative_count)
		combined_with_counts.append(curr_protein)
	return combined_with_counts

def page_and_len(cumulative_count):
	length = (cumulative_count + 3087)% 26700
	pages = (cumulative_count+3087) // 26700
	return ['Pages before the line: ' + str(pages), 'Length of new page before the line: ' + str(length/1000) + ' cm']


