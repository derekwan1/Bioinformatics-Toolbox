def exactMatch(string1, string2, string3):
	length1, length2, length3 = len(string1), len(string2), len(string3)
	minimum = min(length1, length2, length3)
	k = 0
	f = open('ExactMatches.txt', 'w')
	while k<=minimum-1:
		print(str(minimum-k) + ' iterations left to go!')
		if string1[k] == string2[k] ==string3[k] != '-':
			sequence_size = MEM(string1[k:], string2[k:], string3[k:])
			print(string1[k:28], string2[k:28], string3[k:28])
			break
			if sequence_size >= 20:
				f.write('>' + str(k) + '-' + str(k+sequence_size) + '\n')
				f.write(string1[k:k+sequence_size] + '\n')
				k += sequence_size - 1
		k += 1
	f.close()
	return None

def newProcessLoops(filename, f):
	lines = reader(filename)
	string1, string2, string3 = f(do_nothing)
	stringList = [string1, string2, string3]
	k = 0
	f = open('LoopsThatAreBackbone.txt', 'w')
	while k<=len(lines) - 1:
		print(str((len(lines)-1) - k))
		loopStartIndex = isolator(lines[k], 'start')
		try:
			loopEndIndex = isolator(lines[k+2], 'end')
		except IndexError:
			loopEndIndex = len(s1) + 1
		percentage = compare(string1[loopStartIndex:loopEndIndex], string2[loopStartIndex:loopEndIndex], string3[loopStartIndex:loopEndIndex])
		if percentage > 0.76:
		# Find a sequence that has minimum number of dashes
			reference = most_letters(['.', string1, '.', string2, '.', string3])
			f.write('>' + str(loopStartIndex) + '-' + str(loopEndIndex) + '\n')
			f.write(reference + '\n')
		k += 2
	f.close()
	return None

def construction(match_filename, loop_filename):
	matches = reader(match_filename)
	almost_matches = reader(loop_filename)
	k = 0
	f = open('Construction.txt', 'w')
	while matches and almost_matches:
		start1 = isolator(matches[k], 'end')
		start2 = isolator(almost_matches[k], 'end')
		if start1 < start2:
			f.write(matches[k+1])
			matches = matches[2:]
		elif start1 > start2:
			f.write(almost_matches[k+1])
			almost_matches = almost_matches[2:]
	if matches and not almost_matches:
		while k<=len(matches) - 1:
			f.write(matches[k+1])
			k += 2
	elif not matches and almost_matches:
		while k<=len(almost_matches) - 1:
			f.write(almost_matches[k+1])
			k += 2
	f.close()
	return None

def compare(string1, string2, string3):
	reference = most_letters(['.', string1, '.', string2, '.', string3])
	stringList = [string1, string2, string3]
	minimum = 1
	k = 0
	while k<=len(stringList) - 1:
		curr_sequence = stringList[k]
		z, count = 0, 0
		while z <= len(curr_sequence) - 1:
			if curr_sequence[z] == reference[z] != '-':
				count += 1
			z += 1
		x = count/len(curr_sequence)
		if x <= 0.76:
			return 0
		elif x < minimum:
			minimum = x
		k += 1
	return minimum

def isolator(string, start_or_end):
	index = 0
	while index <= len(string) - 1:
		if string[index] == '-':
			if start_or_end == 'start':
				string = string[index+1:]
			else:
				string = string[:index]
				string = string[1:]
		index += 1
	return int(string)

def cut(stringList, cutIndex):
	k = 0
	newList = []
	while k<=len(stringList) - 1:
		newList.append(stringList[k][cutIndex:])
		k += 1
	return newList

def check_gap(string_list, k):
	boolean = False
	for elem in string_list:
		if elem[k] == '-' and boolean == False:
			boolean = True
		elif elem[k] == '-' and boolean == True:
			return True
	return False

def gap_size(string_list):
	count = 0
	string1, string2, string3 = string_list[0], string_list[1], string_list[2]
	k = 0
	minimum = min(len(string1), len(string2), len(string3))
	while k<=minimum - 1:
		if string1[k] == string2[k] ==string3[k] != '-':
			return count
		else:
			count += 1
			k += 1

def MEM(string1, string2, string3):
	minimum = min(len(string1), len(string2), len(string3))
	k = 0
	while k<=minimum-1:
		if string1[k] == string2[k] ==string3[k] != '-':
			k += 1
		else:
			break
	return k

def reader(filename):
    f = open(filename, 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()
    return lines

def filter_gap(gap):
	k = 1
	new_gap = []
	while k<=len(gap) - 1:
		if not all_dashes(gap[k]):
			new_gap.append(gap[k-1])
			new_gap.append(gap[k])
		k += 2
	return new_gap

def all_dashes(sequence):
	for elem in sequence:
		if elem != '-':
			return False
	return True

def similar(gap):
	reference = most_letters(gap)
	k = 1
	while k<=len(gap) - 1:
		curr_sequence = gap[k]
		z, count, denom = 0, 0, len(curr_sequence)
		while z <= len(curr_sequence) - 1:
			if reference[z] == curr_sequence[z]:
				count += 1
			z += 1
		ratio = count/denom
		if ratio <= 0.76:
			return False
		k += 2
	return True

def most_letters(gap):
	k, maximum = 1, 0
	while k<=len(gap) - 1:
		count = 0
		curr_sequence = gap[k]
		for elem in curr_sequence:
			if elem != '-':
				count += 1
		if count > maximum:
			maximum = count
			reference = curr_sequence
		k += 2
	return reference

def create_construction(func1, func2, master_f, match_filename, loop_filename, do_nothing):
	func1(f(do_nothing))
	func2(match_filename, f)
	construction(match_filename, loop_filename)
	return None

def concat(lst):
	string = ''
	for elem in lst:
		if '>' not in elem:
			string += elem
	return string

def string_generator(func, filename1, filename2, filename3):
	lines1 = reader(filename1)
	lines2 = reader(filename2)
	lines3 = reader(filename3)
	string1 = concat(lines1)
	string2 = concat(lines2)
	string3 = concat(lines3)
	return func(string1, string2, string3)
