from math import log

"HELPER FUNCTIONS BELOW THIS LINE"
"_____________________________________________________________________________________________"

	
def reader(filename):
    f = open(filename, 'r')
    lines = [line.strip() for line in f.readlines() if '>' not in line]
    f.close()
    return lines

def GC_count(sequence):
	count_c, count_g = 0, 0
	for letter in sequence:
		if letter == 'G' or letter == 'g':
			count_g += 1
		if letter == 'C' or letter == 'c':
			count_c += 1
	return count_g, count_c


def consecutive_repeats(sequence, max_repeat_length=3):
	"Returns True if there are consecutive repeats of max_repeat_length"
	k = 0
	while k<=len(sequence)-max_repeat_length:
		curr_window = sequence[k:k+max_repeat_length]
		first_letter = curr_window[0]
		k2 = 0
		while k2 <= len(curr_window) - 1:
			curr_letter = curr_window[k2]
			if curr_letter != first_letter:
				break
			if k2 == len(curr_window) - 1 and curr_letter == first_letter:
				return True
			k2 += 1
		k += 1
	return False

def tandem_repeats(sequence, tandem_length = 2, number_of_repeats = 2):
	"Returns True if primer contains number_of_repeats tandem repeats of length tandem_length"
	k = 0
	while k<=len(sequence) - tandem_length:
		curr_tandem = sequence[k:k+tandem_length]
		rest_of_seq = sequence[k+tandem_length:]
		k2, num_repeats = 0, 1
		while k2 <= len(rest_of_seq) - tandem_length:
			if num_repeats >= number_of_repeats:
				return True
			if rest_of_seq[k2:k2+tandem_length] != curr_tandem:
				break
			num_repeats += 1
			k2 += tandem_length
		k += 1
	return False

def GC_beginning(primer, num_GC):
	window = primer[:num_GC]
	g, c = GC_count(window)
	if g + c == num_GC:
		return True

def RC(string):
	k = len(string) - 1
	RC_string = ''
	while k >= 0:
		letter = string[k]
		if letter == 'A' or letter == 'a':
			RC_string += 'T'
		elif letter == 'T' or letter =='t':
			RC_string += 'A'
		elif letter == 'G' or letter == 'g':
			RC_string += 'C'
		elif letter =='C' or letter =='c':
			RC_string += 'G'
		k -= 1
	return RC_string

def check_complementary_bp(nuc1, nuc2):
	if nuc1 == 'A':
		return nuc2 == 'T'
	elif nuc1 == 'T':
		return nuc2 == 'A'
	elif nuc1 == 'G':
		return nuc2 == 'C'
	elif nuc1 == 'C':
		return nuc2 == 'G'

def total_pairings(string1, string2):
	"""Finds the number of pairings between two strings"""
	min_length = min(len(string1), len(string2))
	k, count = 0, 0
	while k<=min_length-1:
		if check_complementary_bp(string1[k], string2[k]):
			count += 1
		k += 1
	return count

def consecutive_pairings(string1, string2):
	"""Finds the maximum number of consecutive pairings between two strings"""
	min_length = min(len(string1), len(string2))
	k, count, consecutive_runs = 0, 0, []
	while k<=min_length-1:
		while k<=min_length-1 and check_complementary_bp(string1[k], string2[k]):
			count += 1
			k += 1
		consecutive_runs.append(count)
		if count == 0:
			k += 1
		count = 0
	return max(consecutive_runs)

def checker(filename, primer):
	"Checks whether the primer only occurs once in the sequence"
	x = reader(filename)
	string = ''
	for elem in x:
		string += elem
	k, count = 0, 0
	while k<=len(string) - 1:
		if string[k:k+len(primer)] == primer:
			count += 1
		k += 1
	return "OK" if count == 1 else "Not unique"

def max_consecutive_pairs(total_pairs, consecutive_pairs):
	"""Returns the maximum consecutive pairs--if there are multiple positions that give the same maximum consecutive pairs, this returns
	the index that gives the largest number of total pairs with the position associated with the maximum consecutive pairs index
	"""
	maximum = max(consecutive_pairs)
	k, indices = 0, []
	while k<=len(consecutive_pairs) - 1:
		if consecutive_pairs[k] == maximum:
			indices.append(k)
		k += 1
	max_total_pairs, max_index = 0, 0
	for elem in indices:
		if total_pairs[elem] > max_total_pairs:
			max_total_pairs = total_pairs[elem]
			max_index = elem
	return 'Shift: ' + str(max_index), max_total_pairs, maximum #Shift, maximum total pairs associated with the position that gives the max consecutive pairs, and maximum consecutive pairs

def amplicon_size(primer, fixed_primer, string):
	if primer in string:
		# This means that primer is the forward primer and fixed_primer is the reversed one
		fixed_primer = RC(fixed_primer)
	else:
		primer = RC(primer)
	k = 0
	while k<=len(string) - 1:
		if string[k:k+len(primer)] == primer:
			primer_start_index = k
		if string[k:k+len(fixed_primer)] == fixed_primer:
			fixed_primer_start_index = k
		k += 1
	# Determine which primer is the farther one
	if primer_start_index > fixed_primer_start_index:
		return primer_start_index - (fixed_primer_start_index + len(fixed_primer))
	else:
		return fixed_primer_start_index - (primer_start_index + len(primer))

"CORE FUNCTIONS BELOW THIS LINE"
"_____________________________________________________________________________________________"


def prelim_primers(sequence, primer_size, min_GC, max_GC):
	"""sequence is the string representation of the sequence in which we want to make a primer
	   primer_size is the length of the desired primer
	   min_GC and max_GC are the lower and upper bound proportion of GC allowed in the primer
	   num_G_C_in_beginning is the number of consecutive GC desired in the beginning of the primer
	   num_G_C_in_end is the number of GC's (don't need to be consecutive) desired in the last 5 nucleotides of the primer

	   As a whole, this function chooses preliminary primers based on GC proportion, size, and GC's in the end
	   May need to pass the possible primers through tandem_repeats and consecutive_repeats
	"""
	k, possible_primers = 0, []
	while k <= len(sequence) - primer_size:
		curr_primer = sequence[k:k+primer_size]
		g, c = GC_count(curr_primer)
		GC_content = (g+c) / primer_size
		if min_GC <= GC_content <= max_GC:
			last_five = curr_primer[primer_size-5:]
			g, c = GC_count(last_five)
			if (g + c) <= 3:
				possible_primers.append(curr_primer)
		k += 1
	return possible_primers

def shift_right(stationary, sliding):
	"""Stationary sequence is stationary; sliding sequence slides to the right relative to stationary
		Returns the number of pairings for each position
	"""
	total_pairs, consecutive_pairs = [], []
	# Shift index is the amount by which the sliding primer is shifted to the right
	# Sliding index is the index of the sliding strand
	shift_index, sliding_index = 0, 0
	while shift_index <= len(stationary):
		# We just add nonsense semicolons to the beginning of the string because total_pairings and 
		# consecutive_pairings only check bp up to the shorter of the two strings of comparison 
		temp_sliding = (';' * shift_index) + sliding
		total_pairs.append(total_pairings(stationary, temp_sliding))
		consecutive_pairs.append(consecutive_pairings(stationary, temp_sliding))
		shift_index += 1
	return total_pairs, consecutive_pairs


def melting_temp_calculator(primer, salt_concentration, mg_concentration, primer_concentration):
	"""Both equations assume that the annealing occurs under the standard conditions of 50 nM primer, 50 mM Na+, and pH 7.0.
		salt_concentration is in MILLI Molar
		mg_concentration is in MILLI Molar
		primer_concentration is in NANO Molar
		Online programs use defaults of ~ 50 mM salt, 0 mM mg, and 200 nM primer

	Sources: https://www.biophp.org/minitools/melting_temperature/demo.php?formula=basic
			http://www.ncbi.nlm.nih.gov/pmc/articles/PMC19045/table/T2/
	"""
	g, c = GC_count(primer)
	a_and_t = len(primer) - g - c
	h, s = 0, 0
	enthalpy_values = {'AA':-7.9, 'AC': -8.4, 'AG':-7.8, 'AT':-7.2, 'CA':-8.5, 'CC':-8, 'CG':-10.6, 'CT': -7.8,
						'GA': -8.2, 'GC': -9.8, 'GG': -8, 'GT': -8.4, 'TA':-7.2, 'TC':-8.2, 'TG':-8.5, 'TT':-7.9}
	entropy_values = {'AA': -22.2, 'AC':-22.4, 'AG':-21, 'AT':-20.4, 'CA':-22.7, 'CC': -19.9, 'CG':-27.2, 'CT':-21,
						'GA': -22.2, 'GC': -24.4, 'GG': -19.9, 'GT': -22.4, 'TA': -21.3, 'TC': -22.2, 'TG': -22.7, 'TT': -22.2
	}

	"Take into account the effect of Mg and salts"
	salt_effect = (salt_concentration / 1000) + ((mg_concentration/1000) * 140)
	s += 0.368 * ((len(primer)-1)*log(salt_effect))

	"Adjustments for the terminal nucleotides"
	first_nucleotide = primer[0]
	if first_nucleotide == 'G' or first_nucleotide == 'C':
		h += 0.1
		s += -2.8
	elif first_nucleotide == 'A' or first_nucleotide == 'T':
		h += 2.3
		s+=4.1
	last_nucleotide = primer[len(primer)-1]
	if last_nucleotide == 'G' or last_nucleotide == 'C':
		h += 0.1
		s += -2.8
	elif last_nucleotide == 'A' or last_nucleotide == 'T':
		h+=2.3
		s+=4.1

	"Take into account the overall nucleotide composition of the primer sequence"
	k = 0
	while k<=len(primer) - 2:
		window = primer[k:k+2]
		h += enthalpy_values[window]
		s += entropy_values[window]
		k += 1
	return ((1000*h)/(s+(1.987*log(primer_concentration/2000000000))))-273.15

def check_dimer(forward, reverse, num_consecutive, num_total):
	"""Input both primers 5->3
	num_consecutive is the number of consecutive bp pairs below which a dimer is considered unlikely
	num_total is the number of total bp pairs below which a dimer is considered unlikely 
	We determine dimerization based on the position of greatest alignment

	FALSE return value means that there is unlikely dimerization
	TRUE return value means that there is likely dimerization
	"""
	reverse_complement = ''
	for elem in reverse:
		reverse_complement = elem + reverse_complement
	"Calculate the longest consecutive stretch of complementary base pairs as well as the total number of complementary base pairs"

	reverse_shift = shift_right(forward, reverse_complement)
	forward_shift = shift_right(reverse_complement, forward)

	# Unpack the tuple
	total_pairs1, consecutive_pairs1 = reverse_shift[0], reverse_shift[1]
	total_pairs2, consecutive_pairs2 = forward_shift[0], forward_shift[1]

	# Grab the maximum pairing position
	reverse_shift = max_consecutive_pairs(total_pairs1, consecutive_pairs1)
	forward_shift = max_consecutive_pairs(total_pairs2, consecutive_pairs2)

	shift1, total_pairs1, consecutive_pairs1 = 'Reverse ' + reverse_shift[0] + ' to the right relative to the forward', reverse_shift[1], reverse_shift[2]
	shift2, total_pairs2, consecutive_pairs2 = 'Forward ' + forward_shift[0] + ' to the right relative to the reverse', forward_shift[1], forward_shift[2]

	if consecutive_pairs1 > consecutive_pairs2:
		return total_pairs1 >= num_total or consecutive_pairs1 >= num_consecutive
	elif consecutive_pairs1 < consecutive_pairs2:
		return total_pairs2 >= num_total or consecutive_pairs2 >= num_consecutive
	else:
		# In the event that both consecutive pairs are equal, use total pairs
		return max(total_pairs1, total_pairs2) >= num_total or consecutive_pairs1 >= num_consecutive

	"Use consecutive pairings as the main criteria for determining the most likely dimer: if they are equal, look at total pairings"
	# if consecutive_pairs1 > consecutive_pairs2:
	# 	return shift1, 'Total Pairings: ' + str(total_pairs1), 'Consecutive Pairings: ' + str(consecutive_pairs1)
	# elif consecutive_pairs1 < consecutive_pairs2:
	# 	return shift2, 'Total Pairings: ' + str(total_pairs2), 'Consecutive Pairings: ' + str(consecutive_pairs2)
	# else:
	# 	if total_pairs1 > total_pairs2:
	# 		return shift1, 'Total Pairings: ' + str(total_pairs1), 'Consecutive Pairings: ' + str(consecutive_pairs1)
	# 	else:
	# 		return shift2, 'Total Pairings: ' + str(total_pairs2), 'Consecutive Pairings: ' + str(consecutive_pairs2)

"PIPELINE FUNCTIONS BELOW THIS LINE"
"_____________________________________________________________________________________________"


def pipeline(filename, primer_size, min_GC, max_GC, num_G_C_in_beginning, tandem_length, number_of_tandem_repeats, max_repeat_length, reverse_complement, salt_concentration, mg_concentration, primer_concentration, fixed_primer, delta_T, num_consecutive, num_total):
	"Returns all possible forward-reverse primer pairs that satisfy Diassess parameters"
	lines = reader(filename)
	string = ''
	for elem in lines:
		string += elem

	Tm1 = melting_temp_calculator(fixed_primer, salt_concentration, mg_concentration, primer_concentration)

	if reverse_complement:
		"""If we passed in the reverse complement, this means that the fixed primer is forward
		   1. Grab the window of 80-150 bp, reverse complement it
		   2. Search for preliminary primers based on overall GC composition, terminal GC's, and size (prelim_primers)
		   3. Compare the melting temperatures
		   4. Check for dimerization
		"""
		k = 0
		while k<=len(string) - 1:
			if string[k:k+primer_size]==fixed_primer:
				window = string[k+80+primer_size:k+151+primer_size]
				window = RC(window)
				break
			k += 1
	else:
		"If we passed in False for reverse_complement, this means that the fixed primer is reverse"
		fixed_primer_copy = RC(fixed_primer)
		k = 0
		while k<=len(string) - 1:
			if string[k:k+len(fixed_primer_copy)] == fixed_primer_copy:
				window = string[k-primer_size-151:k-primer_size-80]
				break
			k += 1

	primers = prelim_primers(window, primer_size, min_GC, max_GC)
	primers = [primer for primer in primers if abs(Tm1 - melting_temp_calculator(primer, salt_concentration, mg_concentration, primer_concentration)) <= delta_T]
	"Enable these dimerization programs if it's too tedious to check manually online"
	if reverse_complement:
		primers = [primer for primer in primers if not check_dimer(fixed_primer, primer, num_consecutive, num_total)]
	else:
		primers = [primer for primer in primers if not check_dimer(primer, fixed_primer, num_consecutive, num_total)]

	#"Optional filtering for beginning GC's, if further filtering is needed"
	#primers = [primer for primer in primers if GC_beginning(primer, num_G_C_in_beginning)]
	#print(len(primers))

	#"Optional filtering for tandem repeats and runs"
	#primers = [primer for primer in primers if not consecutive_repeats(primer, max_repeat_length)]
	#primers = [primer for primer in primers if not tandem_repeats(primer, tandem_length, number_of_tandem_repeats)]
	"These primers are 5 -> 3"
	return [[primer, 'amplicon size: ' + str(amplicon_size(primer, fixed_primer, string))] for primer in primers]

def f():
	"Custom parameters"
	filename = 'parCS87RIW.muscle.fasta'
	num_G_C_in_beginning = 2

	# If fixed primer is reverse, it should be passed in as 5 -> 3
	fixed_primer = 'CGCCTCATAGGCGGAACT'.upper()
	primer_size = len(fixed_primer)

	"50/50 fixed and custom"
	"Set reverse_complement True if fixed_primer is the forward"
	reverse_complement = False

	"Generally fixed parameters"
	num_consecutive = 4 # 4 or more consecutive pairings is considered a dimer
	num_total = 8 # 7 or more total pairings is considered a dimer
	# Note: the above parameters are related in an 'or' fashion. 

	min_GC = 0.5
	max_GC = 0.6

	tandem_length = 2
	number_of_tandem_repeats = 3
	max_repeat_length = 3

	salt_concentration = 50
	mg_concentration = 0
	primer_concentration = 200

	delta_T = 3

	return pipeline(filename, primer_size, min_GC, max_GC, num_G_C_in_beginning, tandem_length, number_of_tandem_repeats, max_repeat_length, reverse_complement, salt_concentration, mg_concentration, primer_concentration, fixed_primer, delta_T, num_consecutive, num_total)


