
def reader(filename):
    f = open(filename, 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()
    return lines

def find_index_of_names(lines):
	k = 0
	indices = []
	while k<=len(lines) - 1:
		if '>' in lines[k]:
			indices.append(k)
		if len(indices) == 2:
			break
		k += 1
	return indices

def separator(lines, indices):
	backbone = ''
	genome = ''
	backbone_index = indices.pop(0)
	genome_index = indices.pop(0)
	for elem in lines[backbone_index+1:genome_index]:
		backbone += elem
	for elem in lines[genome_index+1:]:
		genome += elem
	return [backbone, genome]

def loop_combiner(backbone, genome):
	loop = ''
	k = 0
	while k<=len(backbone)-1:
		if backbone[k] == '-' and genome[k]!= 'n': 
			loop+=genome[k]
		k += 1
	return loop

