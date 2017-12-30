def writer(filename):
	x = reader(filename)
	string = ''
	for elem in x:
		if '>' not in elem:
			string += elem
	f = open('Formatted.fasta', 'w')
	while len(string) >= 1:
		if len(string) >= 60:
			f.write(string[:60] + '\n')
		else:
			f.write(string + '\n')
		string = string[60:]
	f.close()
	return None
	
def reader(filename):
    f = open(filename, 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()
    return lines

def make_leading_strand(filename, origin, terminus, keep_graph_order):
	x = reader(filename)
	string = ''
	for elem in x:
		if '>' not in elem:
			string += elem
	f = open('LeadingStrand.txt', 'w')
	if keep_graph_order == False:
		if origin > terminus:
			front = string[origin:]
			front += string[:terminus]
			back = string[terminus:origin]
		else:
			front = string[origin:terminus]
			back = string[terminus:]
			back += string[:origin]
		RC(back)
		f.write(front)
		string = reader_to_string('RC.txt')
		f.write(string)
	else:
		if origin < terminus:
			first = string[:origin]
			origin_to_terminus = string[origin:terminus]
			second = string[terminus:] 
			RC(origin_to_terminus)
			f.write(first)
			string = reader_to_string('RC.txt')
			f.write(string)
			f.write(second)
		else:
			terminus_to_origin = string[terminus:origin]
			first = string[:terminus]
			second = string[origin:]
			RC(first)
			string = reader_to_string('RC.txt')
			f.write(string)
			f.write(terminus_to_origin)
			RC(second)
			string = reader_to_string('RC.txt')
			f.write(second)
	f.close()
	return None

def reader_to_string(filename):
	x = reader(filename)
	string = ''
	for elem in x:
		if '>' not in elem:
			string += elem
	return string

def RC(string):
	k = len(string) - 1
	g = open('RC.txt', 'w')
	while k >= 0:
		if string[k] == 'A' or string[k] == 'a':
			g.write('T')
		elif string[k] == 'T' or string[k] == 't':
			g.write('A')
		elif string[k] == 'G' or string[k] == 'g':
			g.write('C')
		elif string[k] == 'C' or string[k] == 'c':
			g.write('G')
		k -= 1
	g.close()
	return None