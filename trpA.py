def stringify(filename):
	x = reader(filename)
	a = ''
	for elem in x:
		if '>' not in elem:
			a += elem
	return a

def reader(filename):
    f = open(filename, 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()
    return lines

def starts(lines):
	allStartingPos = []
	copy = lines[0::2]
	for elem in copy:
		k = 33
		while k <=len(elem)-1:
			if elem[k:k+3] == 'pos':
				i = k+4
				startingPos = ''
				while elem[i] != ' ':
					startingPos += elem[i]
					i += 1
				allStartingPos.append(eval(startingPos))
				break
			k += 1
	return allStartingPos

def constructBackbone(lines, allStartingPos, fileNumber):	
	# First, group all 100% overlapping sequences (overlapping in position)
	k, redundantSequences = 0, {}
	copy = lines[1::2]
	for elem in copy:
		currStartIndex = str(allStartingPos[k])

		if currStartIndex in redundantSequences.keys():
			redundantSequences[currStartIndex].append(elem)
		else:
			redundantSequences[currStartIndex] = [elem]
		k+=1

	# Now start compiling all the letter counts for each position in the genome, and use the letter that occurs the most
	positionDict = {}
	printerCount = 0
	for pair in redundantSequences.items():
		startIndex = pair[0]
		i = 0
		outerLimit = len(pair[1]) - 1
		while i <= outerLimit:
			segment = pair[1][i]
			k = 0
			limit = len(segment) - 1
			while k <= limit:
				letter = segment[k]
				backboneIndex = k + eval(startIndex)
				if str(backboneIndex) not in positionDict.keys():
					# Initialize a dictionary for that position that will hold letter counts
					positionDict[str(backboneIndex)] = {}
				letterCounter = positionDict[str(backboneIndex)]
				if letter not in letterCounter.keys():
					letterCounter[letter] = 1

				else:
					letterCounter[letter] += 1

				k+=1
			i += 1
		printerCount += 1
		if printerCount % 10000 == 0:
			print(len(redundantSequences) - printerCount)
	newFileName = 'PosDictionary' + str(fileNumber*100) + '.txt'
	y = open(newFileName, 'w')
	y.write(str(positionDict))
	y.close()

	# Start inserting letters into the backbone in order
	backboneString = ''
	index = 0
	limit = eval(max(positionDict.keys()))
	while index <= limit:
		try:
			letterCounts = positionDict[(str(index))]
			letterToUse = max(letterCounts.items(), key = lambda x: x[1])[0]
			backboneString += letterToUse
		except KeyError:
			print('Missing: ' + str(index))
		print(index)
		index += 1
	newFileName = 'CleanedFile' + str(fileNumber) + '.txt'

	f = open(newFileName, 'w')
	f.write(backboneString)
	f.close()
	return backboneString

def f():
	files = ['ERS153020', 'ERS153045', 'ERS153046.fasta', 'ERS153047', 'ERS351377', 'ERS351383', 'ERS351384', 'ERS351392', 'ERS747489', 'ERS747490', 'ERS747491', 'ERS747492', 'ERS747493', 'ERS747494', 'ERS153021', 'ERS153022', 'ERS351385']
	k = 2
	for elem in files:
		print('Starting file ' + str(k))

		x = reader(elem)
		print('Finished reading for file ' + str(k))

		allStartingPos = starts(x)
		print('Got the starting positions for file ' + str(k))

		backbone = constructBackbone(x, allStartingPos, k)
		print(str(k) + ' has been completed')
		k += 1

def findTrpA():
	files = ['ERS153019', 'ERS153020', 'ERS153045', 'ERS153046', 'ERS153047', 'ERS351377', 'ERS351383', 'ERS351384', 'ERS351392', 'ERS747489', 'ERS747490', 'ERS747491', 'ERS747492', 'ERS747493', 'ERS747494', 'ERS153021', 'ERS153022', 'ERS351385']
	trpAVersions = []
	k = 0
	reference = stringify('trpA Reference.txt')
	while k<=len(files) - 1:
		currFileName = files[k] + '.txt'
		currGenome = stringify(currFileName)
		# Now search through the genome for trpA similarities
		i = 0
		while i<=len(currGenome):
			if currGenome[i:i+12] == 'ATGAGTAAATTA':
				trpAVersions.append(currGenome[i:i+762])
			i += 1
		k += 1
	k = 0
	f = open('Clinical trpA Versions', 'w')
	for elem in trpAVersions:
		f.write('>' + files[k] + '\n' + elem + '\n')
		k += 1
	f.close()




