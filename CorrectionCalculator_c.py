import matplotlib.pyplot as plt

def reader(filename):
    f = open(filename, 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()
    return lines

def concat(lines):
	DNA = ''
	for elem in lines:
		DNA += elem
	return DNA

def number_of_G_and_C(string):
	count_g, count_c = 0, 0
	for elem in string:
		if elem == 'G' or elem == 'g':
			count_g += 1
		elif elem == 'C' or elem == 'c':
			count_c += 1
	return count_g, count_c

def min_skew(x_axis, y_values):
	minimum = min(y_values)
	number_of_minimums, index = 0, 0
	while index < len(y_values):
		if y_values[index] == minimum:
			#number_of_minimums += 1
			return index
			#print(min_index, minimum)
		index += 1
	#if number_of_minimums == 1:
	#	return min_index
	#else:
	#	print('There are two minima')

def max_skew(x_axis, y_values):
	maximum = max(y_values)
	number_of_maximums, index = 0, 0
	while index < len(y_values):
		if y_values[index] == maximum:
			number_of_maximums += 1
			# maximum_index = index
			#print(number_of_maximums)
			return index
		index += 1
	#if number_of_maximums == 1:
	#	return maximum_index
	#else:
	#	print('There are two maxima')
	
def windowCalculator(concatenatedDNA, windowSize):
	all_window_values = []
	while len(concatenatedDNA) >= windowSize:
		curr_window = concatenatedDNA[:windowSize]
		G, C = number_of_G_and_C(curr_window)
		curr_skew = (G-C)/(G+C)
		all_window_values.append(curr_skew)
		concatenatedDNA = concatenatedDNA[windowSize:]
	return all_window_values

def cummulative_plotter(all_window_values, windowSize):
	y_values = []
	length_x_axis = len(all_window_values)
	for curr_index in range(length_x_axis):
		curr_value = all_window_values[curr_index]
		if curr_index != 0:
			curr_value += y_values[curr_index-1]
		y_values.append(curr_value)
	x_axis = [elem*windowSize for elem in range(len(y_values))]
	plt.plot(x_axis, y_values)
	plt.show()
	return 'Minimum: ' + str(min_skew(x_axis, y_values)*windowSize), 'Maximum: ' + str(max_skew(x_axis, y_values)*windowSize)
