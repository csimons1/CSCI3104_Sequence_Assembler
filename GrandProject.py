#==================================================================
#
#					CSCI3104 Grand Project 1
#					De-Novo Genome Assembler
#	
#							Group:
#				Christian Simons, Connor Anthony 
#				 Melissa Pollich, Johnson Cheung
#
#==================================================================


# Imports and Global Constants

import numpy as np # For Tie Breaking

import sys	# For CLIs
import time as time	# Debug/Performance Monitoring
import os	# Debug/Performance Monitoring
clear = lambda: os.system('cls')	# Debug/Performance Monitoring

discardSize = 0 # Set to be equal to 20% of the initial input length.


# Helper Functions
# Christian Simons - 04/09/2019

def metricSetter(n):
# Takes in the length of one of the input reads and
# sets the global discardSize equal to 20% of that length.
# This metric is arbitrary.
#============================================================
	global discardSize
	discardSize = int(n * 0.20)
#============================================================

def sequencePadder(A, B):
# Appends underscores until string sizes are equal in case of 
# reads of different sizes (this will happen).
# Underscores get scrubbed before returning.
#============================================================
	if len(A) < len(B):
		lenDelta = len(B) - len(A)
		for e in range(lenDelta):
			A = A + '_'
	elif len(B) < len(A):
		lenDelta = len(A) - len(B)
		for e in range(lenDelta):
			B = B + '_'
	
	return A,B
#============================================================


def paddingScrubber(sequence):
# Removes the padded underscores appended
# to the end of a read
#============================================================
	scrubDistanceStart = 0
	for i in range(len(sequence)):
		if sequence[i] != '_':
			scrubDistanceStart += 1
			
	if scrubDistanceStart > 0:
		A = sequence[:scrubDistanceStart]

	return A
#============================================================

# Contig Assembler Function
# Christian Simons - 03/11/2019
# Refactored - 04/11/2019

def contigAssembler(A, B, direction):
	# Direction controls the string that gets offset. 'r' -> prepend
    # '-' to A. 'l' -> prepend '-' to B.
    # Code differences between r and l is just a matter of switching
    # the arrays in the code
	# =================================================================
	if direction == 'r':
	# =================================================================
		
		
	# Initialization
	# offset keeps track of how many prepended dummy characters
	# are in the string, in order to simplfy the comparison method.
	# A visual example of a overlapping aligned contigs:
	# A: -----------AATCCCCGTACATGTTGTTA
	# B: GCGCCTTGGCTAATCCCCGT
	#               ^^ Point of Alignment
	# contigPossible flags whether an alignment meets the criteria
	# of a contig. The outside of loop initialization makes no
	# functional sense, but appears to be required in order for
	# execution.
	# =================================================================
		offset = 0
		contigPossible = False
	# =================================================================
	
	
	# Due to a lack of understanding of lambda and toying with
	# python's Y-Combinator, a while loop is used in order to speed
	# up processing by effectively simulating Tail Call Optimization.
	# =================================================================
		while not (contigPossible):
	# =================================================================
		
		
	# Assume all alignments are correct, anytime a non-aligned 
	# nucleotide pair is compared, the flag is set to false.
	# If the flag is false, then the loop (the "tail recursion")
	# continues
	# =================================================================		
			contigPossible = True
	# =================================================================
	
	
	# Base Case for "recursion", if the prepended dashes exceed the 
	# length of the other string, then an alignment is impossible.
	# =================================================================		
			if offset >= len(B):
				return (None, None)
	# =================================================================	
	
	
	# overlapCounter keeps track of how many contiguous overlapping
	# nucleotides exist so far. Only equal index positions need to
	# be checked, due to the offset. Underscores are ignored, as they
	# are appended to strings that are longer than the other string.
	# This prevents index bounding errors. A break is issued when
	# an underscore is reached, because the underscore indicates the
	# end of relevant information in the string, no nucleotides will
	# follow an underscore.
	# =================================================================	
			overlapCounter = 0
			for j in range(offset, len(B)):
				if (A[j] != B[j]) and (A[j] != '_'):
					contigPossible = False
					break
				elif (A[j] == B[j]) and (A[j] != '_'):
					overlapCounter += 1
				if A[j] == '_' or B[j] == '_':
					break
			
			if overlapCounter <= discardSize:
				contigPossible = False
	# =================================================================	
	
	
	# Checks for residual underscore padding, and removes it.
	# The offsets are used for slice points in order to concatenate
	# without doubling the overlap region.
	# =================================================================	
			if contigPossible == True:
				if A[-1] == '_':
					A = paddingScrubber(A)
				elif B[-1] == '_':
					B = paddingScrubber(B)
				return (B[:offset] + A[offset:], overlapCounter)
	# =================================================================	
	
	
	# Parameter adjustment to effectively create tail recursion.
	# =================================================================	
			else:
				A = '-' + A
				offset += 1
	# =================================================================	

	
	elif direction == 'l':
		offset = 0
		contigPossible = False
		
		while not (contigPossible):
			contigPossible = True
			if offset >= len(A):
				return (None, None)
		
			overlapCounter = 0
			for j in range(offset, len(A)):
				if (A[j] != B[j]) and (A[j] != '_'):
					contigPossible = False
					break
				elif (A[j] == B[j]) and (A[j] != '_'):
					overlapCounter += 1
				if A[j] == '_' or B[j] == '_':
					break
				
			if overlapCounter <= discardSize:
				contigPossible = False
			
			if contigPossible == True:
				if A[-1] == '_':
					A = paddingScrubber(A)
				elif B[-1] == '_':
					B = paddingScrubber(B)	
				return (A[:offset] + B[offset:], overlapCounter)
			else:
				B = '-' + B
				offset += 1
			


# Assembler Controller Function

# Christian Simons - 03/25/2019
# assemblerController runs through four stages of processing:
# 1. Assemble Contigs
# 2. Ensure Unique Read Pairs + Prune Suboptimal Contigs
# 3. Recombine Contigs with Reads
# 04/09/2019 - Basic Testing Done, Implementing Base Case

def assemblerController(sequences):
	reads = sequences.copy()
	results = []

	# Stage 1: Assemble Contigs
	# The First step iterates over all reads that are provided from the <file IO> function.
	# As it is a nested for loop, the first read gets compared to each other read after it in sequence.
	# the constructContig function is called at each comparison, performing both a rightward and leftward
	# alignment. Whichever alignment provides the highest overlap value is saved into the results list.
	# While this processing is performed, only the largest overlaps per i -> j pair is saved.
	#=============================================================================================
	
	# Iterate over reads 0,...,n-1. Read n is excluded
	# because it's been compared to all previous reads
	#=============================================================================================
	for i in range(len(reads) - 1):
	#=============================================================================================
		
		
	# Initialize variables
	#=============================================================================================
		overlapMax = float("-inf")
		temp = 0
		maxContig = None
	#=============================================================================================	
	
	
	# Compares read i to reads i + 1,...,n by iterating over reads i + 1,...,n
	#=============================================================================================
		for j in range(i + 1, len(reads)):
	#=============================================================================================		

	
	# Take reads i and j, and pad them such that their lengths are equal
	# Then attempt to construct a contig between the two by running and leftward and 
	# rightward alignment. Store these results.
	#=============================================================================================			
			read0, read1 = sequencePadder(reads[i], reads[j])
			contigResultsR = contigAssembler(read0, read1, 'r')
			contigResultsL = contigAssembler(read0, read1, 'l')
	#=============================================================================================
	
	
	# Determine the maximum overlap between the two contigs. If the overlap of both
	# alignments are equal, break the tie at random uniformly.
	#=============================================================================================
			if (contigResultsR[1] != None) and (contigResultsL[1] != None):
				if contigResultsR[1] > contigResultsL[1]:
					temp = (contigResultsR, i, j)
				elif contigResultsL[1] > contigResultsR[1]:
					temp = (contigResultsL, i, j)
				elif contigResultsR[1] == contigResultsL[1]:
					#Randomly pick?
					x = np.random.choice([1,2])
					if x == 1:
						temp = (contigResultsR, i, j)
					else:
						temp = (contigResultsL, i, j)
			elif (contigResultsR[1] != None) and (contigResultsL[1] == None):
				temp = (contigResultsR, i, j)
			elif (contigResultsL[1] != None) and (contigResultsR[1] == None):
				temp = (contigResultsL, i, j)
			else:
				continue
			if temp[0][1] > overlapMax:
				overlapMax = temp[0][1]
				maxContig = temp
	#=============================================================================================

	
	# Store the contig with maximum overlap in the results array.
	#=============================================================================================	
		if (maxContig != None):
			results.append(maxContig)
	#=============================================================================================	
		
		
	# Stage 2.0: Ensure Unique Read Pairs
	# The second step uses a nested for loop to iterate over the assembled contigs, if there are a pair
	# of contigs that share the same read, the contig with the larger overlap is stored in a new list,
	# and the values of the contig with the smaller overlap has its values set in a way that they are
	# flagged for deletion.
	#
	# Stage 2.1: Prune Suboptimal Contigs
	# Uses a nested for loop that searches for contigs with the None/-1 values set in the
	# previous step and removes them from the list, creating a new list with only contigs of 
	# the largest large overlap values and unique read pairs.
	#=============================================================================================

	finalContigList = []
	if len(results) > 1:
		for s1 in range(len(results)):
			for s2 in range(s1 + 1, len(results)):
				if (results[s1][1] == results[s2][2]) or (results[s1][2] == results[s2][1]) or (results[s1][1] == results[s2][1]) or (results[s1][2] == results[s2][2]):
					if results[s2][0][0] == None:
						continue
					elif (results[s1][0][1] > results[s2][0][1]):
						results[s2] = ((None, -1), -1, -1)
					elif (results[s2][0][1] > results[s1][0][1]):
						results[s1] = ((None, -1), -1, -1)
					elif (results[s1][0][1] == results[s2][0][1]):
						breakTie = np.random.choice([1,2])
						if breakTie == 1:
							results[s2] = ((None, -1), -1, -1)
						else:
							results[s1] = ((None, -1), -1, -1)

	elif len(results) == 1:
		finalContigList.append(results[0])
	
	# Iterate over the resulting contigs, if they were not flagged for deletion, append them to
	# a new list.
	#=============================================================================================
	if len(results) > 1:
		for seq in range(len(results)):
			if results[seq][0][1] != -1:
				finalContigList.append(results[seq])
	
	#=============================================================================================	
	
	
	# Stage 3: Recombine
	# The fourth step takes the list created from step 3, and recombines them with the original read list
	# in a way that creates a new list that can be ran through the assemblerController again.
	#=============================================================================================		
	if len(finalContigList) > 0:
		readsPresent = []
		for i in range(len(finalContigList)):
			readsPresent.append(finalContigList[i][1])
			readsPresent.append(finalContigList[i][2])
		
		readIndices = set(readsPresent)
		readIndices = list(readIndices)
		finalReads = []
	
		for j in range(len(reads)):
			if j not in readIndices:
				finalReads.append(reads[j])
			
		for k in range(len(finalContigList)):
			finalReads.append(finalContigList[k][0][0])
	else:
		finalReads = reads
	#=============================================================================================
	
	
	# Stage 4: Recursion Control:
	# If the finalReads are the same as what was passed in (reads), then no processing could be
	# performed, indicating that assembly is finished.
	# If they are different, more processing needs to be done, and the function is called
	# recursively.
	# This can probably be refactored into a while loop to act as a tail recursive function.
	#=============================================================================================
	if reads == finalReads:
		return finalReads
	else:
		return assemblerController(finalReads)
	#=============================================================================================
	
	
# N50 Calculator Function
# Christian Simons - 04/09/2019

def determineN50(sequences):
	# Initialization
	# =================================================================	
	N50 = ''
	# =================================================================	
	
	
	# Sort the sequences by length in descending order, if there
	# is more than one contig in the final results list.
	# =================================================================	
	if len(sequences) > 1:
		sequences.sort(key = len)
		sequences.reverse()
	# =================================================================	
	
	
	# Create a list of the lengths of the contigs. The index i 
	# represents which contig, and the value assigned to it
	# is the length.
	# =================================================================	
		sequenceLengths = []
		for i in range(len(sequences)):
			sequenceLengths.append(len(sequences[i]))
	# =================================================================	
	
	
	# n50Length is the index position where the N50 needs to stop
	# being made.
	# Initialization
	# =================================================================	
		n50Length = sum(sequenceLengths)
		n50Length = n50Length // 2
		genomeCutOff = 0
		runningN50Length = 0
	# =================================================================	
	
	
	# Iterate over the lengths list, and add each lengths value
	# to the running total, j holds the value of the last contig
	# to be in the N50.
	# =================================================================	
		for j in range(len(sequenceLengths)):
			runningN50Length += sequenceLengths[j]
			if n50Length < runningN50Length:
				genomeCutOff = j
	# =================================================================	
	
	
	# Concatenates the sequences to create the N50
	# =================================================================	
		for k in range(genomeCutOff):
			N50 = N50 + sequences[k]
	# =================================================================	
	
	
	# If there is only one contig in the final results, then that
	# contig has to be the N50.
	# =================================================================	
	else:
		N50 = sequences[0]
	# =================================================================	
	
	
	return N50
	

# File IO
# Christian Simons - 04/09/2019
# File Format:
# 0. Header
# 1. ATCG ... CCGT
# 2. +
# 3. IIII ... IIII
def fileIO(filename):
	file = open(filename, "r")
	reads = []
	
	# Offset the zero index such that the first line of the file is
	# stored in index 1. Then iterate over all lines in the file and
	# store them temporarily in order to parse through and store
	# only the genome reads.
	# =================================================================	
	parseBlock = ["_"]
	for line in file:
		parseBlock.append(line)
	# =================================================================	
	
	
	# Stores every fourth line, starting at the second line, removes \n.
	# This stores every genome without the extra information.
	# =================================================================	
	for pos in range(2, len(parseBlock), 4):
		reads.append(parseBlock[pos].strip())
	# =================================================================	
	
	file.close()
	
	metricSetter(len(reads[0]))
	
	return reads

#jupyter notebook --NotebookApp.iopub_data_rate_limit=1.0e10 # Needed to use fileIO() in jupyter.

# Main Function
# Christian Simons - 04/09/2019
def main():
	return determineN50(assemblerController(fileIO(sys.argv[1])))


	
# Prints and Main Call
# 'C:/Users/Christian/Downloads/rand.500.1.fq'
clear()	# Debug/Performance Monitoring
timeObject_st = time.localtime(None)	# Debug/Performance Monitoring
hour = timeObject_st[3]	# Debug/Performance Monitoring
minute = timeObject_st[4]	# Debug/Performance Monitoring
print("Running...")	# Debug/Performance Monitoring
print("Start Time: " + str(hour) + ':' + str(minute))	# Debug/Performance Monitoring

N50 = main()	# What we've all been waiting for.

timeObject_fn = time.localtime(None)	# Debug/Performance Monitoring
hour_fn = timeObject_fn[3]	# Debug/Performance Monitoring
minute_fn = timeObject_fn[4]	# Debug/Performance Monitoring
print('')	# Debug/Performance Monitoring
print("Finished.")	# Debug/Performance Monitoring
print("Start Time: " + str(hour) + ':' + str(minute))	# Debug/Performance Monitoring
print("Finish Time: " + str(hour_fn) + ':' + str(minute_fn))	# Debug/Performance Monitoring
#print('')	# Debug/Performance Monitoring
#print("N50: " + N50)

# Proof that this thing actually attempted to do something useful.
file = open("rand3500n50.txt", "w")
file.write(N50)
file.close()
print("Written N50 to file.")