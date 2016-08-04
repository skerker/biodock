### give this script a file directory containing isolate directories. Each isolate directory contains .gff files.
#!/usr/bin/python

######################################################################
# Housekeeping
######################################################################
import time
date = time.strftime('%H%M%S_%d%b%Y')
import os
temp_file_dir = './.TEMP_FILES_' + date
if not os.path.exists(temp_file_dir):
	os.makedirs(temp_file_dir)


######################################################################
# Functions
######################################################################
import sys, getopt

# function to collect arguments from the command line
def CLI_args(argv):
	inputfile = ''
	search_term = ''
	outputfile = ''
	try:
		opts, args = getopt.getopt(argv,"hi:t:o:",["ifile=","term=","ofile="])
	except getopt.GetoptError:
		print 'X.py -i <inputfile> -t <search term> -o <outputfile>'
		sys.exit(2)
	if len(opts) == 0:
		print ('No arguments supplied\n')
	for opt, arg in opts:
		if opt == '-h':
			print 'X.py -i <inputfile> -o <outputfile>'
			sys.exit()
		elif opt in ('-i', '--ifile'):
			inputfile = arg
		elif opt in ('-t', '--term'):
			search_term = arg
		elif opt in ('-o', '--ofile'):
			outputfile = arg
	print ('Input file: {}\nSearch term: {}\nOutput file: {}\n' .format(inputfile, search_term, outputfile))
	return inputfile, search_term, outputfile


######################################################################
# Program
######################################################################

if __name__ == "__main__":

	# get the input directory containing the .gff files
	# get the output file to write results to
	inputfile, search_term, outputfile = CLI_args(sys.argv[1:])

	# get dictionary of all .gff files in genome directory (name and file path)
	list_of_files = {}
	for (dirpath, dirnames, filenames) in os.walk(inputfile):
		for filename in filenames:
			if filename.endswith('.gff'):
				list_of_files[filename] = os.sep.join([dirpath, filename])

	# get the number of .gff files
	print ('number of .gff files: {}\n' .format(len(list_of_files)))

	# read each .gff file line by line and look for keyword match
	# count number of files the contain a match
	match_in_file, match_counter, match_file_counter = 0, 0, 0
	import re
	for gff_file in list_of_files:
		# print the name of each file in the outfile
		with open (outputfile, 'a') as outfh:
			outfh.write('{}' .format(gff_file))
			outfh.close()
		# look at each line to see if there is a match to the search term
		with open (list_of_files[gff_file], 'r') as fh:
			for line in fh:
				match = re.search(search_term, line, re.IGNORECASE)
				# if match found, capture what the predicted product is in the .gff file and the start position
				if match:
					match_in_file = 1
					match_counter += 1
					start_position = line.partition('CDS\t')[2]
					match = re.search('product\=.*', line)
					if match:
						with open (outputfile, 'a') as outfh:
							outfh.write('\t{}\t{}' .format(match.group(), start_position.partition('\t')[0]))
							outfh.close()
					else:
						print ('failed to capture predicted protein, investigate sample: {}' .format(gff_file))
						sys.exit()
			if match_in_file == 1:
				with open (outputfile, 'a') as outfh:
					outfh.write('\n')
					outfh.close()
				match_file_counter += 1
				match_in_file = 0
			# if no match in gff file, state in outfile
			elif match_in_file == 0:
				with open (outputfile, 'a') as outfh:
					outfh.write('\tNO MATCH for {}\n' .format(search_term))
					outfh.close()

		fh.close()
	print ('total number of matches: {} \n' .format(match_counter))
	print ('number of files with one or more matches: {} \n' .format(match_file_counter))



######################################################################
# File outputs and tidy
######################################################################

#remove temporary files
import shutil
shutil.rmtree(temp_file_dir)
