### give this script a file directory containing isolate directories. Each isolate directory contains .gff files.
### give this script two search terms (e.g. transposase and lactamase), it will see if any instance of term1 occurs within a set distance of term2 (set to 1000).
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
	search_term1 = ''
	search_term2 = ''
	outputfile = ''
	try:
		opts, args = getopt.getopt(argv,"hi:a:b:o:",["ifile=","term1=","term2=","ofile="])
	except getopt.GetoptError:
		print 'X.py -i <inputfile> -t1 <search term1> -t2 <search term2> -o <outputfile>'
		sys.exit(2)
	if len(opts) == 0:
		print ('No arguments supplied\n')
	for opt, arg in opts:
		if opt == '-h':
			print 'X.py -i <inputfile> -o <outputfile>'
			sys.exit()
		elif opt in ('-i', '--ifile'):
			inputfile = arg
		elif opt in ('-a', '--term1'):
			search_term1 = arg
		elif opt in ('-b', '--term2'):
			search_term2 = arg
		elif opt in ('-o', '--ofile'):
			outputfile = arg
	print ('Input file: {}\nSearch term1: {}\nSearch term2: {}\nOutput file: {}\n' .format(inputfile, search_term1, search_term2, outputfile))
	return inputfile, search_term1, search_term2, outputfile


######################################################################
# Program
######################################################################

if __name__ == "__main__":

	# get the input directory containing the .gff files
	# get the two search terms
	# get the output file to write results to
	inputfile, search_term1, search_term2, outputfile = CLI_args(sys.argv[1:])

	# get dictionary of all .gff files in genome directory (name and file path)
	list_of_files = {}
	for (dirpath, dirnames, filenames) in os.walk(inputfile):
		for filename in filenames:
			if filename.endswith('.gff'):
				list_of_files[filename] = os.sep.join([dirpath, filename])

	# get the number of .gff files
	print ('number of .gff files: {}\n' .format(len(list_of_files)))

	# set up outfile
	with open (outputfile, 'a') as outfh:
		outfh.write('Isolate\t{}\t{}\tDistance between products\tStart position for product A\n' .format(search_term1, search_term2))
		outfh.close()

	# read each .gff file line by line and look for keyword match
	# count number of files the contain a match
	import re
	for gff_file in list_of_files:
		# look at each line to see if there is a match to the first search term
		with open (list_of_files[gff_file], 'r') as fh:
			for line in fh:
				match = re.search(search_term1, line, re.IGNORECASE)
				# if match found, capture what the predicted product is in the .gff file and the start position
				if match:
					termA_line = line.partition('CDS\t')[2]
					termA_product = re.search('product\=.*', line)
					termA_start_position = termA_line.partition('\t')[0]
					# if there is a match and the product has been captured, commence search for second key word
					if termA_product:
						# open a new instance of the gff file and go through line by line, find termB
						with open (list_of_files[gff_file], 'r') as fh2:
							for fh2_line in fh2:
								fh2_match = re.search(search_term2, fh2_line, re.IGNORECASE)
								# if termB match found, capture what the predicted product is in the .gff file and the start position
								if fh2_match:
									termB_line = fh2_line.partition('CDS\t')[2]
									termB_product = re.search('product\=.*', fh2_line)
									termB_start_position = termB_line.partition('\t')[0]

									# now we have start position for termA and termB, see if they are close
									distance_between_products = abs(float(termA_start_position) - float(termB_start_position))
									if (distance_between_products < 1000):
										# print result in the outfile
										with open (outputfile, 'a') as outfh:
													outfh.write('{}\t{}\t{}\t{}\t{}\n' .format(gff_file, termA_product.group(), termB_product.group(), distance_between_products, termA_start_position))
													outfh.close()
									# one the search term B has been found in a line and the location compared to the instance of search term A, the loop continues on to the next line of the file and continues until all instances of search term B have been found in the gff file and compared to the the location of the first instance of search term A in the gff file.
						fh2.close()
					else:
						print ('failed to capture predicted protein for search term A, investigate sample: {}' .format(gff_file))
						sys.exit()
		fh.close()

######################################################################
# File outputs and tidy
######################################################################

#remove temporary files
import shutil
shutil.rmtree(temp_file_dir)
