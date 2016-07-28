#!/usr/bin/env python3
import os, sys, re, getopt


"""
################################################################################
### INFO

This script generates a list of TPM values from a set of count data and an annotation file

version = 0.2
author  = Will Rowe
email = will.rowe@liverpool.ac.uk


################################################################################
### NOTES

Only tested with HTSeq count output and GFF2 input files
This script assumes that the identifers used in the HTSeq count file are the same as the feature_identifier specified in the GFF file. E.G. both using the ID attribute tag.


################################################################################
### BACKGROUND

Extract from TPM blog post (https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/):

	Transcripts per million (TPM) is a measurement of the proportion of transcripts in your pool of RNA.

	Since we are interested in taking the length into consideration, a natural measurement is the rate, counts per base (X_i / L_i).
	As you might immediately notice, this number is also dependent on the total number of fragments sequenced.
	To adjust for this, simply divide by the sum of all rates and this gives the proportion of transcripts i in your sample.
	After you compute that, you simply scale by one million because the proportion is often very small and a pain to deal with.

	rate = counts per base of a transcript
	tpm = (rate for transcript of interest) * (1.0 / total rate for all transcripts) * (10**6)


"""


################################################################################
### FUNCTIONS
"""
feature length function
returns a dicitonary of lengths for each feature in GFF annotation file
"""
def get_feature_lengths(gff_annotation):
	feature_lengths = {}

	# set the gff_attribute_delimitereter and the GFF attribute to use as the feature identifier
	gff_attribute_delimiter = '='
	gff_attribute = 'ID'
	feautre_identifier = gff_attribute + gff_attribute_delimiter

	# open gff file and extract the start, end and feature identifier from each line
	# calculate the lengths of each feature add to the feature_lengths dictionary using the identifier as the key
	with open(gff_annotation) as gtf:
		for line in gtf:
			line = line.split("\t")
			start = int(line[3])
			end = int(line[4])
			gff_feature_attributes = line[8]
			if gff_feature_attributes.count(feautre_identifier) < 1:
				print ('Error: the feature identifier {} can\'t be found' .format(feautre_identifier))
				sys.exit(1)
			else:
				gff_feature_attributes = gff_feature_attributes.split(";")
				for attribute in gff_feature_attributes:
					#get rid of "" if present in the GFF attributes
					if feautre_identifier in attribute:
						feature = attribute.split(gff_attribute_delimiter)[-1].replace("\"", "")
				if feature not in feature_lengths.keys():
					feature_lengths[feature] = end - start
				elif feature in feature_lengths.keys():
					feature_lengths[feature] += end - start

	return feature_lengths


"""
feature count function
returns a dicitonary of counts for each feature in HTSeq count file
"""
def get_feature_counts(htseq_count_output):
	feature_counts = {}

	# open count file and add each count to the feature_counts dictionary using the feature name as the key
	with open(htseq_count_output) as counts:
		for line in counts:
			line = line.split("\t")
			feature = line[0]
			count = int(line[1])
			# don't include extra information from htseq (eg. __no_feature) in the count dictionary
			if re.match("__", feature):
				continue
			if feature not in feature_counts.keys():
				feature_counts[feature] = count
			else:
				print ('Error: the feature {} occured more than once in the HTSeq count file!' .format(feature))
				sys.exit(1)

	return feature_counts


"""
TPM calculation function
generates a .tsv file with following rows: feature name\t	feature length\t	feature count\t	feature TPM value
"""
def TPM_calc(gff_annotation, htseq_count_output, TPM_output_file):
	feature_lengths = get_feature_lengths(gff_annotation)
	feature_counts = get_feature_counts(htseq_count_output)

	# calculate the sum of all rates for the sample
	sum_of_rates = 0
	for feature in feature_counts.keys():
		# check to see if feature in counts file is present in the gff file
		if feature not in feature_lengths.keys():
			print ('Error: the feature {} is not present in the supplied annotation file (but is in the count file)' .format(feature))
			sys.exit(1)
		else:
			sum_of_rates += float(feature_counts[feature]) / feature_lengths[feature]

	# calculate TPM value per feature
	for feature in feature_counts.keys():
		tpm_value = (float(feature_counts[feature]) / feature_lengths[feature]) * (1.0 / sum_of_rates) * (10**6)
		with open (TPM_output_file, 'a') as fh:
			fh.write ('{}\t{}\t{}\t{}\n' .format(feature, feature_lengths[feature], feature_counts[feature], str(tpm_value)))

	return TPM_output_file


################################################################################
### MAIN SCRIPT
if __name__ == "__main__":

    # read commandline arguments and create variables
    A_met, C_met, O_met = (False,)*3
    try:
        myopts, args = getopt.getopt(sys.argv[1:],'a:c:o:')
    except getopt.GetoptError as e:
        print (str(e))
        print ('usage: {} -a [annotation file] -c [count file] -o [output file]' .format(sys.argv[0]))
        sys.exit(2)
    for option, argument in myopts:
        if option == '-a':
            gff_annotation = os.path.abspath(argument)
            try:
                os.path.isfile(gff_annotation)
            except Exception, e:
                print >> sys.stderr, 'Supplied annotation does not appear to be a file . . .'
                print >> sys.stderr, 'Exception: %s' % str(e)
                sys.exit(1)
		A_met = True
        if option == '-c':
            htseq_count_output = os.path.abspath(argument)
            try:
                os.path.isfile(htseq_count_output)
            except Exception, e:
                print >> sys.stderr, 'Supplied count file does not appear to be a file . . .'
                print >> sys.stderr, 'Exception: %s' % str(e)
                sys.exit(1)
		C_met = True
        if option == '-o':
        	TPM_output_file = os.path.abspath(argument)
		O_met = True
    if not any((A_met, C_met, O_met)):
	print ('options a, c and o must all be supplied!\n')
	print ('usage: {} -a [annotation file] -c [count file] -o [output file]' .format(sys.argv[0]))
	sys.exit(2)


    # run the TPM calculation
    try:
        TPM_results = TPM_calc(gff_annotation, htseq_count_output, TPM_output_file)
    except Exception, e:
        print >> sys.stderr, 'Can\'t calculate TPM values . . .'
        print >> sys.stderr, 'Exception: %s' % str(e)
        sys.exit(1)
