#!/usr/bin/python
# this script takes a sequence and uses the GeneMark exon prediction tool to find the number of genes/exons in the input sequence

import sys, getopt, re, os

################################################################################
## Path to GeneMark directory containing the human matrices ('human_x_x.mtx'):
GM_path = '/Users/willrowe/Documents/genemark/'


################################################################################
# FUNCTIONS

## A. GC content calculation (req. for GeneMark)
def GC_content_calculation(sequence):
    # turn sequence string into list - each base is an item in the list
    base_list = list(sequence)
    gc = ['G', 'C', 'g', 'c']
    # list comprehension - if a base in the base_list matches a g/c, add it to a new list
    gc_list = [x for x in base_list if x in gc]
    # find length of gc list and base list, convert to a float value and calculate gc content
    gc_content = (float(len(gc_list)) / float(len(base_list))) * 100
    return gc_content


## B. run GeneMark (gmhmme2)
def GeneMark(input, gc_content, output):
    if gc_content < 43:
        matrix = 'human_00_43.mtx'
    elif gc_content < 49:
        matrix = 'human_43_49.mtx'
    elif gc_content >= 49:
        matrix = 'human_49_99.mtx'
    GM_command = ('gmhmme2 -m {}{} -o {} {}' .format(GM_path, matrix, output, input))
    try:
        os.system( GM_command )
        result = 1
        return result
    except:
        print ('\nGeneMark failed\n')
        sys.exit(2)

## C. GeneMark result parser
def GeneMarkParser(gm_output_file):
    ## parse GeneMark results file to find number of predicted genes and introns:
    result = 0
    try:
        with open (gm_output_file, 'r') as fh:
            file = fh.readlines()
            fh.close()
        ## GeneMark outputs exons on inividual lines with a unique start (whitespace followed by a digit), the following regex will only capture the found exons
        exon_list = [x for x in file if (re.match('\s+\d', x))]
        if len(exon_list) == 0:
            print ('\nGeneMark found no exons\n')
            sys.exit(2)
        for i in exon_list:
            initial_exons = [x for x in exon_list if (re.search('Initial', x))]
            single_exons = [x for x in exon_list if (re.search('Single', x))]
            internal_exons = [x for x in exon_list if (re.search('Internal', x))]
        if len(exon_list) == 1 and len(single_exons) == 1:
            print ('1 gene was found, containing 1 exon -> query sequence likely to be a single gene transcript / cDNA.\n')
            result = 1
        if len(exon_list) == 1 and len(internal_exons) == 1:
            print ('1 gene was found, containing 1 exon -> query sequence likely to be a single gene transcript / cDNA.\n')
            result = 1
        if len(exon_list) > 1 and len(initial_exons) > 1:
            print ('multiple genes were found -> query sequence could be genomic DNA.\n')
        if len(exon_list) > 1 and len(initial_exons) <= 1:
            print ('1 gene was found with multiple exons -> query sequence could be genomic DNA or unspliced transcript.\n')
        check = 'positive'
        return result
    ## use except to catch an error if GeneMark did not complete properly
    except IOError:
        print ('\nGeneMark failed\n')
        sys.exit(2)


################################################################################
# SCRIPT

if __name__ == '__main__':

## parse CL inputs
    try:
        myopts, args = getopt.getopt(sys.argv[1:], 'i:o:')
    except getoptGetoptError as e:
        print (str(e))
        print ('usage: {} -i input -o output' .format(sys.argv[0]))
        sys.exit(2)
    for option, argument in myopts:
        if option == '-i':
            input = argument
        elif option == '-o':
            output = argument
    print ('\nThis script takes a single DNA sequence and produces an output file containing the number of predicted exons.\ninput file:\t{}\noutput file:\t{}\n'  .format(input, output))

## open input, extract sequence from file and store as string
    with open (input, 'r') as fh:
        input_file = fh.readlines()
        fh.close()
    sequence = ''
    for i in input_file:
        if re.match('\>', i):
            sequence_header = i
        else:
            sequence += i.replace('\n','')

## calculate GC content of sequence, then run GeneMark using appropriate GC matrix
    gc_content = GC_content_calculation(sequence)
    GeneMark(input, gc_content, output)
    result = GeneMarkParser(output)