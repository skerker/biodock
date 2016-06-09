#!/usr/bin/env python3
# Takes a SAM file and set of read names. The script creates a new SAM file that contains only the reads from the supplied file (if present)
# Note: need to copy headers over from original SAM file manually
import os, sys, getopt, re

################################################################################
# functions
def function(input1, input2):

    new_sam_file = []
    orig_sam_file = {}
    with open(input2, 'r') as fh2:
        for line in fh2:
            match = re.search('^\S+\s\S+', line)
            if match:
                key = match.group()
                val = line.rstrip('\n')
                orig_sam_file[str(key)] = val

    with open (input1, 'r') as fh:
        inputfile1 = fh.readlines()
        for line in inputfile1:
            line = line.strip('\n')
            if line in orig_sam_file:
                new_sam_file.append(orig_sam_file[line])

    return new_sam_file

################################################################################
# main script
if __name__ == "__main__":

    #read commandline arguments and create variables for input and output file names
    try:
        myopts, args = getopt.getopt(sys.argv[1:],'1:2:o:')
    except getopt.GetoptError as e:
        print (str(e))
        print ('usage: {} -1 read_names -2 sam_file -o output' .format(sys.argv[0]))
        sys.exit(2)
    for option, argument in myopts:
        if option == '-1':
            input1 = argument
        elif option == '-2':
            input2 = argument
        elif option == '-o':
            output = argument

    #print summary of arguments to terminal
    print ('\n\ninput file1:\t{}\ninput file2:\t{}\noutput file:\t{}\n'  .format(input1, input2, output))

    #give input to function and return permutation list
    returned_item = function(input1, input2)



    with open(output, "wt") as out_file:
        for line in returned_item:
            out_file.write('{}\n' .format(line))
        out_file.close()
