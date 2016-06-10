#!/usr/bin/env python3
# This is a basic template that I made to quickly write some python scripts to take input files, do some functions and then create an ouput file

import os, sys, getopt

################################################################################
# functions
def function(input):
    
    with open (input, 'r') as fh:
        inputfile = fh.read()
    
    ### do some stuff ###
    
    return item


################################################################################
# Main script
if __name__ == "__main__":

    #read commandline arguments and create variables for input and output file names
    try:
        myopts, args = getopt.getopt(sys.argv[1:],'i:o:')
    except getopt.GetoptError as e:
        print (str(e))
        print ('usage: {} -i input -o output' .format(sys.argv[0]))
        sys.exit(2)
    for option, argument in myopts:
        if option == '-i':
            input = argument
        elif option == '-o':
            output = argument

    #print summary of arguments to terminal
    print ('\n\ninput file:\t{}\noutput file:\t{}\n'  .format(input, output))

    #give input to function and return output
    returned_item = function(input)



    with open(output, "wt") as out_file:
        out_file.write('returned_item')
        out_file.close()
