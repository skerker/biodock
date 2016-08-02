#!/usr/bin/env python3
import os, sys, getopt, re, glob, subprocess, QC_data, MAP_data, TPM_calculation
from datetime import datetime
start_time = datetime.now()

"""
################################################################################
### INFO

This is the RNA-seq pipeline for the Hinton Lab.

The pipeline steps:
    1   -   QC      -   initial QC of reads, quality based trimming of reads, second QC of reads and taxonomic binning of reads
    2   -   MAP     -   map reads to D23580_liverpool or ST474 Salmonella typhimurium genomes using Bowtie2, filter alignments for MAPQ score >10, generate mapping stats
    3   -   COUNT   -   use HTSeq to count reads per CDS and ncRNA feature in gff file
    4   -   TPM     -   generate TPM values for each feature


version = 0.1
author  = Will Rowe
email = will.rowe@liverpool.ac.uk




################################################################################
### NOTES

this pipeline is designed for the CGR cluster - it creates a workdir in SCRATCH for each job and then transfers results back to BRICK (to the supplied output dir)
input file must be filename.fastq.gz - no checks exist for this yet!
pipeline breaks with no error message if it tries to use more cores than the machine has - no checks exist for this yet either!




################################################################################
### How to run using GNU Parallel:

First, set server list as a variable:
    PARALLEL_HOSTS=ada05,ada06,ada07

Second, if outputing python STDOUT to log file, make sure directory for logs exist or error will be thrown

Run GNU parallel by piping list of files that you want to run the pipeline on into the GNU parallel command:
    ls /pub46/willr/000_HOME/0003_PROJECTS/rocio/00_seq_data/D23580/*.gz | parallel --progress --workdir $PWD -j 10% --delay 2.0 -S $PARALLEL_HOSTS "python /pub46/willr/000_HOME/0005_RNA-SEQ-PIPELINE/02_PIPELINE_FILES/00_RNAseq_pipeline.py -i {} -o $PWD/RESULTS -q Y -m Y -c Y -t Y -r D23 >> $PWD/LOGS/{/.}.LOG"

Breakdown of GNU command:
       -j 10%            -   how much to load the receiving servers (number of cores to use)
       --delay 2.0       -   run with delay between jobs (otherwise there are conflicts with scratch directory naming)   *** REQUIRED ***
       -S                -   server list to use
       --progress        -   prints job progress to STDOUT
       --filter-hosts    -   ignore nodes that are down (this is currently broken in current gnu parallel release)




################################################################################
### TO DO

need to add checks - for fastq.gz input, same file names, number of threads etc.
smarter way to assign threads for each process
add customisation of trimmomatic, mapping etc.


"""



################################################################################
### SET PATHS AND DEFAULTS
# directory name in scratch space
scratch = '/scratch/will___' + start_time.strftime('%H%M%S')

# path to software on cluster
htseq_path = '/usr/local/bin/htseq-count'

# bowtie2 reference indices
bt2_D23_index = '/pub46/willr/000_HOME/0002_REF_DATA/0002_BT2_REFS/D23580_liv/D23580_liv_chrom+4plasmids'
bt2_474_index = '/pub46/willr/000_HOME/0002_REF_DATA/0002_BT2_REFS/474/474_chrom+3plasmids'

# annotation files
gff_D23 = '/pub46/willr/000_HOME/0002_REF_DATA/0003_ANNOTATIONS/D23580_liv/D23580_liv_2016.gff'
gff_474 = ''

# default settings
threads = '20'
run_QC = False
length_cutoff = '20'
run_MAP = False
bt2_index = bt2_D23_index
run_COUNT = False
gff_file = gff_D23
run_TPM = False




################################################################################
### FUNCTIONS
"""
this function creates a file list of all the fastq.gz files in a given directory.
** this is not currently used in the pipeline... **
"""
def create_file_list(input_data_dir):
    print ('*Path to input files\t==>\t{}\n\n*Files that the pipeline will process:' .format(input_data_dir))
    with open('./file_list.txt', 'w') as file_list:
        for filename in glob.glob(os.path.join(input_data_dir, '*.fastq.gz')):
            print ('{}' .format(filename))
            file_list.write('{}\n' .format(filename))
        file_list.close()
        file_list_success = 1
    return file_list_success


"""
this function gets all the scratch dirs in place on the CGR servers and returns the results_sub_dir for subsequent functions
"""
def pipeline_setup(scratch, input_fastq_file, results_dir):
    if not os.path.exists(scratch):
        try:
            os.makedirs(scratch)
        except Exception, e:
            print >> sys.stderr, 'Can\'t create the working directory in scratch . . .'
            print >> sys.stderr, 'Exception: %s' % str(e)
            sys.exit(1)
    try:
        os.chdir(scratch)
    except Exception, e:
        print >> sys.stderr, 'Can\'t cd into the server\'s scratch dir . . .'
        print >> sys.stderr, 'Exception: %s' % str(e)
        sys.exit(1)
    file_basename = os.path.basename(input_fastq_file)
    results_sub_dir = scratch + '/RNA-seq-pipeline___' + file_basename + '___' + start_time.strftime('%H%M%S')
    print ('\t*Directory for output files does not exist, making it now . . .\n\n')
    try:
        os.makedirs(results_sub_dir)
    except Exception, e:
        print >> sys.stderr, 'Failed to create the results subdirectory'
        print >> sys.stderr, 'Exception: %s' % str(e)
        sys.exit(1)

    # get all results files ready
    QC_results_file = scratch + '/' + file_basename + '.___QCresults.csv'
    MAPPING_results_file = scratch + '/' + file_basename + '.___MAPPINGresults.csv'
    COUNT_results_file = scratch + '/' + file_basename + '.___COUNTresults.csv'
    TPM_results_file = scratch + '/' + file_basename + '.___TPMresults.csv'

    return (results_sub_dir, file_basename, QC_results_file, MAPPING_results_file, COUNT_results_file, TPM_results_file)


"""
this function runs HTSeq on the mapped data, it then writes a results file to the output directory (as well as giving info to STDOUT)
"""
def count_data(sorted_filtered_bam, gff_file, results_sub_dir, COUNT_results_file):
    # commands
    count_CDS_features = results_sub_dir + '/CDS_features.count'
    count_ncRNA_features = results_sub_dir + '/ncRNA_features.count'
    htseq_CDS_cmd = htseq_path + ' -f bam -t CDS -i ID -m intersection-nonempty -q ' + sorted_filtered_bam + ' ' + gff_file + ' > ' + count_CDS_features
    htseq_ncRNA_cmd = htseq_path + ' -f bam -t ncRNA -i ID -m intersection-nonempty -q ' + sorted_filtered_bam + ' ' + gff_file + ' > ' + count_ncRNA_features

    # run htseq-count
    try:
        os.system(htseq_CDS_cmd)
        os.system(htseq_ncRNA_cmd)
    except Exception, e:
        print >> sys.stderr, 'Can\'t run htseq-count . . .'
        print >> sys.stderr, 'Exception: %s' % str(e)
        sys.exit(1)

    # combine count files
    filenames = [count_CDS_features, count_ncRNA_features]
    with open(COUNT_results_file, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    return COUNT_results_file



################################################################################
### main script
if __name__ == "__main__":


    # read commandline arguments and create variables
    try:
        myopts, args = getopt.getopt(sys.argv[1:],'i:o:q:l:m:c:t:r:')
    except getopt.GetoptError as e:
        print (str(e))
        print ('usage: {} -i single fastq file (fastq.gz) -o directory for output files -q Y -m Y -c Y -t Y-r D23' .format(sys.argv[0]))
        sys.exit(2)
    for option, argument in myopts:
        if option == '-i':
            input_fastq_file = os.path.abspath(argument)
            try:
                os.path.isfile(input_fastq_file)
            except Exception, e:
                print >> sys.stderr, 'Supplied input does not appear to be a file . . .'
                print >> sys.stderr, 'Exception: %s' % str(e)
                sys.exit(1)
        if option == '-o':
            results_dir = os.path.abspath(argument)
            if not os.path.exists(results_dir):
                try:
                    os.makedirs(results_dir)
                except Exception, e:
                    print >> sys.stderr, 'Can\'t create the supplied results directory . . .'
                    print >> sys.stderr, 'Exception: %s' % str(e)
                    sys.exit(1)
        if option == '-q':
            if argument == 'Y':
                run_QC = True
            elif argument == 'N':
                run_QC = False
            else:
                print >> sys.stderr, 'Please specify Y or N for the q option . . .'
                sys.exit(1)
        # get length cutoff value
        if option == '-l':
            length_cutoff = argument
            try:
                length_cutoff.isdigit()
            except Exception, e:
                print >> sys.stderr, 'The supplied length cut-off value does not appear to be a digit . . .'
                print >> sys.stderr, 'Exception: %s' % str(e)
                sys.exit(1)
        if option == '-m':
            if argument == 'Y':
                run_MAP = True
            elif argument == 'N':
                run_MAP = False
            else:
                print >> sys.stderr, 'Please specify Y or N for the m option . . .'
                sys.exit(1)
        if option == '-c':
            if argument == 'Y':
                run_COUNT = True
            elif argument == 'N':
                run_COUNT = False
            else:
                print >> sys.stderr, 'Please specify Y or N for the c option . . .'
                sys.exit(1)
        if option == '-t':
            if argument == 'Y':
                run_TPM = True
            elif argument == 'N':
                run_TPM = False
            else:
                print >> sys.stderr, 'Please specify Y or N for the t option . . .'
                sys.exit(1)
        if option == '-r':
            if argument == '474':
                bt2_index = bt2_474_index
                gff_file = gff_474
            elif argument == 'D23':
                bt2_index = bt2_D23_index
                gff_file = gff_D23
            else:
                print >> sys.stderr, 'Please specify 474 or D23 for the r option . . .'
                sys.exit(1)


    # print some messages --- could use this to check for multiple same files in the source input dir??
    border = '#################################################################################################'
    print ('{}{}\n## This is the test RNA-seq pipeline for the Vertis data\n{}{}' .format(border, border, border, border))
    print ('\n\n\tMore useful messages to follow here . . .\n\n')


    # set up working dir, results directory etc. and cd into scratch space
    print ('{}\n## Pipeline setup\n{}' .format(border, border))
    print ('\n\n\t*Setting working directory\t==>\t{}\n\n\n\t*Setting results directory\t==>\t{}\n\n' .format(scratch, results_dir))
    try:
        results_sub_dir, file_basename, QC_results_file, MAPPING_results_file, COUNT_results_file, TPM_results_file = pipeline_setup(scratch, input_fastq_file, results_dir)
    except Exception, e:
        print >> sys.stderr, 'Can\'t setup required directories . . .'
        print >> sys.stderr, 'Exception: %s' % str(e)
        sys.exit(1)


    # if QC wanted:
    if run_QC == True:
        print ('{}\n## QC-ing data\n{}' .format(border, border))
        print ('\n\n\t*Input file\t==>\t{}\n\n' .format(input_fastq_file))
        print ('\t*Spawning FastQC job to check raw data . . .\n\n\n\t*Running Trimmomatic (with os.system) to check raw data (length cutoff = {}). . .\n\n' .format(length_cutoff))

        try:
            trimmed_data = QC_data.qc_data(input_fastq_file, length_cutoff, results_sub_dir, file_basename, QC_results_file)
        except Exception, e:
            print >> sys.stderr, 'Can\'t run QC . . .'
            print >> sys.stderr, 'Exception: %s' % str(e)
            sys.exit(1)

        print ('\n\n\t*Finished QC of RNA-seq data')
        print ("\t--- %s seconds ---\n\n" % (datetime.now() - start_time))

    else:
        print ('\n\n\tThe QC option was not selected, skipping QC . . .\n\n')


    # if mapping wanted:
    if run_MAP == True:
        if run_QC == True:
            input_fastq_file = trimmed_data
            MAPPING_results_file += '.trimmed_data.csv'
            COUNT_results_file += '.trimmed_data.csv'
            TPM_results_file += '.trimmed_data.csv'

        print ('{}\n## Mapping data\n{}' .format(border, border))
        print ('\n\n\t*Input file\t==>\t{}\n\n' .format(input_fastq_file))
        print ('\t*Reference\t==>\t{}\n\n' .format(bt2_index))
        print ('\t*Running Bowtie2 with input file and pre-compiled reference index . . .\n\n\n\t\t*Samtools:\n\n')

        try:
            sorted_filtered_bam = MAP_data.map_data(input_fastq_file, bt2_index, results_sub_dir, scratch, MAPPING_results_file)
        except Exception, e:
            print >> sys.stderr, 'Can\'t run mapping . . .'
            print >> sys.stderr, 'Exception: %s' % str(e)
            sys.exit(1)

        print ('\n\n\t*Finished mapping of RNA-seq data')
        print ("\t--- %s seconds ---\n\n" % (datetime.now() - start_time))


        # if counting wanted:
        if run_COUNT == True:
            print ('{}\n## Generating read counts\n{}' .format(border, border))
            print ('\n\n\t*Input file\t==>\t{}\n\n' .format(sorted_filtered_bam))
            print ('\t*Annotation file\t==>\t{}\n\n' .format(gff_file))
            print ('\t*Running HTSeq-count with input file . . .\n\n')

            try:
                read_counts = count_data(sorted_filtered_bam, gff_file, results_sub_dir, COUNT_results_file)
            except Exception, e:
                print >> sys.stderr, 'Can\'t run htseq-count . . .'
                print >> sys.stderr, 'Exception: %s' % str(e)
                sys.exit(1)
        else:
            print ('\n\n\tThe read counting option was not selected, skipping htseq-count . . .\n\n')



        # if tpm values wanted:
        if run_TPM == True:
            print ('{}\n## Performing TPM calculations\n{}' .format(border, border))

            try:
                TPM_result = TPM_calculation.TPM_calc(gff_file, read_counts, TPM_results_file)
            except Exception, e:
                print >> sys.stderr, 'Can\'t run TPM calculation . . .'
                print >> sys.stderr, 'Exception: %s' % str(e)
                sys.exit(1)

            print ('\n\n\t*Finished calculating TPM values')
            print ("\t--- %s seconds ---\n\n" % (datetime.now() - start_time))
        else:
            print ('\n\n\tThe TPM calculation option was not selected, skipping TPM calculation . . .\n\n')



    else:
        print ('\n\n\tThe mapping option was not selected, skipping mapping (and read counting + TPM calculations) . . .\n\n')











### DELETE UNWANTED TEMP FILES - eg. bt2 input etc.


    # transfer files from scratch back to results dir (on brick) and delete scratch dir
    from distutils.dir_util import copy_tree
    copy_tree(scratch, results_dir)
    import shutil
    shutil.rmtree(scratch)
