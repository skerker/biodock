#!/usr/bin/env python3
################################################################################
### What does this script do? At the moment it does a couple of things
#1. creates a file list for running the RNA-seq pipeline with gnu parallel. Give it a directory and it will return a file list of all suitable files within that directory.
#2. runs a set of system commands for QC-ing the RNA-seq data. One file is accepted and each QC step is spawned in parallel.
#3. mapping....
### this pipeline creates a workdir in SCRATCH for each job and then transfers results back to BRICK (to the supplied output dir)
### input file must be filename.fastq.gz !!!!!!

################################################################################
### How to run using GNU Parallel:
# First, set server list as a variable e.g.:
#       PARALLEL_HOSTS=ada05,ada06,ada07

# Second, if outputing python STDOUT to log file, make sure directory for logs exist or error will be thrown

# Run GNU parallel by piping list of files that you want to run the pipeline on into the GNU parallel command:
#       ls /pub46/willr/000_HOME/0003_PROJECTS/rocio/00_seq_data/D23580/*.gz | parallel --progress --workdir $PWD -j 20% --delay 2.0 -S $PARALLEL_HOSTS "python $PWD/00_RNAseq_pipeline.py -i {} -o $PWD/ROCIO_D23_prelim_run -q Y -m Y >> $PWD/ROCIO_D23_LOGS/{/.}.LOG"

#  Breakdown of GNU command:
#       -j 20%            -   how much to load the receiving servers    *** this does not take into account the additional cpus requested by the python pipeline (threads variable * 2)
#       --delay 2.0       -   run with delay between jobs (otherwise there are conflicts with scratch directory naming)   *** REQUIRED ***
#       -S                -   server list to use
#       --progress        -   prints job progress to STDOUT
#       --filter-hosts    -   ignore nodes that are down (this is currently broken in current gnu parallel release)


################################################################################
##### TO DO
### need to add checks - for fastq.gz input, same file names etc.
### smarter way to assign threads for each process
### add customisation of trimmomatic
### add customisation of mapping (eg. change reference etc.)


################################################################################
### load modules and get start time
import os, sys, getopt, re, glob, subprocess
from datetime import datetime
start_time = datetime.now()


################################################################################
### set paths and variable defaults
scratch = '/scratch/will___' + start_time.strftime('%H%M%S')
fastqc_path = '/usr/local/share/FastQC/fastqc'
kraken_path = '/pub46/willr/000_HOME/test_RNA-seq-pipeline/01_BIN/kraken'
kraken_report_path = '/pub46/willr/000_HOME/test_RNA-seq-pipeline/01_BIN/kraken-report'
kraken_db = '/pub46/willr/000_HOME/test_RNA-seq-pipeline/02_INDICES+REFS/minikraken_20141208'
trimm_path = '/pub46/willr/000_HOME/test_RNA-seq-pipeline/01_BIN/trimmomatic-0.36.jar'
bt2_path = '/usr/local/share/bowtie2-2.2.5/bowtie2'
samtools_path = '/usr/local/bin/samtools'

bt2_D23_index = '/pub46/willr/000_HOME/0002_REF_DATA/0002_BT2_REFS/D23580_liv/D23580_liv_chrom+4plasmids'
bt2_474_index = '/pub46/willr/000_HOME/0002_REF_DATA/0002_BT2_REFS/474/474_chrom+3plasmids'

threads = '20'
run_QC = False
run_MAP = False

################################################################################
### functions
## 1. FILE-LIST ** this is not currently used in the pipeline...
# this function creates a file list of all the fastq.gz files in a given directory.
def create_file_list(input_data_dir):
    print ('*Path to input files\t==>\t{}\n\n*Files that the pipeline will process:' .format(input_data_dir))
    with open('./file_list.txt', 'w') as file_list:
        for filename in glob.glob(os.path.join(input_data_dir, '*.fastq.gz')):
            print ('{}' .format(filename))
            file_list.write('{}\n' .format(filename))
        file_list.close()
        file_list_success = 1
    return file_list_success


## 2. SETUP PIPELINE
# this function gets all the scratch dirs in place on the CGR servers and returns the results_sub_dir for subsequent functions.
def pipeline_setup(border, scratch, input_fastq_file, results_dir):
    print ('{}\n## Pipeline setup\n{}' .format(border, border))
    print ('\n\n\t*Setting working directory\t==>\t{}\n\n\n\t*Setting results directory\t==>\t{}\n\n' .format(scratch, results_dir))
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
    return (results_sub_dir)


## 3. QC DATA
# this function runs FastQC and Trimmomatic on the raw fastq data, then runs Kraken and FastQC on the trimmed reads. It then writes a QC results file to the output directory (as well as giving info to STDOUT)
def qc_data(border, input_fastq_file, results_sub_dir, scratch):
    print ('{}\n## QC-ing data\n{}' .format(border, border))
    print ('\n\n\t*Input file\t==>\t{}\n\n' .format(input_fastq_file))
    print ('\t*Spawning FastQC job to check raw data . . .\n\n\n\t*Running Trimmomatic (with os.system) to check raw data . . .\n\n')
    file_basename = os.path.basename(input_fastq_file)
    trimmed_data = scratch + '/trimmed.' + file_basename

    # QC commands
    fastqc_cmd = fastqc_path + ' -t ' + threads + ' ' + input_fastq_file + ' --outdir ' + results_sub_dir
    trimm_cmd = trimm_path + ' SE -threads ' + threads + ' ' + input_fastq_file + ' ' + trimmed_data + ' ' + ' SLIDINGWINDOW:4:20 MINLEN:50' + ' &> ' + results_sub_dir + '/trimm.log'
    kraken_cmd = kraken_path + ' --threads ' + threads + ' --preload --fastq-input --gzip-compressed --db ' + kraken_db + ' ' + trimmed_data + ' | ' + kraken_report_path + ' --db ' + kraken_db + ' > ' + results_sub_dir + '/krakenreport'
    fastqc_cmd_2 = fastqc_path + ' -t ' + threads + ' ' + trimmed_data + ' --outdir ' + results_sub_dir

    # run fastqc subprocess
    processes = []
    with open('%s/fastqc.log' % results_sub_dir, 'w') as fastqc_log:
        p1 = subprocess.Popen(fastqc_cmd, shell=True, stdout=fastqc_log, stderr=fastqc_log)
        processes.append(p1)

    # os.system for trimmomatic
    trimm_log = results_sub_dir + '/trimm.log'
    os.system(trimm_cmd)

    # wait for fastqc subprocesses to complete:
    exit_codes = [p.wait() for p in processes]

    # retrieve QC stats for raw data
    sample_stats = {}
    # retrieve QC stats for raw data - initial fastqc of data
    fastqc_file = file_basename + '_fastqc/fastqc_data.txt'
    fastqc_file = fastqc_file.replace('.fastq.gz','')
    with open('{}/{}' .format(results_sub_dir, fastqc_file), 'r') as fh:
        fastqc_lines = fh.readlines()
        fastqc_base_seq_qual = [line for line in fastqc_lines if re.search('>>Per base sequence quality\t', line)]
        fastqc_len_dist = [line for line in fastqc_lines if re.search('>>Sequence Length Distribution\t', line)]
        fastqc_seq_len = [line for line in fastqc_lines if re.match('Sequence length', line)]
        sample_stats['fastqc'] = fastqc_base_seq_qual, fastqc_len_dist, fastqc_seq_len
    # retrieve QC stats for raw data - trimmomatic summary
    with open('{}/trimm.log' .format(results_sub_dir), 'r') as fh:
        trimm_lines = fh.readlines()
        trimm_summary = [line for line in trimm_lines if re.search(' Surviving: ', line)]
        sample_stats['trimmomatic'] = trimm_summary

    # parse QC stats for raw data
    # parse QC stats for raw data - initial fastqc
    if re.search('pass', str(sample_stats["fastqc"][0])) is not None:
        print ('\t\t*FastQC:\tsequence base quality\t==>\tPASS')
        base_qual = 'PASS'
    else:
        print ('\t\t*FastQC:\tsequence base quality\t==>\tFAIL')
        base_qual = 'FAIL'
    if re.search('pass', str(sample_stats["fastqc"][1])) is not None:
        print ('\t\t*FastQC:\tseq length distribution\t==>\tPASS')
        len_dist = 'PASS'
    else:
        print ('\t\t*FastQC:\tseq length distribution\t==>\tFAIL')
        len_dist = 'FAIL'
    fqc_search = re.compile(r"Sequence length\\t(.*)\\t\\n")
    fastqc_match = re.search(fqc_search, str(sample_stats["fastqc"][2]))
    if fastqc_match:
        print ('\t\t*FastQC:\tinitial sequence length\t==>\t{}' .format(fastqc_match.group(1)))
        seq_length = fastqc_match.group(1)
    else:
        print ('\t\t*FastQC:\tinitial sequence length\t==>\terror with parsing fastqc results\n')
        seq_length = 'ERROR'
    # parse QC stats for raw data - trimmomatic
    trimm_match = re.search('Dropped:\s+\d+\s+\((\S+)\)', str(sample_stats["trimmomatic"]))
    if trimm_match:
        print ('\t\t*Trimmomatic:\tdropped sequences\t==>\t{}' .format(trimm_match.group(1)))
        dropped_seqs = trimm_match.group(1)
    else:
        print ('\t\t*Trimmomatic:\tdropped sequences\t==>\terror with parsing trimmomatic results\n')
        dropped_seqs = 'ERROR'

    # run second fastqc + kraken subprocesses
    print ('\n\n\t*Completed FastQC and Trimmomatic jobs!\n\n\n\t*Running FastQC and and Kraken on trimmed data . . .\n\n')
    processes2 = []
    with open('%s/fastqc2.log' % results_sub_dir, 'w') as fastqc_log2:
        p1 = subprocess.Popen(fastqc_cmd_2, shell=True, stdout=fastqc_log2, stderr=fastqc_log2)
        processes2.append(p1)
    with open('%s/kraken.log' % results_sub_dir, 'w') as kraken_log:
        p2 = subprocess.Popen(kraken_cmd, shell=True, stdout=kraken_log, stderr=kraken_log)
        processes2.append(p2)

    # wait for subprocesses to complete:
    exit_codes = [p.wait() for p in processes2]

    # retrieve QC stats - kraken results for trimmed data
    with open('{}/kraken.log' .format(results_sub_dir), 'r') as fh:
        kraken_lines = fh.readlines()
        kraken_unclassified = [line for line in kraken_lines if re.search('sequences unclassified ', line)]
        sample_stats['kraken'] = kraken_unclassified
        for line in kraken_lines:
            seq_num_match = re.search('\r(\d+) sequences \(.*\) processed in', line)
            if seq_num_match:
                seq_num = seq_num_match.group(1)
                break
            else:
                seq_num = 'Could not parse Kraken file for seq num'

    # retrieve QC stats - fastqc results for trimmed data
    trimmed_basename = os.path.basename(trimmed_data)
    fastqc_file2 = trimmed_basename + '_fastqc/fastqc_data.txt'
    fastqc_file2 = fastqc_file2.replace('.fastq.gz','')
    with open('{}/{}' .format(results_sub_dir, fastqc_file2), 'r') as fh:
        fastqc_lines2 = fh.readlines()
        fastqc_base_seq_qual_trimmed = [line for line in fastqc_lines2 if re.search('>>Per base sequence quality\t', line)]
        fastqc_len_dist_trimmed = [line for line in fastqc_lines2 if re.search('>>Sequence Length Distribution\t', line)]
        fastqc_seq_len_trimmed = [line for line in fastqc_lines2 if re.match('Sequence length', line)]
        sample_stats['fastqc_trimmed'] = fastqc_base_seq_qual_trimmed, fastqc_len_dist_trimmed, fastqc_seq_len_trimmed

    # parse QC stats for results
    # parse QC stats for results - kraken
    kraken_match = re.search('\((\S+)\)', str(sample_stats["kraken"]))
    if kraken_match:
        print ('\t\t*Kraken:\tunclassified sequences\t==>\t{}' .format(kraken_match.group(1)))
        unclassified_reads = kraken_match.group(1)
    else:
        print ('\t\t*Kraken:\tunclassified sequences\t==>\terror with parsing kraken results\n')
        unclassified_reads = 'ERROR'
    # parse QC stats for results - fastqc
    if re.search('pass', str(sample_stats["fastqc_trimmed"][0])) is not None:
        print ('\t\t*FastQC:\tsequence base quality\t==>\tPASS')
        base_qual_trimmed = 'PASS'
    else:
        print ('\t\t*FastQC:\tsequence base quality\t==>\tFAIL')
        base_qual_trimmed = 'FAIL'
    if re.search('pass', str(sample_stats["fastqc_trimmed"][1])) is not None:
        print ('\t\t*FastQC:\tseq length distribution\t==>\tPASS')
        len_dist_trimmed = 'PASS'
    else:
        print ('\t\t*FastQC:\tseq length distribution\t==>\tFAIL')
        len_dist_trimmed = 'FAIL'
    fqc_search = re.compile(r"Sequence length\\t(.*)\\t\\n")
    fastqc2_match = re.search(fqc_search, str(sample_stats["fastqc_trimmed"][2]))
    if fastqc2_match:
        print ('\t\t*FastQC:\ttrimmed sequence length\t==>\t{}' .format(fastqc2_match.group(1)))
        seq_length_trimmed = fastqc2_match.group(1)
    else:
        print ('\t\t*FastQC:\ttrimmed sequence length\t==>\terror with parsing fastqc results\n')
        seq_length_trimmed = 'ERROR'

    # write results to file
    with open('{}/{}.___QCresults.csv' .format(scratch, file_basename), 'w') as fh:
        fh.write('{}\t{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t{}\t\n' .format(file_basename, seq_length, base_qual, len_dist, dropped_seqs, seq_length_trimmed, base_qual_trimmed, len_dist_trimmed, unclassified_reads))

    print ('\n\n\t*Finished QC of RNA-seq data')
    print ("\t--- %s seconds ---\n\n" % (datetime.now() - start_time))
    return trimmed_data


## 3. MAP DATA
# this function runs Bowtie2 on the QCed data (or falls back to the raw fastq data if QC not run), it then writes a results file to the output directory (as well as giving info to STDOUT)
def map_data(border, input_fastq_file, results_sub_dir, scratch):
    print ('{}\n## Mapping data\n{}' .format(border, border))
    print ('\n\n\t*Input file\t==>\t{}\n\n' .format(input_fastq_file))
    print ('\t*Reference\t==>\t{}\n\n' .format(input_fastq_file))
    print ('\t*Running Bowtie2 with input file and pre-compiled reference index . . .\n\n\n\t\t*Samtools:\n\n')

    # commands
    bt2input = scratch + '/bt2.input.fq'
    bt2output_bam = scratch + '/bt2.output.bam'
    flagstats = scratch + '/00_mapping.flagstats'
    idxstats = scratch + '/00_mapping.idxstats'
    gunzip_cmd = 'gzip -cd ' + input_fastq_file + ' > ' + bt2input
    bt2_cmd = bt2_path + ' -x ' + bt2_D23_index + ' -q ' + bt2input + ' --very-sensitive-local -p ' + threads + ' | ' + samtools_path + ' view -bS - > ' + bt2output_bam
    min_map_qual_cmd = samtools_path + ' view -q 30 -b ' + bt2output_bam + ' > ' + bt2output_bam + '.q30_filter.bam'
    sort_bam_cmd = samtools_path + ' sort ' + bt2output_bam + '.q30_filter.bam ' + bt2output_bam + '.q30_filter.sorted'
    index_bam_cmd = samtools_path + ' index ' + bt2output_bam + '.q30_filter.sorted.bam'
    flagstat_cmd = samtools_path + ' flagstat ' + bt2output_bam + '.q30_filter.sorted.bam > ' + flagstats
    idxstats_cmd = samtools_path + ' idxstats ' + bt2output_bam + '.q30_filter.sorted.bam > ' + idxstats

    # prepare input data for bowtie2
    try:
        os.system(gunzip_cmd)
    except Exception, e:
        print >> sys.stderr, 'Can\'t gunzip input data for mapping . . .'
        print >> sys.stderr, 'Exception: %s' % str(e)
        sys.exit(1)

    # run mapping with bowtie2
    try:
        os.system(bt2_cmd)
    except Exception, e:
        print >> sys.stderr, 'Can\'t run mapping command . . .'
        print >> sys.stderr, 'Exception: %s' % str(e)
        sys.exit(1)

    # filter, sort and index bam file
    try:
        os.system(min_map_qual_cmd)
        os.system(sort_bam_cmd)
        os.system(index_bam_cmd)
    except Exception, e:
        print >> sys.stderr, 'Can\'t run filter/sort/index bam file . . .'
        print >> sys.stderr, 'Exception: %s' % str(e)
        sys.exit(1)

    # get mapping stats
    print ('\n\n\t\t*Getting mapping stats . . .\n\n')
    try:
        os.system(flagstat_cmd)
        os.system(idxstats_cmd)
    except Exception, e:
        print >> sys.stderr, 'Can\'t get mapping stats . . .'
        print >> sys.stderr, 'Exception: %s' % str(e)
        sys.exit(1)
    print ('reference\tmapped reads\tunmapped reads')
    with open (idxstats, 'r') as fh:
        idxstats_file = fh.readlines()
        for line in idxstats_file:
            columns = line.split()
	    if re.match('\*', str(columns[0])) is None:
	        print ('{}\t{}\t{}' .format(columns[0], columns[2], columns[3]))
                # write results to file
		file_basename = os.path.basename(input_fastq_file)
                with open('{}/{}.___MAPPINGresults.csv' .format(scratch, file_basename), 'a') as fh:
                    fh.write('{}\t{}\t{}\t' .format(columns[0], columns[2], columns[3]))

    print ('\n\n\t*Finished mapping of RNA-seq data')
    print ("\t--- %s seconds ---\n\n" % (datetime.now() - start_time))
    os.remove(bt2input)
    #os.remove(bt2output_bam)

    return flagstats



################################################################################
### main script
if __name__ == "__main__":


    # read commandline arguments and create variables
    try:
        myopts, args = getopt.getopt(sys.argv[1:],'i:o:q:m:')
    except getopt.GetoptError as e:
        print (str(e))
        print ('usage: {} -i single fastq file (fastq.gz) -o directory for output files -q Y -m Y' .format(sys.argv[0]))
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
        if option == '-m':
            if argument == 'Y':
                run_MAP = True
            elif argument == 'N':
                run_MAP = False
            else:
                print >> sys.stderr, 'Please specify Y or N for the m option . . .'
                sys.exit(1)


    # print some messages --- could use this to check for multiple same files in the source input dir??
    border = '#################################################################################################'
    print ('{}{}\n## This is the test RNA-seq pipeline for the Vertis data\n{}{}' .format(border, border, border, border))
    print ('\n\n\tMore useful messages to follow here . . .\n\n')


    # set up working dir, results directory etc. and cd into scratch space
    try:
        results_sub_dir = pipeline_setup(border, scratch, input_fastq_file, results_dir)
    except Exception, e:
        print >> sys.stderr, 'Can\'t setup required directories . . .'
        print >> sys.stderr, 'Exception: %s' % str(e)
        sys.exit(1)


    # if QC wanted:
    if run_QC == True:
        try:
            trimmed_data = qc_data(border, input_fastq_file, results_sub_dir, scratch)
        except Exception, e:
            print >> sys.stderr, 'Can\'t run QC . . .'
            print >> sys.stderr, 'Exception: %s' % str(e)
            sys.exit(1)
    else:
        print ('\n\n\tThe QC option was not selected, skipping QC . . .\n\n')


    # if mapping wanted:
    if run_MAP == True:
        if run_QC == True:
            input_fastq_file = trimmed_data
        try:
            value = map_data(border, input_fastq_file, results_sub_dir, scratch)
        except Exception, e:
            print >> sys.stderr, 'Can\'t run mapping . . .'
            print >> sys.stderr, 'Exception: %s' % str(e)
            sys.exit(1)
    else:
        print ('\n\n\tThe mapping option was not selected, skipping mapping . . .\n\n')






### DELETE UNWANTED TEMP FILES - eg. bt2 input etc.


    # transfer files from scratch back to results dir (on brick) and delete scratch dir
    from distutils.dir_util import copy_tree
    copy_tree(scratch, results_dir)
    import shutil
    shutil.rmtree(scratch)
