#!/usr/bin/env python3
import os, sys, subprocess, re, getopt


"""
################################################################################
### INFO

This script maps the supplied FASTQ file using the supplied Bowtie2 reference.

version = 0.1
author  = Will Rowe
email = will.rowe@liverpool.ac.uk


################################################################################
### NOTES

This script is not complete:
    need to add customisation of parameters
    need to make stand alone (name==main)


"""


################################################################################
### SET PATHS AND DEFAULTS
# path to software on cluster
bt2_path = '/usr/local/share/bowtie2-2.2.5/bowtie2'
samtools_path = '/usr/local/bin/samtools'

# default settings
threads = '20'
MAPQ_score = '10'


################################################################################
### FUNCTIONS
"""
MAP data
this function runs Bowtie2 on the QCed data (or falls back to the raw fastq data if QC not run)
returns a filtered bam file (MAPQ > 10)
"""
def map_data(input_fastq_file, bt2_index, results_sub_dir, scratch, MAPPING_results_file):
    file_basename = os.path.basename(input_fastq_file)

    # commands
    bt2input = scratch + '/bt2.input.fq'
    bt2output_bam = results_sub_dir + '/bt2.output.bam'
    sorted_filtered_bam = results_sub_dir + '/bt2.output.bam.q' + MAPQ_score + '_filter.sorted.bam'
    flagstats = results_sub_dir + '/00_mapping.flagstats'
    idxstats = results_sub_dir + '/00_mapping.idxstats'
    gunzip_cmd = 'gzip -cd ' + input_fastq_file + ' > ' + bt2input
    bt2_cmd = bt2_path + ' -x ' + bt2_index + ' -q ' + bt2input + ' --very-sensitive-local -p ' + threads + ' | ' + samtools_path + ' view -bS - > ' + bt2output_bam
    sort_bam_cmd_1 = samtools_path + ' sort ' + bt2output_bam + ' ' + bt2output_bam + '.sorted'
    index_bam_cmd_1 = samtools_path + ' index ' + bt2output_bam + '.sorted.bam'
    mapped_reads_cmd = 'samtools view -F 0x904 -c ' + bt2output_bam + '.sorted.bam'
    min_map_qual_cmd = samtools_path + ' view -q ' + MAPQ_score + ' -b ' + bt2output_bam + ' > ' + bt2output_bam + '.q' + MAPQ_score + '_filter.bam'
    sort_bam_cmd = samtools_path + ' sort ' + bt2output_bam + '.q' + MAPQ_score + '_filter.bam ' + bt2output_bam + '.q' + MAPQ_score + '_filter.sorted'
    index_bam_cmd = samtools_path + ' index ' + bt2output_bam + '.q' + MAPQ_score + '_filter.sorted.bam'
    flagstat_cmd = samtools_path + ' flagstat ' + sorted_filtered_bam + ' > ' + flagstats
    idxstats_cmd = samtools_path + ' idxstats ' + sorted_filtered_bam + ' > ' + idxstats

    # prepare input data for bowtie2
    try:
        os.system(gunzip_cmd)
    except Exception, e:
        print >> sys.stderr, 'Can\'t gunzip input data for mapping . . .'
        print >> sys.stderr, 'Exception: %s' % str(e)
        sys.exit(1)

    # run mapping with bowtie2
    try:
        #os.system(bt2_cmd)
        p = subprocess.Popen(bt2_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        err = p.communicate()
    except Exception, e:
        print >> sys.stderr, 'Can\'t run mapping command . . .'
        print >> sys.stderr, 'Exception: %s' % str(e)
        sys.exit(1)

    # sort and index bam file, count mapped reads *** note that this only works for single end data
    try:
        os.system(sort_bam_cmd_1)
        os.system(index_bam_cmd_1)
	p = subprocess.Popen(mapped_reads_cmd, shell=True, stdout=subprocess.PIPE)
	init_mapped_reads, err = p.communicate()
    except Exception, e:
        print >> sys.stderr, 'Can\'t run sort/index original bam file . . .'
        print >> sys.stderr, 'Exception: %s' % str(e)
        sys.exit(1)

    # filter for MAPQ > MAPQ_score, sort and index filtered bam file
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
    print ('Reads aligned to the reference before filtering for MAPQ > 10:\n{}\n\nStats post-filtering:\nreference\tmapped reads\tunmapped reads' .format(init_mapped_reads))
    with open (idxstats, 'r') as fh:
        idxstats_file = fh.readlines()
        for line in idxstats_file:
            columns = line.split()
        # print mapped reads per reference to STDOUT
	    if re.match('\*', str(columns[0])) is None:
	        print ('{}\t{}\t{}' .format(columns[0], columns[2], columns[3]))
            # write mapped reads per reference to file
            with open(MAPPING_results_file, 'a') as fh:
                fh.write('{}\t{}\t{}\t' .format(columns[0], columns[2], columns[3]))
        # append the file name and the number of mapped reads (pre-mapping filter) to the results file
        with open(MAPPING_results_file, 'a') as fh:
            fh.write('{}\t{}\n' .format(file_basename, init_mapped_reads))

    # delete bowtie2 input (unzipped version of script input file)
    os.remove(bt2input)

    return sorted_filtered_bam


################################################################################
### MAIN SCRIPT
#if __name__ == "__main__":
