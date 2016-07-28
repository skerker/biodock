#!/usr/bin/env python3
import os, sys, subprocess, re, getopt


"""
################################################################################
### INFO

This script QCs the supplied FASTQ file. It runs FastQC and Trimmomatic (with a sliding window, quality based trim), it then re-runs FastQC and inspects the reads with Kraken.

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
fastqc_path = '/usr/local/share/FastQC/fastqc'
kraken_path = '/pub46/willr/000_HOME/0005_RNA-SEQ-PIPELINE/01_BIN/kraken'
kraken_report_path = '/pub46/willr/000_HOME/0005_RNA-SEQ-PIPELINE/01_BIN/kraken-report'
kraken_db = '/pub46/willr/000_HOME/0002_REF_DATA/0004_DBs/minikraken_20141208'
trimm_path = '/pub46/willr/000_HOME/0005_RNA-SEQ-PIPELINE/01_BIN/trimmomatic-0.36.jar'

# default settings
threads = '20'

################################################################################
### FUNCTIONS
"""
QC data
this function runs FastQC and Trimmomatic on the raw fastq data, then runs Kraken and FastQC on the trimmed reads.
returns a QC results file and fastq.gz file of trimmed reads
"""
def qc_data(input_fastq_file, results_sub_dir, file_basename, QC_results_file):
    # save trimmed input data as
    trimmed_data = results_sub_dir + '/trimmed.' + file_basename

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
    with open(QC_results_file, 'w') as fh:
        fh.write('{}\t{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t{}\t{}\t\n' .format(file_basename, seq_length, base_qual, len_dist, dropped_seqs, seq_length_trimmed, base_qual_trimmed, len_dist_trimmed, unclassified_reads, seq_num))

    return trimmed_data


################################################################################
### MAIN SCRIPT
#if __name__ == "__main__":
