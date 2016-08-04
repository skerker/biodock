#!/usr/bin/python
# this script takes a cDNA sequence, translates the peptides in each open reading frame (ORF), finds the ORF with the longest uninterupted peptide, and then blasts this to the ncbi protein database in order to give a putative annotation

## CHANGES TO MAKE: the reverse complement function (A) is not needed as the program is for cDNA (single stranded), not DNA sequences. i.e. only 3 ORFs on a single strand need to be checked.

import sys, getopt, re

################################################################################
## Minimum protein length cut off ( number of amino acids in peptide chain):
min_protein_length = 100

## BLAST settings:
E_VALUE_THRESH = 0.00001
IDENTITY_THRESH = 100

################################################################################
# FUNCTIONS

## A. reverse complement (unnecessary as only 3 ORF on cDNA are required)
def reverse_complement(dna):
    lookup = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    try:
        return ''.join([lookup[base] for base in reversed(dna)])
    except KeyError:
        print ('could not perform reverse complement')
        sys.exit(2)

## B. amino acid lookup
def codon_lookup(triplet):
    protein = None
    if len(triplet) == 3 and triplet in CODON_TABLE_1:
        amino_acid = CODON_TABLE_1[triplet]
    if len(triplet) != 3:
        print ('could not perform translation: encountered a codon in an ORF that was not 3 bases long ({})' .format(triplet))
        sys.exit(2)
    return amino_acid

## C. translation of open reading frame (ORF): cDNA sequence -> amino acid string (from the first start codon to the first stop codon encountered)
def translation(orf):
    aa_string = ''
    protein_length = 0
    seq_length = len(orf)
    ## from the input ORF, take 3 bases at a time and find if they code for a start codon.
    ## use the step argument of the range function to remain in the ORF (progress 3 bases on each iteration)
    for i in range(0, seq_length, 3):
        triplet_query = orf[i:i+3]
        if len(triplet_query) != 3:
            break
        amino_acid = codon_lookup(triplet_query)
        if amino_acid and amino_acid == 'M':
            ## translate all codons into amino acids and push to aa_string until stop codon is reached
            found_stop = False
            for triplet in range(i, seq_length, 3):
                aa = codon_lookup(orf[triplet:triplet+3])
                if not aa:
                    break
                if aa == 'Stop':
                    found_stop = True
                    break
                aa_string += aa
            ## if aa_string is greater than X aminoacids, return protein sequence
            if found_stop and len(aa_string) >= min_protein_length:
                protein_length = len(aa_string)
                break
            ## else, return flag and change protein length to 0 so that the checking step (line 96) halts the program if no proteins are found/stringency too high
            if found_stop and len(aa_string) < min_protein_length:
                protein_length = len(aa_string)
                aa_string = 'translated peptide chain shorter than ' + str(min_protein_length) + ' amino acids\n'
                protein_length = 0
                break
    if protein_length == 0:
        aa_string = 'not found'
    return (aa_string, protein_length)

## D. translation for each ORF (using functions A, B and C).
## 6 ORFs not actually necessary as this script is for single stranded cDNA, not DNA sequences
def protein_finder(cDNA, output):
    orf1 = cDNA.upper()
    orf2 = cDNA.upper()[1:]
    orf3 = cDNA.upper()[2:]
    orf4 = reverse_complement(cDNA.upper())
    orf5 = orf4[1:]
    orf6 = orf5[2:]
    orfs = (orf1, orf2, orf3, orf4, orf5, orf6)
    protein_seqs = ''
    protein_lengths = []
    counter = 1
    for orf in orfs:
        protein_seq, protein_length = translation(orf)
        protein_seqs += ('\n>fasta entry: ORF {}\n' .format(counter))
        protein_seqs += protein_seq
        counter += 1
        protein_lengths.append(protein_length)
    ## perform some checks
    if sum(protein_lengths) == 0:
        print ('No proteins were translated: could try reducing the min_protein_length value')
        sys.exit(2)
    print ('The longest protein was translated from ORF {}\n\nNow comparing the longest protein to the NCBI nr protein database using remote BLAST . . .\n' .format(protein_lengths.index(max(protein_lengths))+1))
    with open (output, 'w') as fh:
        fh.write(protein_seqs)
        fh.close()
    return protein_seqs, protein_lengths

## E. BLAST longest translated protein against nr protein database, use top hit to annotate our protein. Using BioPython here to make blast search easier.
def protein_blast_search(output, protein_lengths):
    from Bio import SeqIO
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
    ## reopen output file (multifasta of protein seqs ---> is this the best way?
    sequences = [seq_record.seq for seq_record in SeqIO.parse(output, "fasta")]
    protein_query = sequences[protein_lengths.index(max(protein_lengths))]
    ## blast with qblast(blastp), use SearchIO to read and write the results, then find if there is a (top) hit.
    result_handle = NCBIWWW.qblast("blastp", "nr", protein_query, expect= E_VALUE_THRESH, i_thresh= IDENTITY_THRESH)
    from Bio import SearchIO
    record = SearchIO.read(result_handle, 'blast-xml')
    #SearchIO.write(record, 'record.xml', 'blast-xml') # debug - file write not necessary
    try:
        ref_accession, ref_description, ref_length = (record[0].accession), (record[0].description), (record[0].seq_len)
    except IndexError:
        print ('\nBLAST failed - try different stringency settings\n\n')
        sys.exit(2)
    return (ref_accession, ref_description, ref_length)

# codon table 1 from NCBI
CODON_TABLE_1 = {
    'TTT': 'F',     'CTT': 'L',     'ATT': 'I',     'GTT': 'V',
    'TTC': 'F',     'CTC': 'L',     'ATC': 'I',     'GTC': 'V',
    'TTA': 'L',     'CTA': 'L',     'ATA': 'I',     'GTA': 'V',
    'TTG': 'L',     'CTG': 'L',     'ATG': 'M',     'GTG': 'V',
    'TCT': 'S',     'CCT': 'P',     'ACT': 'T',     'GCT': 'A',
    'TCC': 'S',     'CCC': 'P',     'ACC': 'T',     'GCC': 'A',
    'TCA': 'S',     'CCA': 'P',     'ACA': 'T',     'GCA': 'A',
    'TCG': 'S',     'CCG': 'P',     'ACG': 'T',     'GCG': 'A',
    'TAT': 'Y',     'CAT': 'H',     'AAT': 'N',     'GAT': 'D',
    'TAC': 'Y',     'CAC': 'H',     'AAC': 'N',     'GAC': 'D',
    'TAA': 'Stop',  'CAA': 'Q',     'AAA': 'K',     'GAA': 'E',
    'TAG': 'Stop',  'CAG': 'Q',     'AAG': 'K',     'GAG': 'E',
    'TGT': 'C',     'CGT': 'R',     'AGT': 'S',     'GGT': 'G',
    'TGC': 'C',     'CGC': 'R',     'AGC': 'S',     'GGC': 'G',
    'TGA': 'Stop',  'CGA': 'R',     'AGA': 'R',     'GGA': 'G',
    'TGG': 'W',     'CGG': 'R',     'AGG': 'R',     'GGG': 'G'
}


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
    print ('\nThis script takes a single DNA sequence and produces an output file containing the protein sequences found in each reading frame.\ninput file:\t{}\noutput file:\t{}\n'  .format(input, output))

## open input, extract cDNA sequence from file and store as string
    with open (input, 'r') as fh:
        input_file = fh.readlines()
        fh.close()
    sequence = ''
    for i in input_file:
        if re.match('\>', i):
            sequence_header = i
        else:
            sequence += i.replace('\n','')

## find proteins in each ORF, temporarily write proteins to outfile
    protein_seqs, protein_lengths = protein_finder(sequence, output)

## send translated protein to protein blast search function and return the top hit accession, description and reference length --> save this as putative annotation of the protein encoded by our cDNA
    ref_accession, ref_description, ref_length = protein_blast_search(output, protein_lengths)
    print ('{}\n{}\n{}\n\n' .format(ref_accession, ref_description, ref_length))
    with open (output, 'w') as fh:
        fh.write('Top hit for translated protein:\nAccession:\t{}\nDescription:\t{}\nLength:\t{}' .format(ref_accession, ref_description, ref_length))
        fh.close()