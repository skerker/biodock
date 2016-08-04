#!/usr/bin/python
### Python script to meet CBI-1-C2 (+C1?) competencies
# CBI-1-C2: Take a DNA sequence and use standard bioinformatic tools to locate within a genome, annotate and infer function, including gene prediction, transcription factor (TF) analysis, splice-site boundaries, potential for copy number variants (CNVs).
# Script assumptions: query is a. cDNA, b. single fasta entry c. standard genetic code 

import os, sys, getopt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


######################################################################
# Housekeeping
######################################################################
temp_file_dir = './.TEMP_FILES'
if not os.path.exists(temp_file_dir):
	os.makedirs(temp_file_dir)
import time
date = time.strftime('%d-%b-%Y')


######################################################################
# Get arguments and load in query sequence
######################################################################
#store input and output file names
ifile = ''
ofile = ''

#read commandline arguments
try:
	myopts, args = getopt.getopt(sys.argv[1:],"i:o:")
except getopt.GetoptError as e:
	print (str(e))
	print ("usage: %s -i input -o output" % sys.argv[0])
	sys.exit(2)

#o == option, a == argument
for o, a in myopts:
	if o == '-i':
		ifile=a
	elif o == '-o':
		ofile=a

#print summary of arguments to terminal
print ('\n\nSEQUENCE ANALYSIS\n\ninput file:\t{}\noutput file:\t{}\n'  .format(ifile, ofile))

#read in fasta file and print summary of sequence to terminal
input_sequence = SeqIO.read(ifile, "fasta", IUPAC.unambiguous_dna)
print ('sequence header:\t{}\nsequence length:\t{} bases\n' .format(input_sequence.id, len(input_sequence)))


######################################################################
# Basic sequence analysis/manipulation and record keeping
######################################################################
#create SeqRecord to add data to
sequence_record = SeqRecord(input_sequence)
sequence_record.annotations['data_file_division'] = 'PLN'
sequence_record.annotations['date'] = date.upper()
sequence_record.id = "CBI-1"
sequence_record.name = "unknown"
sequence_record.description = "Sequence record for competency CBI-1."

#transcribe query sequence to mRNA
DNA = input_sequence.seq
sequence_record.seq = DNA
mRNA = DNA.transcribe()
sequence_record.mRNA = mRNA

#translate to protein
protein = mRNA.translate(table=1)
sequence_record.protein = protein


######################################################################
# Locate our query sequence within the human genome
######################################################################
#blast our query sequence against a local copy of RefSeqGene (human)
print ('BLASTing sequence against human RefSeq gene database . . .')
from Bio.Blast.Applications import NcbiblastnCommandline
blastn = NcbiblastnCommandline(query=ifile, db='/home/rowew/DATA/RefSeqGene_human/RefSeqGene', out='{}/query_to_refseqgene.xml' .format(temp_file_dir), outfmt=5)
stdout, stderr = blastn()
with open ('{}/query_to_refseqgene.xml' .format(temp_file_dir), 'r') as blastn_fh:
	from Bio.Blast import NCBIXML
	blast_record = NCBIXML.read(blastn_fh)

#grab the GI of the top blast hit and use it to look up the gene ID
hit_def = blast_record.alignments[0].hit_def
print ('top hit:\t{}\n\n' .format(hit_def))
import re
GI = re.match( r'gi\|(\d+)', hit_def, re.M|re.I)
GI = GI.group(1)
import urllib.request
elink = urllib.request.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=gene&id={}' .format(GI))

#use geneID to fetch the gene record for the top blast hit (putative homologue for our query sequence)
from bs4 import BeautifulSoup
soup = BeautifulSoup(elink, 'html.parser')
ID = soup.elinkresult.linksetdb.link.id
geneID = re.match( r'\<id\>(\d+)', str(ID), re.M|re.I)
geneID = geneID.group(1)
efetch = urllib.request.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={}&rettype=xml' .format(geneID))
soup = BeautifulSoup(efetch, 'html.parser')

#use the gene record to print the location information for our sequence
soup_ref_locus = soup.find('gene-ref_locus')
ref_locus = re.search( r'\>(\w+)\<', str(soup_ref_locus), re.M|re.I)
ref_locus = ref_locus.group(1)
soup_ref_maploc = soup.find('gene-ref_maploc')
ref_maploc = re.search( r'\>(\S+)\<', str(soup_ref_maploc), re.M|re.I)
ref_maploc = ref_maploc.group(1)
print ('gene locus information:\n\t{}\n\t{}\n\n' .format(ref_locus, ref_maploc))

#update sequence record for our sequence
sequence_record.name = ('{}' .format(ref_locus))



######################################################################
# SNVs
######################################################################
#blast our query sequence against a local copy of RefSeqRNA (human)
print ('BLASTing sequence against human RefSeq RNA database . . .')
from Bio.Blast.Applications import NcbiblastnCommandline
blastn = NcbiblastnCommandline(query=ifile, db='/home/rowew/DATA/RefSeqRNA_human/RefSeqRNA', out='{}/query_to_refseqrna.xml' .format(temp_file_dir), outfmt=5)
stdout, stderr = blastn()
with open ('{}/query_to_refseqrna.xml' .format(temp_file_dir), 'r') as blastn_fh:
        from Bio.Blast import NCBIXML
        blast_record = NCBIXML.read(blastn_fh)

#grab the GI of the top blast hit and download copy of the RefSeq RNA sequence
hit_def = blast_record.alignments[0].hit_def
print ('top hit:\t{}\n\n' .format(hit_def))
GI = re.match( r'gi\|(\d+)', hit_def, re.M|re.I)
GI = GI.group(1)
from Bio import Entrez
Entrez.email = 'wpmr2@cam.ac.uk'
handle = Entrez.efetch(db='nucleotide', id='{}' .format(GI), rettype='fasta')
refseq_sequence = SeqIO.read(handle, "fasta", IUPAC.unambiguous_dna)

#create a multifasta of our sequence + the RefSeq RNA sequence
sequences = []
sequences.append(refseq_sequence)
sequences.append(input_sequence)
with open ('{}/alignment_input.fa' .format(temp_file_dir), 'w') as fh:
        SeqIO.write(sequences ,fh, 'fasta')

#align the query to the reference
from Bio.Align.Applications import MuscleCommandline
muscle = MuscleCommandline(input='{}/alignment_input.fa' .format(temp_file_dir), out='{}/alignment.fa' .format(temp_file_dir))
stdout, stderr = muscle()

#view alignment and identify SNVs
from Bio import AlignIO
SNV = 0
SNV_counter = 0
alignment = AlignIO.read('{}/alignment.fa' .format(temp_file_dir), "fasta")
for i in range(0,len(alignment[1].seq)):
	if alignment[0,i] != alignment[1,i]:
		if alignment[0,i] != "-" and alignment[1,i] != "-":
			SNV = SNV + 1
			SNV_counter = SNV_counter + 1
			print ('SNV at position:', i, alignment[0,i], '>', alignment[1,i], SNV)
	else:
		SNV = 0
if SNV_counter == 0: print ('No SNVs found in query sequence (when compared to most similar RefSeq RNA (human)')


######################################################################
# Transcription factor analysis
######################################################################


######################################################################
# Splice site boundaries
######################################################################


######################################################################
# Copy number variants
######################################################################



######################################################################
# File outputs and tidy
######################################################################
#print sequence record to outfile
with  open( ofile, 'w') as f:
        SeqIO.write(sequence_record, f, "genbank")

#remove temporary files
import shutil
shutil.rmtree(temp_file_dir)


######################################################################
# Notes
######################################################################

#make BLAST result handling more robust - take in to account no results/ low scores
#repeat seq record work (since file loss) - add features and TF sites
