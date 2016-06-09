#!/usr/bin/env sh
# This script compares the truncated BAM files generated by the Segemehl pipeline to a reference BAM file (eg. Bowtie2).
# It counts how many reads Segemehl mapped that the other mapper did not, at each Segemehl iteration.
# The script needs bedtools and fastq2fasta.pl in the path

for bam in /DIR_MOUNT/SCRATCH/different_mapping_parameters/segemehl/BAMS/*.bam
    do
        echo $bam
        filename=$(basename "$bam")
        bedtools bamtofastq -i $bam -fq ./$filename.fastq
        fastq2fasta.pl $filename.fastq > $filename.fa
        rm $filename.fastq
        perl -ne 'while(/(^>\S+)_/g){print "$1\n";}' $filename.fa > $filename.headers
        sort $filename.headers > $filename.sorted.headers
        diff $filename.sorted.headers bt2.q30.headers.sorted.txt > comparison.headers
        perl -ne 'while(/^<\s\S+/g){print "$&\n";}' comparison.headers > $filename.only.headers
        wc -l $filename.only.headers >> UNIQUE_COUNT.txt
        rm *.headers
    done
