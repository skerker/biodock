#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

# input format like this:
#file name       Vertis  Pool    sample name     Barcode
#120824_SN132_A_L002_GQC-75-ACAGTG_R1.fastq.gz   3       1       D23_LSP ACAGTG
#120824_SN132_A_L002_GQC-75-ATCACG_R1.fastq.gz   3       1       D23_EEP ATCACG

my $max_len = 75;
my $prefix = '';
my $no_read_count = '';
my $restart = '';
&GetOptions(
    'max_len=i' => \$max_len,
    'no_read_count' => \$no_read_count,
    'restart' => \$restart,
    'prefix=s' => \$prefix,
    );

my $out_dir = '';

my @bams = (<*.bam>);

my %out = ();

while (<>) {
    my $reads = '';
    my $perc = '';
    next if (/^#/);
    chomp;
    my @F = split /\t/, $_;
    unless (-e $F[0]) {
	warn "Not found: $F[0]";
	next;
    } else {

	print "$_ (".(localtime)."):\n";

	if ($no_read_count) {
	    push @{$out{"$prefix$F[-2]"}}, '', '', '';
	} else {
	    if ($F[0] =~ /bz2$/) {
		$_ = `bunzip2 -c "$F[0]" | wc`;
	    } elsif ($F[0] =~ /gz$/) {
		$_ = `gunzip -c "$F[0]" | wc`;
	    } else {
		warn "wrong ending: '$F[0]'";
		next;
	    }

	    my @h = split;
	    $reads = $h[0] / 4;
	    print "# Reads: ".(&commify($reads))."\n";
	    push @{$out{"$prefix$F[-2]"}}, $reads, '', '';
	}

	unless (-s "$prefix$F[-2]_18.bam") {
#	    system "time bunzip2 -c $F[0] | perl -e 'while (<>) {if (/^\\\@HWI/) {print; \$_ = <>; print; }}' > $prefix$F[-2]_$max_len.fa";
#	    system "time perl -p -i -e 's/^\\\@HWI/>HWI/' $prefix$F[-2]_$max_len.fa; grep -c '>' $prefix$F[-2]_$max_len.fa";
	    if (-s "$prefix$F[-2]_$max_len.fa"
		or
		-s "$prefix$F[-2]_$max_len.fa.gz") {
		$reads = `cat \"reads.$F[-2].txt\"`;
                chomp $reads;
                $reads =~ s/ sequences read.+//;
	    } else {
		if ($F[0] =~ /bz2$/) {
		    system "time bunzip2 -c $F[0] | sian_fastq2fasta.pl -trim_polya -max_len $max_len -keep_n > \"$prefix$F[-2]_$max_len.fa\" 2> \"reads.$F[-2].txt\"";
		    $reads = `cat \"reads.$F[-2].txt\"`;
		    chomp $reads;
		    $reads =~ s/ sequences read.+//;
		} elsif ($F[0] =~ /gz$/) {
		    system "time gunzip -c $F[0] | sian_fastq2fasta.pl -trim_polya -max_len $max_len -keep_n > \"$prefix$F[-2]_$max_len.fa\" 2> \"reads.$F[-2].txt\"";
		    $reads = `tail -1 \"reads.$F[-2].txt\"`;
		    chomp $reads;
		    $reads =~ s/ sequences read.+//;
		} else {
		    warn "wrong ending: '$F[0]'";
		    next;
		}
	    }

	    unless (-s "$prefix$F[-2]_$max_len.fa") {
		if (-s "$prefix$F[-2]_$max_len.fa.gz") {
		    system "gunzip \"$prefix$F[-2]_$max_len.fa.gz\"";
		} else {
		    warn "No $prefix$F[-2]_$max_len.fa found - skipping...\n";
		    next;
		}
	    }
	    my $job = "length_pipeline.pl \"$prefix$F[-2]_$max_len.fa\" > \"segemehl_$F[-2].log\" 2> \"segemehl_$F[-2].err\"";
	    if ($restart) {
		$job .= ' -restart';
	    }
	    print STDERR "starting $job on ".(localtime)."\n";
	    system $job;
	    system "rm \"$prefix$F[-2]_$max_len.fa.gz\"" if (-e "$prefix$F[-2]_$max_len.fa.gz");
	    system "gzip -9 \"$prefix$F[-2]_$max_len.fa\"&";
	}

	# extract unique hits and combine all into one file
	my $out_sam = $out_dir."$prefix$F[-2]_segemehl_uniq.sam";
	my $mapped_total = 0;
	unless (-s $out_sam) {
	    my @bams = <*.bam>;
	    foreach (@bams) {
		next unless (-s $_);
		next unless (/^$prefix$F[-2]_(\d+).bam$/);
		print "extracting unique hits from $_\n";
		unless ($restart) {
		    my $stats = `samtools view "$_" | cut -d "\t" -f 1 | sort -u | wc`;
		    $stats =~ s/^\s+//;
		    my @h = split /\s+/, $stats;
		    $mapped_total += $h[0];
		}
		system "samtools view \"$_\" | grep \"NH:i:1\t\" >> \"$out_sam\"";
	    }
	}

	if ($reads) {
	    $perc = sprintf "%.1f", $mapped_total / $reads * 100;
	    push @{$out{"$prefix$F[-2]"}}, $mapped_total, $perc;

	    print "Total mapped for $prefix$F[-2]: ".(&commify($mapped_total))." ($perc%)\n";
	}
	# convert into sorted and indexed BAM file
	my $out_bam = $out_dir."$prefix$F[-2]_segemehl_uniq.bam";
	my $index = $out_dir."$prefix$F[-2]_segemehl_uniq.bam.sorted.bam.bai";
	unless (-s $index) {
	    print "converting $out_sam to $out_bam\n";
	    system "time samtools view -bt /DIR_MOUNT/SCRATCH/low-input-RNA-experiment/segemehl/ST4-74.fa.fai \"$out_sam\" > \"$out_bam\" 2> \"$out_bam.err\"";
	    print "Sorting and indexing $out_bam\n";
	    if (system "time samtools sort \"$out_bam\" \"$out_bam.sorted\"; time samtools index \"$out_bam.sorted.bam\"") {
		die "Problem with 'time samtools sort \"$out_bam\" \"$out_bam.sorted\"; time samtools index \"$out_bam.sorted.bam\"'";
	    }
	}



	unlink $out_bam;
	unlink $out_sam;

	my @stats = `samtools flagstat "$out_bam.sorted.bam"`;
	my $unique = '';
	foreach (@stats) {
	    if (/^([1-9]\d+) .+ mapped/) {
		$unique = $1;
		last;
	    }
	}

	if ($reads) {
	    $perc = sprintf "%.1f", $unique / $reads * 100;
	    my $perc2 = sprintf "%.1f", $unique / $mapped_total * 100;
	    print "Uniquely mapped for $prefix$F[-2]: ".(&commify($unique))." ($perc% of all, $perc2% of mapped reads)\n";
	    push @{$out{"$prefix$F[-2]"}}, $unique, $perc, $perc2;
	}
    }
}

foreach my $lib (sort keys %out) {
    my $out = join "\t", @{$out{$lib}};
    print "$lib\t$out\n";
}

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}
