#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my $max_len = 51;
my $trim_polya = '';
my $tag = '';
my $keep_n = '';
&GetOptions(
    'trim_polya' => \$trim_polya,
    'max_len=i' => \$max_len,
    'tag=s' => \$tag,
    'keep_n' => \$keep_n,
    );

my $tag_len = length($tag);

#@HWI-ST863:229:C20GRACXX:5:1101:2623:1983 1:N:0:ATCACGAT
#CTCGTCCGACGTCACCTAGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGAAAAAAACCCGCCCAAACCCCCCACCCAAACCCACCCCCTTTC
#+
#@CCFDFFFHHGDHGHIJGIJIJJIIIIJJJHFDDDDDDDDDDDDD#################################################
#@HWI-ST863:229:C20GRACXX:5:1101:2819:1953 1:N:0:ATCACGAT
#CTCGTCCGACGTCACCTAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#+

my $i = 0;
my %stats = ();
while (<>) {
    chomp;
    if (/^\@(.+)/
	or
	/^\@(.+)/) {
	$i++;
	my $header = $1;
	$header =~ s/ ef\:\d.+//;
	my $seq = <>;
	my $h2 = <>;
	my $qual = <>;

	if ($tag) {
	    if ($seq =~ /^$tag(.+)/) {
		$seq = $1;
		$stats{$tag}++;
	    } else {
		my $tag2 = substr $seq, 0, $tag_len, '';
		$stats{$tag}++;
	    }
	}

	# remove N's from beginning of sequence
	$seq =~ s/^N+// unless ($keep_n);

	if ($trim_polya) {
	    $seq =~ s/A+$//i;
	}

	my $len = length($seq);
	if ($len < $max_len) {
	    warn "only $len bp left for $header\n";
	    chomp $seq;
#	    next;
	} else {
	    $seq = substr $seq, 0, $max_len;
	}
	print ">$header\n$seq\n";
    } else {
	die "Wrong header symbol: $_";
    }
}

print STDERR "$i sequences read in\n";


foreach my $tag (sort { $stats{$a} <=> $stats{$b} } keys %stats) {
    print STDERR "$tag => $stats{$tag}\n";
}
