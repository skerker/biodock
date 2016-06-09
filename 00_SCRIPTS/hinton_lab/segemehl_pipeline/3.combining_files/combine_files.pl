#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

my $add_header = '';
my $add_empty = '';
my $keep_header = '';
my $keep_order = '';
my $whole_row = '';
my $filter_only = '';

&GetOptions(
    'keep_order' => \$keep_order,
    'keep_header' => \$keep_header,
    'filter_only' => \$filter_only,
    'add_empty' => \$add_empty,
    'add_header' => \$add_header,
    'whole_row' => \$whole_row,
    );

my $ref = shift;

my %add = ();
my %fields = ();
my %names = ();
my %headers = ();
my @order = ();
foreach my $file (@ARGV) {
    my $name = $file;
    $name =~ s/\.txt//;
#    $name = 'V2_'.$name unless ($name =~ /^V\d_/);
    $names{$name}++;
    push @order, $name;
    open (IN, $file)
	or die "Can't read file '$file': $!";
    my $i = 0;
    while (<IN>) {
	chomp;
	my ($id, $val, @rest) = split /\t/, $_;
	if ($keep_header
	    and
	    $i == 0) {
	    $headers{$name} = join "\t", ($val, @rest);
	    $i++;
	    next;
	}
	$id =~ s/\"//g;
	if ($whole_row) {
	    $val = join "\t", ($val,@rest);
	}
	$add{$id}{$name} = $val;
	my @fields = split /\t/, $val;
	my $fields = scalar @fields;
	$fields{$name} = $fields;
    }
    close IN;
}

my @names = sort keys %names;
if ($keep_order) {
    @names = @order;
}

open (IN, $ref)
    or die "Can't open '$ref': $!";

my $header = <IN>;
chomp $header;
$header =~ s/\cM//;
if ($add_header) {
    close IN;
    open (IN, $ref);
    my $name = $ref;
    $name =~ s/\.txt//;
#    $name = 'V2_'.$name unless ($name =~ /^V\d_/);
    $header = "Gene name\t$name";
} 

foreach (@names) {
    if ($keep_header) {
	$header .= "\t$headers{$_}";
    } else {
	$header .= "\t$_" unless ($filter_only);
    }
}
print "$header\n";

while (<IN>) {
    chomp;
    s/\cM//;
    $_ =~ /^(.+?)\t/;
    my $id = $1;
    $id =~ s/\"//g;
    unless (defined $add{$id}) {
	if ($add_empty) {
	    foreach my $name (@names) {
		my @add = ();
		foreach (1..$fields{$name}) {
		    push @add, '';
		}
		$add{$id}{$name} = "".(join "\t", @add);
	    }
	} else {
	    warn "No $id - skipping\n";	
	    next;
	}
    }
    foreach my $name (@names) {
	my $val = '';
	if (defined $add{$id}{$name}) {
	    $val = $add{$id}{$name};
	}
	$_ .= "\t$val" unless ($filter_only);
    }
    print "$_\n";
}

