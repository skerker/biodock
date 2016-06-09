#!/usr/bin/perl
use warnings FATAL=>qw(all);
#use v5.10;
use File::Basename;
#use subs qw(safesys);
use Getopt::Long;

my $string="Now running @ARGV";
my $length=80<length $string ? 80 : length $string;
print STDERR join("\n", '-'x $length, $string, '-'x $length), "\n";
#system @_ and die "Error with command @_: $!";

my $restart = '';

&GetOptions(
    'restart' => \$restart,
    );


for (@ARGV) {
    open my $in, $_ or die "Can't open $_: $!";
    /_(\d+).fa/;
    my $base=basename $_, "_$1.fa";
    my $length=$1;
    while ($length>=18) {
        my $next=$length-3;
        print STDERR `date`;
	unless (-s "${base}_$length.bam") {
	    my $job = "segemehl.x -A 100 -i /DIR_MOUNT/SCRATCH/low-input-RNA-experiment/segemehl/ST474.idx -d /DIR_MOUNT/SCRATCH/low-input-RNA-experiment/segemehl/ST4-74.fa -q \"${base}_$length.fa\" -u \"$base.ump\" --threads 4 | samtools import /DIR_MOUNT/SCRATCH/low-input-RNA-experiment/segemehl/ST4-74.fa.fai - \"${base}_$length.bam\" 2> \"$base.$length.err\" > \"$base.$length.log\"";
	    print STDERR "$job\n";
#	    system $job and die "Error with command $job: $!";
	    system $job;
	    unless (-s "${base}_$length.bam"
		    or
		    -s "$base.ump") {
		die "Error with command $job: $!";
	    }
	}
	if ($restart) {
	    unless (-s "$base.ump") {
		if (-s "${base}_$next.fa") {
		    system "mv \"${base}_$next.fa\" \"$base.ump\"";
		} else {
		    warn "no \"${base}_$next.fa\"  - skipping";
		    next;
		}
	    }
	}
        open my $ump, "$base.ump" or die $!;
        if ($next>=18) {
            open my $trunc, ">${base}_$next.fa";
            while (<$ump>) {
                chomp;
		s/\ ef\:\d.+//;
                $_=substr $_, 0, $next if not /^>/;
                print $trunc "$_\n";
            }
        }
        $length=$next;
    }
}

sub safesys {
    $ENV{PATH} = '/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/local/software/biopieces/bp_bin:/opt/local/software/HomerTools/bin:/opt/local/software/ea-utils.1.1.2-537:/opt/local/software/ncbi-blast-2.2.30+/bin:/opt/local/software/idba-1.1.1/bin:/opt/local/software/MUMmer3.23:/opt/local/software/cd-hit:/opt/local/software/fastx_toolkit_0.0.13/bin:/opt/local/software/HomerTools/cpp/:/opt/local/software/HomerTools/.//bin/:/opt/local/software/biopieces/bp_bin:/opt/local/software/microbiomeutil-r20110519/ChimeraSlayer/:/opt/local/software/microbiomeutil-r20110519/AmosCmp16Spipeline:/opt/local/software/microbiomeutil-r20110519/NAST-iEr:/opt/local/software/microbiomeutil-r20110519/TreeChopper:/opt/local/software/microbiomeutil-r20110519/WigeoN:/opt/local/software/EMIRGE-master:/opt/local/software/Ray-2.3.1/ray-build:/MOUNTED_FILES/software/segemehl_0_0_9_3';
    delete @ENV{'IFS', 'CDPATH', 'ENV', 'BASH_ENV'};
    my $string="Now running @_";
    my $length=80<length $string ? 80 : length $string;
    print STDERR join("\n", '-'x $length, $string, '-'x $length), "\n";
    system @_ and die "Error with command @_: $!";
}
