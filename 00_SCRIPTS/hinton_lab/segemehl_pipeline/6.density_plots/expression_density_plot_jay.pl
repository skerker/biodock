use warnings;
use strict;
use Getopt::Long;

my %filter = ();
my $filter = '';
my $prefix = '';
my $x_range = '';
my $y_range = '0,0.2';

&GetOptions(
    'filter=s' => \$filter,
    'prefix=s' => \$prefix,
    'y_range=s' => \$y_range,
    'x_range=s' => \$x_range,
    );

if ($y_range ne '') {
    $y_range = ", ylim=c($y_range)";
}
if ($x_range ne '') {
    $x_range = ", xlim=c($x_range)";
}

if ($filter) {
    foreach (split /, ?/, $filter) { #/) {
	$filter{$_}++;
    }
}

my $header = <>;
chomp $header;
$header =~ s/\cM//;
my @headers = split /\t/, $header;

my %skip = ('min', 1,
	    'max', 1,
	    'mean', 1,
	    'gene', 1,
);

my %cat = ();
my $i = -1;
my @pos = ();
foreach (@headers) {
    $i++;
    if (/(.+) (\w+)$/) {

	next if (defined $skip{$1});

	if (%filter) {
	    next unless (defined $filter{$1});
	}

	$cat{$2}{$1} = $i;

	push @pos, $i;

    }
}

my %vals = ();
$i = -1;
while (<>) {
    $i++;
    chomp;
    my @h = split /\t/, $_;
    foreach (@pos) {
	my $val = $h[$_];
	if ($val <= 0) {
	    $val = 1;
	}
	$val = log($val) / log(2);
	$vals{$i}{$_} = $val;
    }
}
my $lines = $i;

foreach my $cat (sort keys %cat) {
    my $out = "$prefix$cat.pdf";
    my $R = "Rcmd.$prefix$cat";
    my $data = "$prefix$cat.input";
    open (DATA, ">$data")
	or die;
    my @samples = ();
    foreach my $i (0..$lines) {
	my @out = ();
	foreach my $sample (sort keys %{$cat{$cat}}) {
	    push @samples, $sample if ($i == 0);
	    my $col = $cat{$cat}{$sample};
	    push @out, $vals{$i}{$col};
	}
	print DATA "".(join "\t", @out)."\n";
    }
    close DATA;

    open (R, ">$R")
	or die;
    print R "
library(KernSmooth)
x = read.table('$data', sep='\\t')
attach(x);
fhat <- bkde(V1)
pdf('$out')
plot (fhat, xlab='log2 $cat', ylab='Density function', type='l', col=1, main='gene expression profile ($cat values)'$x_range$y_range)
";

    # vertical lines at 2 (1) and 10 (3.32)
    print R "abline(v=1, col='grey')\n";
    print R "abline(v=3.32, col='grey')\n";

    my $sample = shift @samples;
    my @legend = ($sample);
    my $colour = 1;
    my $line = 1;
    my @colours = ($colour);
    my @lines = ($line);

    if (@samples == 1) {
	$line = 2;
    }

    foreach my $sample (@samples) {
	push @legend, $sample;
	push @colours, $colour;
	push @lines, $line;
	$colour++;

	$line++ if ($colour % 8 == 0);
	print R "
fhat <- bkde(V$colour)    
points (fhat, type='l', lty = $line, col=$colour)
";
    } 
    
    print R "legend('topright', c('".(join "','", @legend)."'), lty = c(".(join ",", @lines)."), col = c(1 : ".(scalar @colours)."))\n";

    print R "dev.off()\n";
    close R;
    
    system "R CMD BATCH $R";
}
