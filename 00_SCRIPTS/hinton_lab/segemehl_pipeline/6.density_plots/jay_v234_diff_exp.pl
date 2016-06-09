use warnings;
use strict;
use Getopt::Long;


my $contrasts = '';
my $method = 'deseq';
my $annotation_file = '';
my $gene_length_file = '';
my $average_method = 'mean';
my $take_log = '';
my $rounding_ratio = '';
my $rounding_abs = '';
my $tcc_iteration = 1;
my $tcc_sample_size = 100;
my $tpm_vs_tcc = '';
my $tpm_comparator = 'TCC';
my $select_gene = '';
my $output_format = '';
my $read_length = 50;
my $universal_length = '';
my $remove_vertis = '';

&GetOptions(
    'remove_vertis' => \$remove_vertis,
    'read_length=i' => \$read_length,
    'universal_length=i' => \$universal_length,
    'method=s' => \$method,
    'gene=s' => \$select_gene,
    'tpm_vs_tcc' => \$tpm_vs_tcc,
    'tpm_comparator=s' => \$tpm_comparator,
    'tcc_iteration=i' => \$tcc_iteration,
    'tcc_sample_size=i' => \$tcc_sample_size,
    'rounding_ratio|rounding=s' => \$rounding_ratio,
    'rounding_abs=s' => \$rounding_abs,
    'take_log|log' => \$take_log,
    'average_method=s' => \$average_method,
    'annotation_file=s' => \$annotation_file,
    'gene_length_file=s' => \$gene_length_file,
    'contrasts=s' => \$contrasts,
    'output_format=s' => \$output_format,
    );

if ($method eq 'all') {
    $method = join ',', ('raw', 'RPKM', 'TPM', 'TPM_Wagner', 'TPM_SQRT', 'DESeq', 'edgeR', 'TCC', 'TCC_fast');
}

my %ann = ();
my @ann = ();
my $ann_size = 0;
if ($annotation_file) {
    open (IN, $annotation_file)
	or die;
    my $header = <IN>;
    $header =~ s/\cM//;
    chomp $header;
    @ann = split /\t/, $header;
    shift @ann;
    while (<IN>) {
	chomp;
	s/\cM//;
	my @h = split /\t/, $_;
	my $id = shift @h;
#	unless (@h) {
#	    @h = ('');
#	}
	$ann{$id} = join "\t", @h;
	if (@h > $ann_size) {
	    $ann_size = @h;
	}
    }
    close IN;
}

my %lengths = ();
if ($gene_length_file) {
    open (IN, $gene_length_file)
	or die;
    while (<IN>) {
	chomp;
	s/\cM//;
	my ($gene, $len, @rest) = split /\t/, $_;
	$lengths{$gene} = $len;
    }
    close IN;
}
if ($universal_length) {
    foreach (sort %lengths) {
	$lengths{$_} = $universal_length;
    }
}

my $header = <>;
$header =~ s/\cM//;
my %headers = &get_headers($header);

if ($contrasts eq 'all') {
    my @libs = ();
    foreach (sort keys %headers) {
	next unless (/^V\d/);
	if ($remove_vertis) {
	    s/V\d+[_ ]?//;
	}
	push @libs, $_;
    }
    my $lib = shift @libs;
    $contrasts = '';
    foreach (@libs) {
	$contrasts .= "$lib vs $_,";
    }
    chop $contrasts;
}

my @contrasts = ();
my %conditions = ();
my %libs = ();
my %all_libs = ();
my $duplicates = 0;
foreach my $contrast (split /, ?/, $contrasts) { #/) {
    push @contrasts, $contrast;

    my @conds = (split / vs /, $contrast);
    @{$conditions{$contrast}} = @conds;
    foreach (@conds) {
	%{$libs{$contrast}{$_}} = ();
    }

    foreach my $head1 (sort keys %headers) {
	my $head = $head1;
	if ($remove_vertis) {
            $head =~ s/V\d+[_ ]?//;
        }
	foreach my $cond (@{$conditions{$contrast}}) {
	    if ($head =~ /$cond$/
		or
		$head =~ /^$cond/) {
		if (defined $all_libs{$head}) {
		    if (defined $all_libs{$head}{$contrast}) {
			warn "Library $head fits for more than one condition in $contrast!\n";
		    }
		}
		warn "Using library $head for condition $cond in $contrast\n";
		$all_libs{$head}{$contrast} = $cond;
#		if (defined $libs{$contrast}{$cond}) {
#		    $duplicates++;
#		}
		$libs{$contrast}{$cond}{$head}++;
	    }
	}
    }
    foreach (@conds) {
	unless (%{$libs{$contrast}{$_}}) {
	    die "No library found for condition $_ in $contrast!";
	}
	my $libs = scalar keys %{$libs{$contrast}{$_}};
	if ($libs > 1) {
	    $duplicates = $libs if ($libs > $duplicates);
	}
    }
}

warn "Found up to $duplicates replicates!\n";

# extract raw data
my %in = ();
my %raw_lib_size = ();
while (<>) {
    chomp;
    s/\cM//;
    my @h = split /\t/, $_;
    my $id = $h[0];
    $id =~ s/\"//g;
    if ($select_gene) {
	next unless ($id eq $select_gene);
    }
    foreach my $lib (sort keys %all_libs) {
	my $col = $headers{$lib};
	my $val = $h[$col];
	$in{'raw'}{$id}{$lib} = $val;
	$raw_lib_size{$lib} += $val;
    }
}

# calculate factor for TPM_Wagner if required
my %skip = ();
my %T_W = ();
if ($method =~ /TPM_Wagner/) {
    foreach my $id (sort keys %{$in{'raw'}}) {
	unless (defined $lengths{$id}) {
	    warn "No gene length for $id!";
	    $skip{$id}++;
	    next;
	}
	my $len = $lengths{$id};
	
	foreach my $lib (sort keys %all_libs) {
	    my $val = $in{'raw'}{$id}{$lib};
	    $T_W{$lib} += $val / $len;
	}
    }
    foreach my $lib (sort keys %T_W) {
	warn "TWLIB\t$lib\t$T_W{$lib}\n";
    }
}

# calculate factor for TPM if required
my %T_L = ();
if ($method =~ /TPM/) {
    foreach my $id (sort keys %{$in{'raw'}}) {
	unless (defined $lengths{$id}) {
	    warn "No gene length for $id!\n";
	    $skip{$id}++;
	    next;
	}
	my $len = $lengths{$id};
	
	foreach my $lib (sort keys %all_libs) {
	    my $val = $in{'raw'}{$id}{$lib} / $raw_lib_size{$lib};
	    $T_L{$lib} += $val / $len;
	}
    }
    foreach my $lib (sort keys %T_L) {
	warn "TLIB\t$lib\t$T_L{$lib}\n";
    }
}


# carry out selected methods of normalisation
my %out = ();
foreach my $type (split /,/, $method) { #/) {
    &normalise($type);
}


# calculate TPM size factors for normalisation of raw libraries:
if (defined $in{'TPM'}) {
    my $tpm_factor_file = "tpm_factors_$read_length.txt";
    my $tpm_factor_file_wide = "tpm_factors_wide_$read_length.txt";
    open (TPM, ">$tpm_factor_file")
	or die;
    open (TPM_WIDE, ">$tpm_factor_file_wide")
	or die;
    print STDERR "Printing TPM factors to $tpm_factor_file and $tpm_factor_file_wide, using read length of $read_length\n";
    my @tpm = ();
    my @factor = ();
    foreach my $lib (sort keys %raw_lib_size) {
	my $raw = $raw_lib_size{$lib};
	my $tpm = 0;
	foreach my $gene (keys %{$in{'TPM'}}) {
	    my $val = $in{'TPM'}{$gene}{$lib};
	    my $len = $lengths{$gene};
	    $tpm += $val * $len;
	}
	$tpm /= $read_length;
	my $factor = $tpm / $raw;
	print TPM "$lib\t$factor\n";
	push @tpm, $lib;
	push @factor, $factor;
    }
    close TPM;
    print TPM_WIDE "".(join "\t", @tpm)."\n";
    print TPM_WIDE "".(join "\t", @factor)."\n";
    close TPM_WIDE;
}


# print output:
# 1. gene, annotation and length
# 2. all raw and derived values (rpkm, tpm, DESeq)
# 3. the ratios for each method
my @header = ('gene', @ann);
if ($gene_length_file) {
    push @header, 'gene length';
}
foreach ('raw', 'RPKM', 'TPM', 'TPM_Wagner', 'TPM_SQRT', 'DESeq', 'edgeR', 'TCC', 'TCC_fast') {
    my %done = ();
    if (defined $in{$_}) {
	foreach my $contrast (@contrasts) {
	    foreach my $cond (@{$conditions{$contrast}}) {
		foreach my $lib (sort keys %{$libs{$contrast}{$cond}}) {
		    next if (defined $done{$lib});
		    push @header, "$lib $_";
		    $done{$lib}++;
		}
	    }
	}
	push @header, "min $_\tmax $_\tmean $_";
    }
}

my @header2 = @header;
foreach my $type (split /,/, $method) { #/) {
    foreach my $contrast (@contrasts) {
	my $ratio = 'log2-ratio';
	unless ($take_log) {
	    $ratio = 'ratio';
	}
	push @header2, "$contrast $type $ratio";
    }
}
if ($tpm_vs_tcc) {
    foreach my $contrast (@contrasts) {
	push @header2, "$contrast TPM vs $tpm_comparator";
    }
}

foreach my $contrast (@contrasts) {
    my $ratio = 'log2-ratio';
    unless ($take_log) {
	$ratio = 'ratio';
    }
    foreach my $type (split /,/, $method) { #/) {	
	push @header, "$contrast $type $ratio";
    }

    if ($tpm_vs_tcc) {
	push @header, "$contrast TPM vs $tpm_comparator";
    }
}

if ($output_format) {
    print "".(join "\t", @header2)."\n";
} else {
    print "".(join "\t", @header)."\n";
}

my %tpm_vs_tcc = ();
foreach my $gene (sort keys %out) {
    next if (defined $skip{$gene});
    my @out = ($gene);
    if (%ann) {
	if (defined $ann{$gene}) {
	    push @out, $ann{$gene};
	} else {
	    my $ann = join "\t", ('' x $ann_size);
	    push @out, $ann;
	}	    
    }
    if (%lengths) {
	push @out, $lengths{$gene};
    }

    foreach my $type ('raw', 'RPKM', 'TPM', 'TPM_Wagner', 'TPM_SQRT', 'DESeq', 'edgeR', 'TCC', 'TCC_fast') {
	my %done = ();
	if (defined $in{$type}) {
	    my @vals = ();
	    foreach my $contrast (@contrasts) {
		foreach my $cond (@{$conditions{$contrast}}) {
		    foreach my $lib (sort keys %{$libs{$contrast}{$cond}}) {
			next if (defined $done{$lib});
			push @out, $in{$type}{$gene}{$lib};
			push @vals, $in{$type}{$gene}{$lib};
			$done{$lib}++;
		    }
		}
	    }
	    @vals = sort { $a <=> $b } @vals;
	    my $mean = sprintf "%.1f", &mean(@vals);
	    push @out, "$vals[0]\t$vals[-1]\t$mean";
	}
    }

    
    if ($output_format) {
	foreach my $type (split /,/, $method) { #/) {
	    foreach my $contrast (@contrasts) {
		push @out, $out{$gene}{$contrast}{$type};
	    }
	}
	if ($tpm_vs_tcc) {
	    foreach my $contrast (@contrasts) {
		my $val = '';
		if (defined $out{$gene}{$contrast}{'TPM'}
		    and
		    defined $out{$gene}{$contrast}{$tpm_comparator}
		    and
		    $out{$gene}{$contrast}{'TPM'} ne ''
		    and
		    $out{$gene}{$contrast}{$tpm_comparator} ne '') {
		    $val = $out{$gene}{$contrast}{'TPM'} - $out{$gene}{$contrast}{$tpm_comparator};
		}
		push @out, $val;
		$tpm_vs_tcc{$gene}{$contrast} = $val unless ($val eq '');
	    }
	}

    } else {

	foreach my $contrast (@contrasts) {
	    foreach my $type (split /,/, $method) { #/) {
		push @out, $out{$gene}{$contrast}{$type};
	    }
	    if ($tpm_vs_tcc) {
		my $val = '';
		if ($out{$gene}{$contrast}{'TPM'} ne ''
		    and
		    $out{$gene}{$contrast}{$tpm_comparator} ne '') {
		    $val = $out{$gene}{$contrast}{'TPM'} - $out{$gene}{$contrast}{$tpm_comparator};
		}
		push @out, $val;
		$tpm_vs_tcc{$gene}{$contrast} = $val unless ($val eq '');
	    }
	}
    }

    unless (defined $out[1]) {
	die;
    }
    print "".(join "\t", @out)."\n";
}

if (%tpm_vs_tcc) {
    my %skip = ();
    my $out_file = 'tpm_vs_tcc.txt';
    open (OUT, ">$out_file")
	or die;
    my @header = ('gene', @contrasts);
    print OUT "".(join "\t", @header)."\n";
    foreach my $gene (sort keys %tpm_vs_tcc) {
	my $contrasts = scalar keys %{$tpm_vs_tcc{$gene}};
	if ($contrasts < @contrasts) {
	    $skip{$contrasts}++;
	    next;
	}
	my @out = ($gene);
	foreach my $contrast (@contrasts) {
	    push @out, sprintf "%.2f", $tpm_vs_tcc{$gene}{$contrast};
	}
	print OUT "".(join "\t", @out)."\n";
    }
    close OUT;
    foreach my $key (sort keys %skip) {
	print STDERR "$key -> $skip{$key}\n";
    }
}
	    
sub normalise {
    my $method = shift;
    
    if ($method eq 'raw') {
	foreach my $gene (sort keys %{$in{'raw'}}) {
	    next if (defined $skip{$gene});
	    foreach my $c (@contrasts) {
		$out{$gene}{$c}{$method} = &get_ratio($gene, $method, $c);
	    }
	}
    } elsif ($method eq 'DESeq') {
	
	# write text file with raw counts
	# use R to carry out DESeq analysis
	# read results
	my $in = "vals_$method.txt";
	my $norm = "norm_$method.txt";
	my $sf_file = 'size_factor_DESeq.txt';
	unlink $sf_file;
	# set up headers for data file
	my @conds = ();
	my @libs = ();
	my %done = ();
	foreach my $c (@contrasts) {
	    foreach my $cond (@{$conditions{$c}}) {
		foreach my $lib (sort keys %{$libs{$c}{$cond}}) {
		    next if (defined $done{$lib});
		    push @libs, $lib;
		    push @conds, $cond;
		    $done{$lib}++;
		}
	    }
	}
	open (OUT, ">$in")
	    or die;
	print OUT "gene\t".(join "\t", @libs)."\n";
	foreach my $gene (sort keys %{$in{'raw'}}) {
	    next if (defined $skip{$gene});
	    my @out = ($gene);
	    foreach my $lib (@libs) {
		my $val = $in{'raw'}{$gene}{$lib};
		push @out, $val;
	    }
	    print OUT "".(join "\t", @out)."\n";
	}
	close OUT;

	my $conds = "".(join '", "', @conds);
	my $estimate = 'cds <- estimateDispersions( cds, method="blind", sharingMode="fit-only", fitType = "local" )';
	if ($duplicates) {
	    $estimate = 'cds <- estimateDispersions( cds )';
	}
	my $R = "library(DESeq)
count = read.table('$in' , header=TRUE, sep = \"\t\")
rownames(count) = count\$gene
conds <- factor( c( \"$conds\" ) )
cds <- newCountDataSet( count[,-1], conds )
cds <- estimateSizeFactors( cds )
$estimate
sf = sizeFactors( cds )
write.table(sf, '$sf_file', quote=FALSE, sep=\"\t\")
x = counts( cds, normalized=TRUE )
write.table(x, '$norm', quote=FALSE, sep=\"\t\")
";
	foreach my $c (@contrasts) {
	    my $c_deseq = $c;
	    unlink "$c.$method.txt";
#	    # switch conditions and add quotes
#	    $c_deseq =~ s/(.+) vs (.+)/$2", "$1/;
	    $c_deseq =~ s/ vs /", "/;
	    $R .= "res <- nbinomTest( cds, \"$c_deseq\" )
write.table( res, file=\"$c.$method.txt\", quote=FALSE, sep=\"\t\")
";
	}    

	my $R_cmd = "R.DESeq.cmd";
	open (R, ">$R_cmd")
	    or die;
	print R $R;
	close R;
	
	system "R CMD BATCH $R_cmd";
	
# read in normalised data
	open (IN, $norm)
	    or die "Can't read normalisation data '$norm'";
	my $head = <IN>;
	while (<IN>) {
	    chomp;
	    my @h = split /\t/, $_;
	    my $gene = shift @h;
	    foreach my $lib (@libs) {
		my $derived = shift @h;
		if ($rounding_abs) {
                    $derived = sprintf "%.${rounding_abs}f", $derived; 
                    $derived =~ s/\.0+$//; 
                } 
                $in{$method}{$gene}{$lib} = $derived; 
	    }
	}
	close IN;

	# read in results from contrasts
	foreach my $c (@contrasts) {
	    my $file = "$c.$method.txt";
	    open (IN, $file )
		or die "Can't read from $file!";
#id	baseMean	baseMeanA	baseMeanB	foldChange	log2FoldChange	pval	padj
#1	4.5S	117091.894796328	27012.5563377282	207171.233254928	7.66944196856978	2.93912161052288	0.81916825653591	1
	    my $header = <IN>;
	    my %headers = &get_headers($header);
	    my $col = $headers{'foldChange'};
	    if ($take_log) {
		$col = $headers{'log2FoldChange'};
	    }
	    while (<IN>) {
		chomp;
		my @h = split /\t/, $_;
		shift @h;
		my $gene = $h[0];
		my $ratio = $h[$col];
		if ($rounding_ratio
		    and
		    $ratio ne ''
		    and
		    $ratio ne 'NA') {
		    $ratio = sprintf "%.${rounding_ratio}f", $ratio;
		    $ratio =~ s/\.0+$//;
		}
		$out{$gene}{$c}{$method} = $ratio;
	    }
	    close IN;
	}
    } elsif ($method eq 'edgeR') {
	
	# write text file with raw counts
	# use R to carry out DESeq analysis
	# read results
	my %tmm = ();
	my $in = "vals_$method.txt";
	my $norm = "norm_$method.txt";
	my $sf_file = 'size_factor_edgeR.txt';
	unlink $sf_file;
	# set up headers for data file
	my @conds = ();
	my @libs = ();
	my %done = ();
	foreach my $c (@contrasts) {
	    foreach my $cond (@{$conditions{$c}}) {
		foreach my $lib (sort keys %{$libs{$c}{$cond}}) {
		    next if (defined $done{$lib});
		    push @libs, $lib;
		    push @conds, $cond;
		    $done{$lib}++;
		}
	    }
	}
	my $entries = 0;
	open (OUT, ">$in")
	    or die;
	print OUT "gene\t".(join "\t", @libs)."\n";
	foreach my $gene (sort keys %{$in{'raw'}}) {
	    next if (defined $skip{$gene});
	    my @out = ($gene);
	    $entries++;
	    foreach my $lib (@libs) {
		my $val = $in{'raw'}{$gene}{$lib};
		push @out, $val;
	    }
	    print OUT "".(join "\t", @out)."\n";
	}
	close OUT;

	my $conds = "".(join '", "', @conds);
	my $R = "library(edgeR)
x <- read.delim(\"$in\",row.names=\"gene\")
group <- factor(c( \"$conds\" ))
y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)
sf = y\$samples
write.table(sf, '$sf_file', quote=FALSE, sep=\"\t\")
z = cpm(y)
write.table(z, \"cpm.txt\", sep=\"\t\", quote=FALSE)
";

	if ($duplicates) {
	    $R .= "y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
";
	    if (@contrasts > 1) {		
		foreach my $c (@contrasts) {
		    my $c_deseq = $c;
		    unlink "$c.$method.txt";
		    # switch conditions and add quotes
		    $c_deseq =~ s/(.+) vs (.+)/$2", "$1/;
#	    $R .= "res <- nbinomTest( cds, \"$c_deseq\" )
#write.table( res, file=\"$c.$method.txt\", quote=FALSE, sep=\"\t\")
#";
		}
	    } else {
		my $c = $contrasts[0];
		unlink "$c.$method.txt";
		$R .= "tt = topTags(et, n=$entries, sort.by=\"none\")
write.table( tt, file=\"$c.$method.txt\", quote=FALSE, sep=\"\t\")
";
	    }
	}

	my $R_cmd = "R.edgeR.cmd";
	open (R, ">$R_cmd")
	    or die;
	print R $R;
	close R;
	
	system "R CMD BATCH $R_cmd";
	
# calculate in normalised data
	open (IN, $sf_file)
	    or die "Can't open $sf_file!";
#group   lib.size        norm.factors
#V2_ESP  1       6122475 1.6315083136883
#V2_LSP  2       2720478 0.612929760522844	    
	my $head = <IN>;
	while (<IN>) {
	    chomp;
	    my ($lib, $group, $size, $sf) = split /\t/, $_;
	    $tmm{$lib} = $sf;
	}
	close IN;

	foreach my $gene (sort keys %{$in{'raw'}}) {
	    next if (defined $skip{$gene});
	    unless (defined $lengths{$gene}) {
		warn "No length for $gene\n";
		$skip{$gene}++;
		next;
	    }
	    my $len = $lengths{$gene};

	    foreach my $lib (@libs) {
		my $factor = 1/$tmm{$lib};
		
		my $val = $in{'raw'}{$gene}{$lib};
		my $derived = $val / ($len/1000) / ($raw_lib_size{$lib}/$factor/1000000);
		if ($rounding_abs) {
                    $derived = sprintf "%.${rounding_abs}f", $derived; 
                    $derived =~ s/\.0+$//; 
                } 
                $in{$method}{$gene}{$lib} = $derived; 
	    }

	    # calculate contrasts:
	    foreach my $c (@contrasts) {
		$out{$gene}{$c}{$method} = &get_ratio($gene, $method, $c);
	    }
	}

    } elsif ($method eq 'TCC'
	     or
	     $method eq 'TCC_fast') {
	
	# write text file with raw counts
	# use R to carry out DESeq analysis
	# read results
	my %tmm = ();
	my $in = "vals_$method.txt";
	my $norm = "norm_$method.txt";
	my $sf_file = "size_factor_$method.txt";
	unlink $sf_file;
	# set up headers for data file
	my @conds = ();
	my @libs = ();
	my %done = ();
	foreach my $c (@contrasts) {
	    foreach my $cond (@{$conditions{$c}}) {
		foreach my $lib (sort keys %{$libs{$c}{$cond}}) {
		    next if (defined $done{$lib});
		    push @libs, $lib;
		    push @conds, $cond;
		    $done{$lib}++;
		}
	    }
	}
	open (OUT, ">$in")
	    or die;
	print OUT "gene\t".(join "\t", @libs)."\n";
	foreach my $gene (sort keys %{$in{'raw'}}) {
	    next if (defined $skip{$gene});
	    my @out = ($gene);
	    foreach my $lib (@libs) {
		my $val = $in{'raw'}{$gene}{$lib};
		push @out, $val;
	    }
	    print OUT "".(join "\t", @out)."\n";
	}
	close OUT;

	my $conds = "".(join '", "', @conds);
	my $R = "library(TCC)
x <- read.delim(\"$in\",row.names=\"gene\")
group <- factor(c( \"$conds\" ))
tcc <- new(\"TCC\", x, group)
";
	if ($duplicates) {
	    if ($method =~ /_fast/) {
		$R .= "tcc <- calcNormFactors(tcc, norm.method = \"tmm\", test.method = \"edger\", iteration = $tcc_iteration, FDR = 0.1, floorPDEG = 0.05)\n";
	    } else {
		$R .= "tcc <- calcNormFactors(tcc, norm.method = \"tmm\", test.method = \"bayseq\", iteration = 1, samplesize = $tcc_sample_size)\n";
	    }
	} else {
	    $R .= "tcc <- calcNormFactors(tcc, norm.method = \"deseq\", test.method = \"deseq\", iteration = 1, FDR = 0.1, floorPDEG = 0.05)\n";
	}
	$R .= "sf = tcc\$norm.factors
write.table(sf, '$sf_file', quote=FALSE, sep=\"\t\")
";
	$R .= "normalized.count <- getNormalizedData(tcc)
write.table(normalized.count, '$norm', quote=FALSE, sep=\"\t\")
";

	my $R_cmd = "R.$method.cmd";
	open (R, ">$R_cmd")
	    or die;
	print R $R;
	close R;
	
	system "R CMD BATCH $R_cmd";
	
	my $error = `tail -1 $R_cmd.Rout`;
	if ($error =~ / halted/) {
	    die "Error in R execution for $method ($R_cmd.Rout)!\n";
	}

# read in normalised data
	open (IN, $norm)
	    or die;
	my $head = <IN>;
	while (<IN>) {
	    chomp;
	    my @h = split /\t/, $_;
	    my $gene = shift @h;
	    foreach my $lib (@libs) {
		my $derived = shift @h;
		if ($rounding_abs) {
                    $derived = sprintf "%.${rounding_abs}f", $derived; 
                    $derived =~ s/\.0+$//; 
                } 
                $in{$method}{$gene}{$lib} = $derived; 
	    }
	}
	close IN;


	# generate ratios from normalised data	
	foreach my $gene (sort keys %{$in{$method}}) {
	    next if (defined $skip{$gene});
	    # calculate contrasts:
	    foreach my $c (@contrasts) {
		$out{$gene}{$c}{$method} = &get_ratio($gene, $method, $c);
	    }
	}

    } elsif ($method eq 'RPKM') {
	foreach my $gene (sort keys %{$in{'raw'}}) {
	    next if (defined $skip{$gene});
	    unless (defined $lengths{$gene}) {
		warn "No length for $gene\n";
		$skip{$gene}++;
		next;
	    }
	    my $len = $lengths{$gene};
	    foreach my $lib (sort keys %all_libs) {
		my $val = $in{'raw'}{$gene}{$lib};
		my $derived = $val / ($len/1000) / ($raw_lib_size{$lib}/1000000);
		if ($rounding_abs) {
		    $derived = sprintf "%.${rounding_abs}f", $derived;
		    $derived =~ s/\.0+$//;
		}
		$in{$method}{$gene}{$lib} = $derived;
	    }
	    
	    foreach my $c (@contrasts) {
		$out{$gene}{$c}{$method} = &get_ratio($gene, $method, $c);
	    }
	}
    } elsif ($method eq 'TPM_Wagner') {
	foreach my $gene (sort keys %{$in{'raw'}}) {
	    next if (defined $skip{$gene});
	    unless (defined $lengths{$gene}) {
		warn "No length for $gene\n";
                $skip{$gene}++;
                next;
	    }
	    my $len = $lengths{$gene};
	    foreach my $lib (sort keys %all_libs) {
		my $val = $in{'raw'}{$gene}{$lib};
		my $derived = $val * 1000000 / ($len * $T_W{$lib});
		if ($rounding_abs) {
		    $derived = sprintf "%.${rounding_abs}f", $derived;
		    $derived =~ s/\.0+$//;
		}
		$in{$method}{$gene}{$lib} = $derived;
	    }
	    
            foreach my $c (@contrasts) {
		$out{$gene}{$c}{$method} = &get_ratio($gene, $method, $c);
	    }
	}
    } elsif ($method eq 'TPM_SQRT') {
	foreach my $gene (sort keys %{$in{'raw'}}) {
	    next if (defined $skip{$gene});
	    unless (defined $lengths{$gene}) {
		warn "No length for $gene\n";
                $skip{$gene}++;
                next;
	    }
	    my $len = $lengths{$gene};
	    foreach my $lib (sort keys %all_libs) {
		my $val = $in{'raw'}{$gene}{$lib};
		my $derived = ($val / $raw_lib_size{$lib}) * 1000000 / ($len * $T_L{$lib});
		$derived = sqrt($derived);
		if ($rounding_abs) {
		    $derived = sprintf "%.${rounding_abs}f", $derived;
		    $derived =~ s/\.0+$//;
		}
		$in{$method}{$gene}{$lib} = $derived;
	    }

            foreach my $c (@contrasts) {
		$out{$gene}{$c}{$method} = &get_ratio($gene, $method, $c);
	    }
	}
    } elsif ($method eq 'TPM') {
	foreach my $gene (sort keys %{$in{'raw'}}) {
	    next if (defined $skip{$gene});
	    unless (defined $lengths{$gene}) {
		warn "No length for $gene\n";
                $skip{$gene}++;
                next;
	    }
	    my $len = $lengths{$gene};
	    foreach my $lib (sort keys %all_libs) {
		my $val = $in{'raw'}{$gene}{$lib};
		my $derived = ($val / $raw_lib_size{$lib}) * 1000000 / ($len * $T_L{$lib});
		if ($rounding_abs) {
		    $derived = sprintf "%.${rounding_abs}f", $derived;
		    $derived =~ s/\.0+$//;
		}
		$in{$method}{$gene}{$lib} = $derived;
	    }

            foreach my $c (@contrasts) {
		$out{$gene}{$c}{$method} = &get_ratio($gene, $method, $c);
	    }
	}
    }
}


sub average {
    unless (@_) {
	return '';
    }

    if ($average_method eq 'mean') {
	return &mean(@_);
    } else {
	return &median(@_);
    }
}

sub mean {
    my @values = @_;
    my $sum = 0;
    foreach (@values) {
	$sum += $_;
    }
    return $sum / @values;
}

sub median {

    unless (@_) {
	return '';
    }

    if (@_ == 1) {
	return $_[0];
    }

    my @order = sort { $a <=> $b } @_;
    my $md = '';

    if ($#order % 2 == 0){
        $md = $order[$#order / 2];
    } else {
	my $index1 = int($#order / 2);
	my $index2 = int($#order / 2) + 1;
	my $low = $order[$index1];
	my $high = $order[$index2];
	$md = $low + ($high - $low)/2;
    }
    return $md;
}


sub get_ratio {
    my $gene = shift;
    my $type = shift;
    my $contrast = shift;
    
    unless (defined $conditions{$contrast}) {
	die "$gene, $type, $contrast";
    }

    my %ratio = ();
    my @conds = @{$conditions{$contrast}};
    foreach my $cond (@conds) {
	my @vals = ();
	foreach my $lib (sort keys %{$libs{$contrast}{$cond}}) {
	    my $val = $in{$type}{$gene}{$lib};
	    push @vals, $val;
	}
	my $average = &average(@vals);
	$ratio{$cond} = $average;
    }
    my $ratio = '';
    if ($ratio{$conds[0]} ne ''
	and
	$ratio{$conds[0]} != 0) {
	if ($ratio{$conds[1]} ne '') {
	    $ratio = $ratio{$conds[1]} / $ratio{$conds[0]};
	}
    }
    if ($take_log) {
	if ($ratio ne ''
	    and
	    $ratio != 0) {
	    $ratio = log($ratio) / log(2);
	} else {
	    $ratio = '';
	}
    }
    if ($rounding_ratio
	and
	$ratio ne '') {
	$ratio = sprintf "%.${rounding_ratio}f", $ratio;
	$ratio =~ s/\.0+$//;
    }
    return $ratio;
}


sub get_headers {
    my %out = ();
    my $header = shift;
    chomp $header;
    my $count = 0;
    my @tmp = split /\t/, $header;
    foreach (@tmp) {
	$out{$_} = $count;
	$count++;
    }
    return %out;
}
