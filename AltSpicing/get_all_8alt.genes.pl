#! /usr/bin/perl  -w
use strict;
my $dir=shift;
$dir=~s/\/$//;

my $gene=(split(/\//,$dir))[-1];
my @files=glob "$dir/*.xls.C";
my %hash=();

foreach my $file (@files){
	chomp($file);
#	print "$file\n";
     my	$temp=(split (/\//,$file))[-1];
	my ($sample,$type)=(split (/\./,$temp))[0,1];
	open IN, $file;
	while (<IN>){
		chomp;
		my @line=split (/\t/,$_);
		my $alt="$sample\t$line[-3]\t$line[-2]\t$line[-1]";
		$hash{$type}{$alt}=$_;
	}
	close IN or die;
#	print "$sample\t$type\n";
}

open OUT, ">$gene.xls";
foreach my $type(keys %hash){ 
	print OUT "$type:\n";
	foreach my $pos (sort {$hash{$type}{$b} cmp  $hash{$type}{$a}} keys %{$hash{$type}}){
		my $sam =(split (/\t/,$pos))[0];
		print OUT "$hash{$type}{$pos}\t$sam\n";
	}
}
