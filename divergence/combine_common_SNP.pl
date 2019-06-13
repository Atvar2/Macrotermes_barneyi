#! /usr/bin/perl  -w
use strict;

my $repeat1=shift;
my $repeat2=shift;
my $repeat3=shift;
my $outfile=shift;

my %sample1=();
my %sample2=();
my %sample3=();

my %pos=();

&read_file(\$repeat1,\%sample1,\%pos) && print STDERR "read file $repeat1 completely\n";
&read_file(\$repeat2,\%sample2,\%pos) && print STDERR "read file $repeat2 completely\n";
&read_file(\$repeat3,\%sample3,\%pos) && print STDERR "read file $repeat3 completely\n";

open OUT, ">$outfile";
foreach my $pos (sort keys %pos){
	if(exists($sample1{$pos}) && exists($sample2{$pos}) && exists($sample3{$pos})){
		print OUT $pos{$pos},"\n";
	}
	else{
		next;
	}
}



sub read_file{
	my $file=shift;
	my $psam=shift;
	my $ppos=shift;
	open IN, $$file;
	while (<IN>)	{
		chomp;
		next if (/^\#/);
		my @line=split /\t/;
		my $pos="$line[0]\t$line[1]\t$line[3]\t$line[4]";
		$psam->{$pos}=$_;
		$ppos->{$pos}=$_;
	}
	close IN or die;
}
