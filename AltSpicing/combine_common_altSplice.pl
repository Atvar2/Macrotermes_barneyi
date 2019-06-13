#! /usr/bin/perl  -w
use strict;

my $repeat1=shift;
my $repeat2=shift;
my $repeat3=shift;
my $type=shift;    #panbie shifou neihanzitiaoyue
my $outfile=shift;

my %sample1=();
my %sample2=();
my %sample3=();

my %pos=();

&read_file(\$repeat1,\%sample1,\%pos,$type) && print STDERR "read file $repeat1 completely\n";
&read_file(\$repeat2,\%sample2,\%pos,$type) && print STDERR "read file $repeat2 completely\n";
&read_file(\$repeat3,\%sample3,\%pos,$type) && print STDERR "read file $repeat3 completely\n";

open OUT, ">$outfile";
foreach my $pos (sort keys %pos){
	if(exists($sample1{$pos}) && exists($sample2{$pos}) && exists($sample3{$pos})){
		print OUT  "$sample1{$pos}\n";
	}
	else{
		next;
	}
}



sub read_file{
	my $file=shift;
	my $psam=shift;
	my $count=shift;
	my $t=shift;
	open IN, $$file;
	while (<IN>)	{
		chomp;
		next if ($_!~/^Mnat/);
		my @line=split /\t/;
		if ($t=~/Y/){
			if ($line[-1]=~/,/){
			my @pos=split (/,/,$line[-1]);
			foreach  my $last (@pos){
				my $pos="$line[0]\t$last";
				$psam->{$pos}= $_;
				$count->{$pos}++;
				}
			}
			else{
				my $pos="$line[0]\t$line[-1]";
                                $psam->{$pos}= $_;
                                $count->{$pos}++;
			}

		}
		else{
			my $pos="$line[0]\t$line[-2]";
			$psam->{$pos}= $_;
                        $count->{$pos}++;						
		}
	}
	close IN or die;
}
