#! /usr/bin/perl  -w
use strict;

my $repeat1=shift;
my $repeat2=shift;
my $repeat3=shift;
my $repeat4=shift;
my $repeat5=shift;
#my $type=shift;    #panbie shifou neihanzitiaoyue
my $outfile=shift;

my %sample1=();
my %sample2=();
my %sample3=();
my %sample4=();
my %sample5=();

my %pos=();

&read_file(\$repeat1,\%sample1,\%pos) && print STDERR "read file $repeat1 completely\n";
&read_file(\$repeat2,\%sample2,\%pos) && print STDERR "read file $repeat2 completely\n";
&read_file(\$repeat3,\%sample3,\%pos) && print STDERR "read file $repeat3 completely\n";
&read_file(\$repeat4,\%sample4,\%pos) && print STDERR "read file $repeat3 completely\n";
&read_file(\$repeat5,\%sample5,\%pos) && print STDERR "read file $repeat3 completely\n";
my ($mps,$MPS,$mpw,$MPW,$N)=(1,2,3,4,5);
my ($a,$b,$c,$d,$e,
$ab,$ac,$ad,$ae,$bc,$bd,$be,$cd,$ce,$de,
$abc,$abd,$abe,$acd,$ace,$ade,$bcd,$bce,$bde,$cde,
$abcd,$acde,$bcde,$abcde)=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);

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
