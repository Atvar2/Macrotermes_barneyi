#! /usr/bin/perl -w

use strict;
use File::Basename;

unless (@ARGV == 3) {
	print "\n\tperl $0 <Q20 dir> <fa_merge list><out dir>\n";
	exit (1);
}

my $snp_dir = $ARGV[0];
my $merge_list = $ARGV[1];
my $outdir=$ARGV[2];
`mkdir $outdir` unless (-d $outdir);
my %hList = ();	#{chr} => [start1, end1, start2, end2...]
my %hScaf = ();	#{chr}{start} => scaffold
my %hIndex = ();
my %chr=();
###read list and build hash### for iTools' result
open LIST, $merge_list or die "$!";
while (my $line = <LIST>) {
	chomp $line;
	next if ($line=~/^##/);
	my @info = split /\s+/, $line;
	push @{$hList{$info[3]}}, $info[4],$info[5];
#	push @{$chrscf{$info[3]}},$info[0];
	$hIndex{$info[3]}=1;
	$hScaf{$info[3]}{$info[4]} = $info[0];
}
close LIST;

###read snp files and trace back###
my @snp_files = `find $snp_dir -name "*.snp.vcf.xls"`;
chomp @snp_files;

foreach my $snp_file (@snp_files) {
	print "$snp_file\n";
        my $basename = basename $snp_file;
        my $sample_name = (split /\./, $basename)[0];
	##next unless (exists $hList{$chr_name});
	if ( exists $hList{$sample_name}){
#	  `cp $snp_file $outdir/$basename.new`;
	}else{
        open SNP, " < $snp_file" or die "$!";
        print "processing $sample_name\n";
        open OUT, ">$outdir/$sample_name.newsnp.vcf.xls" or die "$!";
	foreach my $chr_name(sort keys %hIndex){
	       $hIndex{$chr_name} = 0;
	}
        while (my $line = <SNP>) {
                chomp $line;
				if ($line =~ /^#/){print OUT "$line\n";next;}
                my @info = split /\t/, $line;
                my $chr = $info[0];
                my $pos = $info[1];
                for (my $index = $hIndex{$chr}; $index < (scalar @{$hList{$chr}})/2; ++$index) {   # how many scaffold in one chr
                        if ($pos < $hList{$chr}[2*$index]) {
                              $hIndex{$chr} = $index;
                                last;
                        }
                        elsif (($pos >= $hList{$chr}[2*$index]) and ($pos <= $hList{$chr}[2*$index+1])) {
                                $info[0] = $hScaf{$chr}{$hList{$chr}[2*$index]};
				$info[1] =$pos - $hList{$chr}[2*$index] + 1;
				my $newline = join("\t", @info);
				print OUT "$newline\n";
                                $hIndex{$chr} = $index;   
				#the SNP files need to be sorted. Then the circle will go on one by one.
				#in my opinion, you should not set this if your files are not in order.
                                last;
                        }
                        else {
                                next;   #always cicle until the position was in the range of chr. note the order
                        }
                } #end for
        } #end while
	}
        close OUT;
        close SNP;
}
