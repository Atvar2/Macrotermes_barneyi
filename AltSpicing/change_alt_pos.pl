#! /usr/bin/perl -w

use strict;
use File::Basename;

unless (@ARGV == 3) {
	&usage;
	exit (1);
}

my $alt_dir = $ARGV[0];
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
my @snp_files = `find $alt_dir -name "*.Novel_TU.xls"`;
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
        open OUT, ">$outdir/$sample_name.Novel_TU.xls" or die "$!";
	foreach my $chr_name(sort keys %hIndex){
	       $hIndex{$chr_name} = 0;
	}
#	my $title=<SNP>;
        while (my $line = <SNP>) {
                chomp $line;
		if ($line =~ /^Novel_TU_ID/){print OUT "$line\n";next;}
                my @info = split /\t/, $line;
                my $chr = $info[1];
                my @pos = split(/\,/,$info[-1]);
		$info[-1]="";
                for (my $index = 0; $index < (scalar @{$hList{$chr}})/2; ++$index) {   # how many scaffold in one chr
		  foreach my $pos (@pos){
                        if ($pos < $hList{$chr}[2*$index]) {
                              $hIndex{$chr} = $index;
                                next;
                        }
                        elsif (($pos >= $hList{$chr}[2*$index]) and ($pos <= $hList{$chr}[2*$index+1])) {
                                $info[1] = $hScaf{$chr}{$hList{$chr}[2*$index]};
				$info[-1] .=($pos - $hList{$chr}[2*$index] + 1).",";
#				my $newline = join("\t", @info);
#				print OUT "$newline\n";
                                $hIndex{$chr} = $index;   
				#the SNP files need to be sorted. Then the circle will go on one by one.
				#in my opinion, you should not set this if your files are not in order.
                                next;
                        }
                        else {
                                next;   #always cicle until the position was in the range of chr. note the order
                        }
		    } #end foreach
		} #end for
		    $info[-1]=~s/\,$//;
		   print OUT join("\t",@info),"\n";
        } #end while
	}
        close OUT;
        close SNP;
}

sub usage{
	print <<USAGE;
	usage:
	perl $0 <alt dir> <fa_merge list><out dir>;
	description:
	This software was used for transing chr pos to scf pos
	in transscript analysis of alt-splicing.
	detail:
	20150730 chenjunhu\@geneomics.cn
USAGE
}
	
	
