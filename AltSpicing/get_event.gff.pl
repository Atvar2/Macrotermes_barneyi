#! /usr/bin/perl  
use strict;
use Data::Dumper;
my $dir=shift;
my $type=shift;
#print " $dir\t$type\n";
my $i=0;my %hash=();
my @files=glob "$dir/*.$type.GeneDiffSpliceFilter.xls";
foreach my $file(@files){
#	print "$file\n";
	chomp($file);
	my $out;
	my $outfile=(split (/\/|\./,$file))[-4];
	$outfile .=".$type.gff";
	open OUT, ">$outfile";
	$out=\*OUT;
	if($type eq 'A3SS'){
		&read_A3SS($file,$out);
	}
	elsif ($type eq 'A5SS'){
		&read_A5SS($file,$out);
	}
	elsif($type eq 'RI'){
		&read_RI($file,$out);
	}
	else{
		&read_SE($file,$out);
	}
	close OUT;
}
print Dumper(%hash);
sub  read_A3SS{
	my $file=shift;
	my $out=shift;
	%hash=();
#	print "$file\n";
#	open OUT,">$out";
	$i=0;
	open IN, $file;
	while (<IN>){
		chomp;
		if(/^Mnat_/){
			my @t=split (/\t/,$_);
			my ($g_start,$g_end)=(0,0);
			if(exists($hash{$t[0]})){
				$hash{$t[0]}++;
			}
			else{
				$hash{$t[0]}=1;
				$i=0;
			}
			$t[0]="$t[0].$hash{$t[0]}";
			if ($t[2] eq '+'){
				$g_start=$t[7];
				$g_end=$t[4];
			}
			elsif($t[2] eq '-'){
				$g_start=$t[3];
                                $g_end=$t[8];
			}
			#	my $h=$$out;
				my $up_exon="$t[7]\t$t[8]";
				my $long_exon="$t[3]\t$t[4]";
				my $short_exon="$t[5]\t$t[6]";
				print $out "$t[1]\t$type\tgene\t$g_start\t$g_end\t.\t$t[2]\t.\tID=$t[0];Name=$t[0]\n";
				print  $out  "$t[1]\t$type\tmRNA\t$g_start\t$g_end\t.\t$t[2]\t.\tID=$t[0].A;Parent=$t[0]\n";
				print  $out "$t[1]\t$type\tmRNA\t$g_start\t$g_end\t.\t$t[2]\t.\tID=$t[0].B;Parent=$t[0]\n";
				print   $out  "$t[1]\t$type\texon\t$up_exon\t.\t$t[2]\t.\tID=$t[0].A.up;Parent=$t[0].A\n";
				print   $out  "$t[1]\t$type\texon\t$long_exon\t.\t$t[2]\t.\tID=$t[0].A.coreAndExt;Parent=$t[0].A\n";
				print   $out  "$t[1]\t$type\texon\t$up_exon\t.\t$t[2]\t.\tID=$t[0].B.up;Parent=$t[0].B\n";
				print  $out  "$t[1]\t$type\texon\t$short_exon\t.\t$t[2]\t.\tID=$t[0].B.core;Parent=$t[0].B\n";
		}
		next;
	}
	close IN or die;
}
sub  read_A5SS{
	my $file=shift;
	my $out=shift;
	$i=0;
	%hash=();
        open IN, $file;
        while (<IN>){
                chomp;
                if(/^Mnat_/){
                        my @t=split (/\t/,$_);
			if(exists($hash{$t[0]})){
                                $i++;
                        }
                        else{
                                $hash{$t[0]}=1;
				$i=1;
                        }
                        $t[0]="$t[0].$i"  unless ($i==0);
                        my ($g_start,$g_end)=(0,0);
                        if ($t[2] eq '+'){
                                $g_start=$t[3];
                                $g_end=$t[8];
                        }
                        elsif($t[2] eq '-'){
                                $g_start=$t[7];
                                $g_end=$t[4];
                        }

                                my $up_exon="$t[7]\t$t[8]";
                                my $long_exon="$t[3]\t$t[4]";
                                my $short_exon="$t[5]\t$t[6]";
                                print $out "$t[1]\t$type\tgene\t$g_start\t$g_end\t.\t$t[2]\t.\tID=$t[0];Name=$t[0]\n";
                                print $out "$t[1]\t$type\tmRNA\t$g_start\t$g_end\t.\t$t[2]\t.\tID=$t[0].A;Parent=$t[0]\n";
                                print $out  "$t[1]\t$type\tmRNA\t$g_start\t$g_end\t.\t$t[2]\t.\tID=$t[0].B;Parent=$t[0]\n";
				print $out  "$t[1]\t$type\texon\t$long_exon\t.\t$t[2]\t.\tID=$t[0].A.coreAndExt;Parent=$t[0].A\n";
                                print $out  "$t[1]\t$type\texon\t$up_exon\t.\t$t[2]\t.\tID=$t[0].A.dn;Parent=$t[0].A\n";
				print $out  "$t[1]\t$type\texon\t$short_exon\t.\t$t[2]\t.\tID=$t[0].B.core;Parent=$t[0].B\n";
                                print $out  "$t[1]\t$type\texon\t$up_exon\t.\t$t[2]\t.\tID=$t[0].B.dn;Parent=$t[0].B\n";
                }
		next;
        }
        close IN or die;
}		
sub  read_RI{
	my $file=shift;
	my $out=shift;
	%hash=();
	$i=0;
        open IN, $file;
        while (<IN>){
		chomp;
		 if(/^Mnat_/){
                        my @t=split (/\t/,$_);
			if(exists($hash{$t[0]})){
                                $i++;
                        }
                        else{
                                $hash{$t[0]}=1;
				$i=1;
                        }
                        $t[0]="$t[0].$i"  unless ($i==0);
			my ($g_start,$g_end)=(0,0);
				$g_start=$t[3];
				 $g_end=$t[4];
                                my $up_exon="$t[5]\t$t[6]";
                                my $IR_exon="$t[3]\t$t[4]";
                                my $dn_exon="$t[7]\t$t[8]";
                                print $out "$t[1]\t$type\tgene\t$g_start\t$g_end\t.\t$t[2]\t.\tID=$t[0];Name=$t[0]\n";
                                print $out  "$t[1]\t$type\tmRNA\t$g_start\t$g_end\t.\t$t[2]\t.\tID=$t[0].A;Parent=$t[0]\n";
                                print $out  "$t[1]\t$type\tmRNA\t$g_start\t$g_end\t.\t$t[2]\t.\tID=$t[0].B;Parent=$t[0]\n";
                                print $out "$t[1]\t$type\texon\t$IR_exon\t.\t$t[2]\t.\tID=$t[0].A.withRI;Parent=$t[0].A\n";
                                print $out  "$t[1]\t$type\texon\t$up_exon\t.\t$t[2]\t.\tID=$t[0].B.up;Parent=$t[0].B\n";
                                print $out "$t[1]\t$type\texon\t$dn_exon\t.\t$t[2]\t.\tID=$t[0].B.dn;Parent=$t[0].B\n";
		}
	}
	close IN or die;
}
sub read_SE{
	my $file=shift;
	my $out=shift;
	%hash=();
	$i=0;
        open IN, $file;
        while (<IN>){
                chomp;
                 if(/^Mnat_/){
			 my @t=split (/\t/,$_);
			 if(exists($hash{$t[0]})){
				
                                $i++;
                        }
                        else{
				$i=1;
                                $hash{$t[0]}=1;
                        }
                        $t[0]="$t[0].$i"  unless ($i==0);
                         my ($g_start,$g_end)=(0,0);
                         $g_start=$t[5];
                         $g_end=$t[8];
		         my $up_exon="$t[5]\t$t[6]";
			 my $SE_exon="$t[3]\t$t[4]";
			 my $dn_exon="$t[7]\t$t[8]";
			print $out "$t[1]\t$type\tgene\t$g_start\t$g_end\t.\t$t[2]\t.\tID=$t[0];Name=$t[0]\n";
			print $out  "$t[1]\t$type\tmRNA\t$g_start\t$g_end\t.\t$t[2]\t.\tID=$t[0].A;Parent=$t[0]\n";
			print $out "$t[1]\t$type\tmRNA\t$g_start\t$g_end\t.\t$t[2]\t.\tID=$t[0].B;Parent=$t[0]\n";
			print $out "$t[1]\t$type\texon\t$up_exon\t.\t$t[2]\t.\tID=$t[0].A.up;Parent=$t[0].A\n";
			print $out "$t[1]\t$type\texon\t$SE_exon\t.\t$t[2]\t.\tID=$t[0].A.SE;Parent=$t[0].A\n";
			print $out "$t[1]\t$type\texon\t$dn_exon\t.\t$t[2]\t.\tID=$t[0].A.dn;Parent=$t[0].A\n";
			print $out  "$t[1]\t$type\texon\t$up_exon\t.\t$t[2]\t.\tID=$t[0].B.up;Parent=$t[0].B\n";
			print $out "$t[1]\t$type\texon\t$dn_exon\t.\t$t[2]\t.\tID=$t[0].B.dn;Parent=$t[0].B\n";
		}
	}
	close IN or die;
}


