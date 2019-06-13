#! /usr/bin/perl  -w
use strict;
 
if (@ARGV != 2){
	&usage;
	die "You should provide two parameters: position and genome_fasta!!!$!\n";
}
my %hash=();
my $pos=shift;
my $fa=shift; 
&read_pos($pos);
open IN, $fa;
$/=">"; <IN>; $/="\n";
while (<IN>){
  chomp;
  my  $output1;my $output2;
  my $chr=$1 if(/^(\S+)/);
  $/=">";
  my $seq=<IN>;
  chomp $seq;
  $seq=~s/\s+|\>//g;
  $/="\n";
  if (exists($hash{$chr})){
	foreach my $gene (keys %{$hash{$chr}}){
		my $strand=$hash{$chr}{$gene}{'strand'};
		my @iso=@{$hash{$chr}{$gene}{iso2}};
		my $iso;my $iso1;
		if (exists($hash{$chr}{$gene}{iso1})){
			my @iso1=@{$hash{$chr}{$gene}{iso1}};
			for (my $i=0; $i <@iso1; $i++){
                           $iso1 .= substr($seq,$iso1[$i][0]-1, $iso1[$i][1] - $iso1[$i][0] + 1);
			}
			my $mark="$gene [IC] locus=$chr:$iso1[0][0]:$iso1[-1][1]:$strand";
			$iso1 = Complement_Reverse($iso1) if($strand eq '-');
                        Display_seq(\$iso1);
			print ">$mark\n$iso1" if ($iso1=~/[ATCGN]/ig);
                }

		for (my $i=0; $i <@iso; $i++){
		           $iso .= substr($seq,$iso[$i][0]-1, $iso[$i][1] - $iso[$i][0] + 1);
		}	
		my $mark="$gene [IS] locus=$chr:$iso[0][0]:$iso[-1][1]:$strand";
		$iso = Complement_Reverse($iso) if($strand eq '-');
                        Display_seq(\$iso);
		print ">$mark\n$iso";
		}

	}
}

sub read_pos{
        my $pos=shift;
        open F2, $pos;
        while (<F2>){
                chomp;
                s/\r//g;
                my @t=split (/\s+/,$_);
                next if (/^\s+/);
                if($t[0] eq 'MXE' || $t[0] eq 'A3SS' || $t[0] eq 'A5SS'){
                        $hash{$t[2]}{$t[1]}{'strand'}=$t[3];
                        for(my $i=4;$i < scalar(@t) ;$i+=2){
                                if ($i==4){
                                        push @{$hash{$t[2]}{$t[1]}{'iso1'}},[$t[$i],$t[$i+1]];
                                        }
                                        if ($i==6){
                                                push @{$hash{$t[2]}{$t[1]}{'iso2'}},[$t[$i],$t[$i+1]];
                                        }
                                        if ($i>=8){
                                                 push @{$hash{$t[2]}{$t[1]}{'iso1'}},[$t[$i],$t[$i+1]];
                                                 push @{$hash{$t[2]}{$t[1]}{'iso2'}},[$t[$i],$t[$i+1]];
                                        }
                                }
                        }
                elsif($t[0] eq 'RI' ){
                         $hash{$t[2]}{$t[1]}{'strand'}=$t[3];
                         for(my $i=4;$i < scalar(@t) ;$i+=2){
                                        push @{$hash{$t[2]}{$t[1]}{'iso2'}},[$t[$i],$t[$i+1]];
                	}
		}
		elsif($t[0] eq 'SE'){
			$hash{$t[2]}{$t[1]}{'strand'}=$t[3];
                         for(my $i=6;$i < scalar(@t) ;$i+=2){
				push @{$hash{$t[2]}{$t[1]}{'iso2'}},[$t[$i],$t[$i+1]];
			}
        	}
	}
        close  F2  or die;
 foreach my $scf (keys %hash){
                foreach my $gene (keys %{$hash{$scf}}){
                        my $gene_p=$hash{$scf}{$gene};
                        my @iso1 = sort {$a->[0] <=> $b->[0]} @{$gene_p->{'iso1'}};
                        my @iso2 = sort {$a->[0] <=> $b->[0]} @{$gene_p->{'iso2'}};
                        $gene_p->{iso1} = \@iso1;
                        $gene_p->{iso2} = \@iso2;
                }
        }
}

sub usage{
	print <<USAGE;
	This program was used for get alt sequence basing on rmats result.
	The position file should be as the format:
        SE      Mnat_00428      scaffold9       -       494141  495035  491740  491849  504246  504428
	=============================================================================================	
	usage: perl  $0 position  genome_fasta
	=============================================================================================
USAGE
}

sub Complement_Reverse{
        my $seq=shift;
        $seq=~tr/AGCTagct/TCGAtcga/;
        $seq=reverse($seq);
        return $seq;

}
sub Display_seq{
        my $seq_p=shift;
        my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
        my $disp;

        $$seq_p =~ s/\s//g;
        for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
                $disp .= substr($$seq_p,$i,$num_line)."\n";
        }
        $$seq_p = ($disp) ?  $disp : "\n";
}


