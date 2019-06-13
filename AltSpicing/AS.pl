#!/usr/bin/perl

=head1 program discription
	AS.pl
	
=head1 Command-line Option
	--psl		the psl File
	--junc		the junction File
	--id		the id File(geneID	trID)
	--dc		the depth and cov File
	--out 		the output dir
	--key		the keyname for the output file
	--num		the AS events number, 4 or 7
	--d			the depth for the posCoverage
				(default=2)
	--S			the max intron size
				(default=100000)
	--C			the coverage percentage for Retained Intron
				(default=90)
	--ip		the minimal seq-depth for Retained Intron compared to adjacent exons
				(default=15)
	--I			the idFile is true or false(T represents True, F represents False)
				(default=T)
	--help			
	
=head1 Usage
	perl AS.pl -help

=head1 Author
	Guo guangwu Email: guoguangwu@genomics.org.cn
	He zengquan Email: hezhengquan@genomics.org.cn
 	Zhang Jinbo Email: zhangjinbo@genomics.org.cn

=cut

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

my ($pslFile, $juncFile, $idFile, $depth_covFile);
my ($intron_size_max, $intron_cov_min, $intron_depth_min, $depth, $extend);
my ($keyName, $outDir, $I, $as_num);
my ($help);

GetOptions(
	"psl=s"=>\$pslFile,
	"junc=s"=>\$juncFile,
	"id=s"=>\$idFile,
	"I=s"=>\$I,
	"key=s"=>\$keyName,
	"dc=s"=>\$depth_covFile,
	"out=s"=>\$outDir,
	"num=i"=>\$as_num,
	"S=i"=>\$intron_size_max,
	"C=i"=>\$intron_cov_min,
	"ip=i"=>\$intron_depth_min,
	"e=i"=>\$extend,
	"d=i"=>\$depth,
	"help"=>\$help,
);

$as_num ||= 4;
$depth = 2 if (!defined $depth);
$extend = 5 if (!defined $extend);
$intron_size_max = 100000 if (!defined $intron_size_max);
$intron_cov_min = 90 if (!defined $intron_cov_min);
$intron_cov_min = $intron_cov_min / 100;
$intron_depth_min = 15 if (!defined $intron_depth_min);
$intron_depth_min = $intron_depth_min / 100;
$I ||= "T";
$outDir ||= "./";
$outDir = abs_path($outDir);
&make_Dir($outDir);

die `pod2text $0` if($help);
die `pod2text $0` if(!defined $pslFile || !defined $depth_covFile || !defined $keyName);
die `pod2text $0` if(!defined $juncFile || !defined $idFile);

######################### end input ################################
print STDERR "ip:$intron_depth_min\n";
print STDERR "d:$depth\n";
print STDERR "C:$intron_cov_min\n";
print STDERR "S:$intron_size_max\n";
print STDERR "e:$extend\n";

####################################################################
my (%junctions, %tr_gene, %junctions_a3a5);
my $sk_line = "SK";
my $mxe = "MXE";
my $a3ss = "A3SS";
my $a5ss = "A5SS";
my $AFE = "AFE";
my $ALE = "ALE";
my $RI = "RI";
my (%sk_uniq_exon, %mxe_uniq_exon);
my (%junc_right_biggest, %junc_left_smallest, %Exon_Edge_Left, %Exon_Edge_Right);
my (%a3_uniq_exon, %a5_uniq_exon, %a3a5_uniq_junc, %result_out_a3a5);
my (%afe_uniq, %ale_uniq, %result_afe, %result_ale);
my (%exons_depthes, %exons_covs);
my (%uniq_RI_exon, %uniq_RI_junction, %pos_depth_cov);
my (%overlap_genes, %overlap_introns);
my $pslFile_sort = $outDir."/".$keyName.".psl_sort";

############################### data pro ####################################
system("sort -k 14,14 -k 16n,16n $pslFile > $pslFile_sort");
&read_ID($idFile, \%tr_gene);
&read_JUNC($juncFile, \%junctions);
&AS();
system("rm $pslFile_sort");

############################## subroutine ###################################

sub AS{
	########### prapare ################
	&read_JUNC_A3A5($juncFile, \%junctions_a3a5);
	&read_PSL($pslFile_sort, \%junc_right_biggest, \%junc_left_smallest, \%result_out_a3a5);
	if($I eq "T"){
		&pro_overlap_genes();
	}
	&read_depth_cov( $depth_covFile, \%exons_depthes, \%exons_covs);
	###########  init   ################
	my (%result_out_mxe, %result_out_sk);
	my $sk_outFile = $outDir."/".$keyName.".SkippedExon.xls";
	my $a3_outFile = $outDir."/".$keyName.".A3SS.xls";
	my $a5_outFile = $outDir."/".$keyName.".A5SS.xls";
	my $RI_outFile = $outDir."/".$keyName.".RetainedIntron.xls";
	my $afe_outFile = $outDir."/".$keyName.".AlternativeFirstExon.xls";
	my $ale_outFile = $outDir."/".$keyName.".AlternativeLastExon.xls";
	my $mxe_outFile = $outDir."/".$keyName.".MutuallyExclusiveExon.xls";
	
	open  A5OUT, ">".$a5_outFile or die "Can't open A5_outFile : $a5_outFile";
	print A5OUT "3\)\tevents of alternative 5' splice site\n";
	print A5OUT "\t\t\t\t\t\tevent\n";
	print A5OUT "gene\tchromosome\tstrand\tconstitutive exon\t5' exon\t alternative 5' exon\n";
	
	open  A3OUT, ">".$a3_outFile or die "Can't open A3_outFile : $a3_outFile";
	print A3OUT "4\)\tevents of alternative 3' splice site\n";
	print A3OUT "\t\t\t\t\t\tevent\n";
	print A3OUT "gene\tchromosome\tstrand\tconstitutive exon\t3' exon\t alternative 3' exon\n";
	
	open  SK, ">".$sk_outFile or die "Can't open sk_outFile : $sk_outFile";
	print SK "1\) events of exon skipping\n";
	print SK "\t\t\t\t\t\tevent\n";
	print SK "gene\tchromsome\tstrand\tconstitutive exon\tinclusive exon\tconstitutive exon\n";
	
	open  RIOUT, ">".$RI_outFile or die "Can't open RI_outFile : $RI_outFile";
	print RIOUT "2\)\tevents of intron retaintion\n";
	print RIOUT "\t\t\t\t\t\tevent\n";
	print RIOUT "gene\tchromosome\tstrand\tretained intron\n";
	
	if($as_num == 7){
		open  MXEOUT, ">".$mxe_outFile or die "Can't open mxeFile : $mxe_outFile";
		print MXEOUT "5\)\tevents of mutually exclusive exon\n";
		print MXEOUT "\t\t\t\t\t\tevnet\n";
		print MXEOUT "gene\tchromsome\tstrand\tconstitutive exon\tinclusive exon\tconstitutive exon\n";
	
		open  AFE, ">".$afe_outFile or die "Can't open file : $afe_outFile";
		print AFE  "6\)\tevents of alternative first exon\n";
		print AFE  "\t\t\t\t\t\tevent\n";
		print AFE  "gene\tchromsome\tstrand\tconstitutive exon\tfirst exon\talternative first exon\n";

		open  ALE, ">".$ale_outFile or die "Can't open file : $ale_outFile";
		print ALE  "7\)\tevents of alternative last exon\n";
		print ALE  "\t\t\t\t\t\tevent\n";
		print ALE  "gene\tchromsome\tstrand\tconstitutive exon\tlast exon\talternative first exon\n";
	}
	
	
	open PSL, $pslFile or die "Can't open pslFile : $pslFile";
	print "Finding AS ...\n";
	while(<PSL>){
		chomp;
		my @temp = split('\t', $_);
		next if(@temp < 21);
		my $chr        = $temp[13];
		my $gene       = $temp[9];
		my $strand     = $temp[8];
		my $gene_start = $temp[15];
		my $gene_end   = $temp[16];
		my @sizes      = split( /,/, $temp[18] );
		my @starts     = split( /,/, $temp[20] );
		next if(@sizes != @starts);
		next if(!exists $tr_gene{$gene});
		my @exonlist   = ();
		my %result     = ( $sk_line, "", $mxe, "", $a3ss, "", $a5ss, "", $ALE, "", $AFE, "", $RI, "");
		for ( my $i = 0 ; $i < @starts ; $i++ ) {
			push @exonlist, [ $starts[$i] + 1, $starts[$i] + $sizes[$i] ];
		}
		
		## sk 
		&Skipped_Exon($chr, $strand, \%result, @exonlist);
		if ( $result{$sk_line} ne "" ) {
			my @arr = split( '\+', $result{$sk_line} );
			foreach my $sk (@arr) {
				my ($exons, $junction) = split(";", $sk);
				if(!exists $result_out_sk{$chr."-".$junction}){
					$result_out_sk{$chr."-".$junction} = $tr_gene{$gene}."\t".$chr."\t".$strand."\t".$exons;
				}
				else
				{
					my $tmp = $result_out_sk{$chr."-".$junction};
					my @arr_new = split('\t', $exons);
					my @arr_old = split('\t', $tmp);
					foreach(my $i = 0; $i < @arr_new; $i++){
						$arr_old[$i+3] .= ";".$arr_new[$i];
					}
					$result_out_sk{$chr."-".$junction} = join("\t", @arr_old);
				}
			}
		}
		## mxe
		&MXE($chr, $strand, \%result, @exonlist) if($as_num == 7);
		if ( $result{$mxe} ne "" ) {
			my @arr = split( '\+', $result{$mxe} );
			foreach my $mxe_line (@arr) {
				my ($exons, $junctions) = split(";", $mxe_line);
				if(!exists $result_out_mxe{$chr."-".$junctions}){
					$result_out_mxe{$chr."-".$junctions} = $tr_gene{$gene}."\t".$chr."\t".$strand."\t".$exons;
				}else{
					my $tmp = $result_out_mxe{$chr."-".$junctions};
					my @arr_new = split('\t', $exons);
					my @arr_old = split('\t', $tmp);
					foreach(my $i = 0; $i < @arr_new; $i++){
						$arr_old[$i+3] .= ";".$arr_new[$i];
					}
					$result_out_mxe{$chr."-".$junctions} = join("\t", @arr_old);
				}
			}
		}
		# a3 a5 
		&A3orA5SS($chr, \%result, $strand, @exonlist);
		my ($a5ss_tmp, $a3ss_tmp) = ($strand eq '+') ? ($result{$a5ss}, $result{$a3ss}) : ($result{$a3ss}, $result{$a5ss});
		if($a5ss_tmp ne ""){
			my @arr = split('\+', $a5ss_tmp);
			foreach my $a5 (@arr){
				print A5OUT $tr_gene{$gene}."\t".$chr."\t".$strand."\t";
				if($strand eq "+"){
					my @strand_out = split('\t', $a5);
					$a5 = $strand_out[1]."\t".$strand_out[0]."\t".$strand_out[2];
				}
				print A5OUT $a5."\n";
			}
		} 
		if($a3ss_tmp ne ""){
			my @arr = split('\+', $a3ss_tmp);
			foreach my $a3 (@arr){
				print A3OUT $tr_gene{$gene}."\t".$chr."\t".$strand."\t";
				if($strand eq "-"){
					my @strand_out = split('\t', $a3);
					$a3 = $strand_out[1]."\t".$strand_out[0]."\t".$strand_out[2];
				}
				print A3OUT $a3."\n";
			}
		}
		# RI
		&RetainedIntron( $chr, $strand, \%result, @exonlist );
		if ( $result{$RI} ne "" ) {
			my @out = split( ";", $result{$RI} );
			foreach my $Retain (@out) {
				print RIOUT $tr_gene{$gene}."\t".$chr."\t".$strand."\t";
				print RIOUT $Retain."\n";
			}
		}
		# afe ale
		&AFEorALE($chr, $strand, \%result, @exonlist) if($as_num == 7);
		my ($AFE_tmp, $ALE_tmp) = ($strand eq "+") ? ($result{$AFE}, $result{$ALE}) : ($result{$ALE}, $result{$AFE});
		if($AFE_tmp ne ""){
			my ($afe_exons, $afe_junctions) = split(";", $AFE_tmp);
			if(!exists $result_afe{$chr."-".$afe_junctions}){
				$result_afe{$chr."-".$afe_junctions} = $tr_gene{$gene}."\t".$chr."\t".$strand."\t".$afe_exons;
			}else{
				my $tmp_line = $result_afe{$chr."-".$afe_junctions};
				my @arr_new = split('\t', $afe_exons);
				my @arr_old = split('\t', $tmp_line);
				for(my $i = 0; $i < @arr_new; $i++){
					$arr_old[$i+3] .= ";".$arr_new[$i];
				}
				$result_afe{$chr."-".$afe_junctions} = join("\t", @arr_old);
			} 
		}
		if($ALE_tmp ne ""){
			my ($ale_exons, $ale_junctions) = split(";", $ALE_tmp);
			if(!exists $result_ale{$chr."-".$ale_junctions}){
				$result_ale{$chr."-".$ale_junctions} = $tr_gene{$gene}."\t".$chr."\t".$strand."\t".$ale_exons;
			}else{
				my $tmp_line = $result_ale{$chr."-".$ale_junctions};
				my @arr_new = split('\t', $ale_exons);
				my @arr_old = split('\t', $tmp_line);
				for(my $i = 0; $i < @arr_new; $i++){
					$arr_old[$i+3] .= ";".$arr_new[$i];
				}
				$result_ale{$chr."-".$ale_junctions} = join("\t", @arr_old);
			}
		}
	} # end while
	close(PSL);
	
	foreach my $out (keys %result_out_sk){
		print SK $result_out_sk{$out}."\n";
	}
	foreach my $out (keys %result_out_mxe){
		print MXEOUT $result_out_mxe{$out}."\n" if($as_num == 7);
	}
	foreach my $out(keys %result_afe){
		print AFE $result_afe{$out}."\n" if($as_num == 7);
	}
	foreach my $out(keys %result_ale){
		print ALE $result_ale{$out}."\n" if($as_num == 7);
	}
	close(SK);
	close(A3OUT);
	close(A5OUT);
	close(RIOUT);
	close(MXEOUT) if($as_num == 7);
	close(AFE) if($as_num == 7);
	close(ALE) if($as_num == 7);
	close(PSL);
	print "done with finding AS\n";
}



sub AFEorALE{
	my ($chr, $strand, $result_sub, @exons) = @_;
	my %before_hash = ();
	my %after_hash = ();
	my ($FrontS, $FrontE, $MiddleS, $MiddleE, $AfterS, $AfterE);
	for(my $b = 0; $b <= @exons-2; $b++){
		for(my $a = $b+1; $a <= @exons-1; $a++){
			($FrontS, $FrontE) = ($exons[$b]->[0], $exons[$b]->[1]);
			($AfterS, $AfterE) = ($exons[$a]->[0], $exons[$a]->[1]);
			if(exists $junctions{$chr}{"$FrontE-$AfterS"}){
				# before -> ....
				$before_hash{"$FrontS-$FrontE"}{"$AfterS-$AfterE"} = 1;
				# ...... <- after
				$after_hash{"$AfterS-$AfterE"}{"$FrontS-$FrontE"} = 1;
			}
		}
	}
	my $len = @exons;
	my $flag_afe = 0;
	my $flag_ale = 0;
	for(my $i = 0; $i <= @exons-1; $i++){
		if($flag_afe != 0 && $flag_ale != 0){
			last;
		}
		($FrontS, $FrontE) = ($exons[$i]->[0], $exons[$i]->[1]);
		($AfterS, $AfterE) = ($exons[$len-1-$i]->[0], $exons[$len-1-$i]->[1]);
		my $before = "$FrontS-$FrontE";
		my $after = "$AfterS-$AfterE";
		## AFE 1,2... -> n
		my $junctions_afe_edge = "";
		if(exists $after_hash{$after} && $flag_afe == 0){
			my $tmp = $after_hash{$after};
			my @arr = (); 
			if($strand eq "+"){
				@arr = sort keys %$tmp;
			}elsif($strand eq "-"){
				@arr = reverse sort keys %$tmp;
			}
			my $flag = 0;
			my $afe_line = "";
			my @afe_arr = ();
			if(@arr >= 2){
				$junctions_afe_edge = (split("-", $after))[0];
				foreach my $list(@arr){
					if(!exists $after_hash{$list}){
						$afe_line .= $list."\t";
						$junctions_afe_edge .= (split("-", $list))[1];
						$list =~ s/\-/\t/;
						push(@afe_arr, $list);
						$flag++;
					}
				}
			}
			if($flag >= 2){
				##
				if(exists $afe_uniq{"$chr\t$after\t$afe_line"}
				|| exists $ale_uniq{"$chr\t$after\t$afe_line"}){
					next;
				}
				if($strand eq "+"){
					$afe_uniq{"$chr\t$after\t$afe_line"} = 1;
				}elsif($strand eq "-"){
					$ale_uniq{"$chr\t$after\t$afe_line"} = 1;
				}
				$flag_afe = 1;  #
				##
				##
				$afe_line =~ s/\t$//;
				$result_sub->{$AFE} .= $after."\t".$afe_line.";".$junctions_afe_edge;
			}
		}
		
		## ALE 1 <- ... n-1, n
		my $junctions_ale_edge = "";
		if(exists $before_hash{$before} && $flag_ale == 0){
			my $tmp = $before_hash{$before};
			my @arr = ();
			if($strand eq "-"){
				@arr = sort keys %$tmp;
			}elsif($strand eq "+"){
				@arr = sort keys %$tmp;
			}
			my $flag = 0;
			my @ale_arr = ();
			my $ale_line = "";
			if(@arr >= 2){
				$junctions_ale_edge = (split("-", $before))[1];
				foreach my $list (@arr){
					if(!exists $before_hash{$list}){
						$ale_line .= $list."\t";
						$junctions_ale_edge .= (split("-", $list))[0];
						$list =~ s/\-/\t/;
						push(@ale_arr, $list);
						$flag++
					}
				}
				if($flag >= 2){
					##
					if(exists $afe_uniq{"$chr\t$before\t$ale_line"}
					|| exists $ale_uniq{"$chr\t$before\t$ale_line"}){
						next;
					}
					if($strand eq "+"){
						$ale_uniq{"$chr\t$before\t$ale_line"} = 1;
					}elsif($strand eq "-"){
						$afe_uniq{"$chr\t$before\t$ale_line"} = 1;
					}
					
					$flag_ale = 1; #
					##
					##
					$ale_line =~ s/\t$//;
					$result_sub->{$ALE} .= $before."\t".$ale_line.";".$junctions_ale_edge;
				}
			}	
		}	
	}
}

sub MXE{
	my ( $chr, $strand, $result_sub, @exons ) = @_;
	my ( %mxe_hash);
	my ( $FrontS, $FrontE, $MiddleS, $MiddleE, $AfterS, $AfterE );
	for ( my $F = 0 ; $F <= @exons - 3 ; $F++ ) {
		( $FrontS, $FrontE ) = ( $exons[$F]->[0], $exons[$F]->[1] );
		for ( my $M = $F + 1 ; $M <= @exons - 2 ; $M++ ) {
			( $MiddleS, $MiddleE ) = ( $exons[$M]->[0], $exons[$M]->[1] );
			next if(!exists $junctions{$chr}{"$FrontE-$MiddleS"});
			for ( my $A = $M + 2 ; $A <= @exons - 1 ; $A++ ) {
				( $AfterS, $AfterE ) = ( $exons[$A]->[0], $exons[$A]->[1] );
				
				next if(!exists $junctions{$chr}{"$MiddleE-$AfterS"});	
				$mxe_hash{"$FrontS-$FrontE"}{"$AfterS-$AfterE"}{"$MiddleS-$MiddleE"} = 1;
				
			}
		}
	}
	## MXE
	foreach my $first ( sort ( keys %mxe_hash ) ) {
		foreach my $second ( sort ( keys %{ $mxe_hash{$first} } ) ) {
			my $tmp = $mxe_hash{$first}{$second};
			my @arr = sort( keys %$tmp );
			next if ( @arr != 2 );
			my $mxe_exons = "";
			my $mxe_junctions = "";
			
			$mxe_exons = $first."\t";
			$mxe_junctions = (split("-", $first))[-1];
			
			my $tmp_1 = (split("-", $arr[0]))[-1];
			my $tmp_2 = (split("-", $arr[1]))[0];
			next if(exists $junctions{$chr}{"$tmp_1-$tmp_2"});
			
			foreach my $list (@arr){
				if($list =~ "-"){
					$mxe_junctions .= "-".$list;
					$mxe_exons .= $list."\t"; 
				}
			}
			$mxe_junctions .= "-".(split("-", $second))[0];
			$mxe_exons .= $second;
			
			if(!exists $mxe_uniq_exon{$mxe_exons}){
				$result_sub->{$mxe} .= $mxe_exons.";".$mxe_junctions."+";
				$mxe_uniq_exon{$mxe_exons} = 1;
			}
		}
	}
}

sub Skipped_Exon{
	my ( $chr, $strand, $result_sub, @exons ) = @_;
	my ( $FrontS, $FrontE, $MiddleS, $MiddleE, $AfterS, $AfterE );
	# for single
	for ( my $F = 0 ; $F <= @exons - 3 ; $F++ ) {
		
		( $FrontS, $FrontE ) = ( $exons[$F]->[0], $exons[$F]->[1] );
		
		for ( my $M = $F + 1 ; $M <= @exons - 2 ; $M++ ) {
			
			( $MiddleS, $MiddleE ) = ( $exons[$M]->[0], $exons[$M]->[1] );
			next if(!exists $junctions{$chr}{"$FrontE-$MiddleS"});
			
			for ( my $A = $M + 1 ; $A <= @exons - 1 ; $A++ ) {
				
				( $AfterS, $AfterE ) = ( $exons[$A]->[0], $exons[$A]->[1] );
				next if(!exists $junctions{$chr}{"$MiddleE-$AfterS"});
				next if(!exists $junctions{$chr}{"$FrontE-$AfterS"});		
				
				if(!exists $sk_uniq_exon{"$chr-$FrontS-$FrontE-$MiddleS-$MiddleE-$AfterS-$AfterE"}){
					
					$sk_uniq_exon{"$chr-$FrontS-$FrontE-$MiddleS-$MiddleE-$AfterS-$AfterE"} = 1;	
					my $sk_exons = "$FrontS-$FrontE\t$MiddleS-$MiddleE\t$AfterS-$AfterE";
					my $sk_junctions = "$FrontE-$MiddleS-$MiddleE-$AfterS";
					$result_sub->{$sk_line} .= $sk_exons.";".$sk_junctions."+";
					
				}
			}
		}
	}
	
	# for multiple
	my $flag_middle = 0;
	my $sk_exons = "";
	my $sk_junctions = "";
	
	my( $MiddleS_1, $MiddleE_1, $MiddleS_2, $MiddleE_2 );
	for ( my $F = 0; $F <= @exons - 4; $F++){
		
		( $FrontS, $FrontE ) = ( $exons[$F]->[0], $exons[$F]->[1] );
		
		for( my $A = $F + 3; $A <= @exons - 1; $A++){
			
			( $AfterS, $AfterE ) = ( $exons[$A]->[0], $exons[$A]->[1] );
			next if( !exists $junctions{$chr}{"$FrontE-$AfterS"});
			
			for( my $M = $F+1; $M <= $A - 2; $M++){
				
				( $MiddleS_1, $MiddleE_1 ) = ( $exons[$M]->[0], $exons[$M]->[1] );
			    ( $MiddleS_2, $MiddleE_2 ) = ( $exons[$M+1]->[0], $exons[$M+1]->[1]);
			    
				if( exists $junctions{$chr}{"$MiddleE_1-$MiddleS_2"}){
					if($flag_middle == 0){
						( $MiddleS, $MiddleE ) = ( $MiddleS_1, $MiddleE_2 );
						$flag_middle = 1;
						$sk_exons = "$MiddleS_1-$MiddleE_1,$MiddleS_2-$MiddleE_2";
					}else{
						( $MiddleS, $MiddleE ) = ( $MiddleS , $MiddleE_2 );
						$sk_exons = $sk_exons.","."$MiddleS_2-$MiddleE_2";
					}
				}else{
					if($flag_middle == 1){
						if(exists $junctions{$chr}{"$FrontE-$MiddleS"} 
						&& exists $junctions{$chr}{"$MiddleE-$AfterS"}){
							$sk_exons = "$FrontS-$FrontE\t".$sk_exons."\t$AfterS-$AfterE";
							$sk_junctions = "$FrontE-$MiddleS-$MiddleE-$AfterS";
							$result_sub->{$sk_line} .= $sk_exons.";".$sk_junctions."+";
						}
						$flag_middle = 0;
					}
				}
			}
			if($flag_middle == 1){
				if(exists $junctions{$chr}{"$FrontE-$MiddleS"} 
				&& exists $junctions{$chr}{"$MiddleE-$AfterS"}){
					$sk_exons = "$FrontS-$FrontE\t".$sk_exons."\t$AfterS-$AfterE";
					$sk_junctions = "$FrontE-$MiddleS-$MiddleE-$AfterS";
					$result_sub->{$sk_line} .= $sk_exons.";".$sk_junctions."+";
				}
				$flag_middle = 0;
			}
		}
	}
}

sub A3orA5SS{
	my ($chr, $result_sub, $strand, @exons) = @_;
	my ($FrontS, $FrontE, $AfterS, $AfterE);
	for (my $j = 0; $j <= @exons - 2; $j++){
		($FrontS, $FrontE) = ($exons[$j]->[0], $exons[$j]->[1]);
		($AfterS, $AfterE) = ($exons[$j+1]->[0], $exons[$j+1]->[1]);
		if(exists $junctions_a3a5{$chr}{$FrontE}{$AfterS})
		{
			# -1+ ->2  A5SS
			my $t5 = $junctions_a3a5{$chr}{$AfterS};
			my @a5 = keys %$t5;
			if(@a5 > 1){
				foreach my $per(@a5){
					if( $FrontS < $per && $per != $FrontE && $per < $AfterS 
					&&   exists $junc_left_smallest{$chr}{$FrontE}
					&&   $per < $junc_left_smallest{$chr}{$FrontE}
					&&  !exists $a5_uniq_exon{"$chr-$FrontS-$FrontE,$AfterS-$AfterE,$FrontS-$per"}
					&&  !exists $a3_uniq_exon{"$chr-$FrontS-$FrontE,$AfterS-$AfterE,$FrontS-$per"}
					&&  !exists $a3a5_uniq_junc{$chr.$FrontE.$per.$AfterS}
					)
					{
						my $uniq = join("", sort{$a <=> $b}($FrontE, $per, $AfterS));
						if(exists $a3a5_uniq_junc{$uniq}){
							next;
						}else{
							$a3a5_uniq_junc{$uniq} = 1;
						}
						#
						if(exists $Exon_Edge_Left{$chr}{$per}){
							my @lefts = split(";", $Exon_Edge_Left{$chr}{$per});
							my $flag_same = 0;
							for(my $pp = 0; $pp < @lefts; $pp++){
								if($lefts[$pp] < $FrontE ){
									$flag_same = 1;
									last;
								}
							}
							next if($flag_same == 0);
						}

						$a5_uniq_exon{"$chr-$FrontS-$FrontE,$AfterS-$AfterE,$FrontS-$per"} = 1;
					    $a3_uniq_exon{"$chr-$FrontS-$FrontE,$AfterS-$AfterE,$FrontS-$per"} = 1;
					    
						my $first = (split('\t', $result_out_a3a5{$chr}{"$FrontE-$AfterS"}))[0];
						my $per_line = "";
						my @tmp_per = split(";", $first);{
							foreach my $ele_FrontS(@tmp_per){
								$ele_FrontS = (split("-", $ele_FrontS))[0];
								$per_line .= "$ele_FrontS-$per;";
							}
						}
						$per_line =~ s/\;$//;
						$result_sub->{$a5ss} .= $result_out_a3a5{$chr}{"$FrontE-$AfterS"}."\t".$per_line."+";
						
					}
				}
			}
			# 2 <- +1-  A3SS
			my $t3 = $junctions_a3a5{$chr}{$FrontE};
			my @a3 = keys %$t3;
			if(@a3 > 1){
				foreach my $per(@a3){
					if( $FrontE < $per && $per != $AfterS && $per < $AfterE 
					&&   exists $junc_right_biggest{$chr}{$AfterS}
					&&   $per > $junc_right_biggest{$chr}{$AfterS}
					&&  !exists $a5_uniq_exon{"$chr-$FrontS-$FrontE,$AfterS-$AfterE,$FrontS-$per"}
					&&  !exists $a3_uniq_exon{"$chr-$FrontS-$FrontE,$AfterS-$AfterE,$FrontS-$per"}
					&&  !exists $a3a5_uniq_junc{$chr.$FrontE.$per.$AfterS}
					)
					{
						my $uniq = join("", sort{$a <=> $b}($FrontE, $per, $AfterS));
						if(exists $a3a5_uniq_junc{$uniq}){
							next;
						}else{
							$a3a5_uniq_junc{$uniq} = 1;
						}
						#
						if(exists $Exon_Edge_Right{$chr}{$per}){
							my @rights = split(";", $Exon_Edge_Right{$chr}{$per});
							my $flag_same = 0;
							for(my $pp = 0;$pp < @rights; $pp++){
								if($rights[$pp] > $AfterS){
									$flag_same = 1;
									last;
								}
							}
							next if($flag_same == 0);
						}

						$a5_uniq_exon{"$chr-$FrontS-$FrontE,$AfterS-$AfterE,$FrontS-$per"} = 1;
					    $a3_uniq_exon{"$chr-$FrontS-$FrontE,$AfterS-$AfterE,$FrontS-$per"} = 1;
						
						my $second = (split('\t', $result_out_a3a5{$chr}{"$FrontE-$AfterS"}))[1];
						my @tmp_per = split(";", $second);
						my $per_line = "";
						foreach my $ele_AfterE(@tmp_per){
							$ele_AfterE = (split("-", $ele_AfterE))[1];
							$per_line .= "$per-$ele_AfterE;";
						}
						$per_line =~ s/\;$//;
						$result_sub->{$a3ss} .= $result_out_a3a5{$chr}{"$FrontE-$AfterS"}."\t".$per_line."+";
					}
				}
			}	
		}
	}
}

sub read_PSL{
	my ($file, $right_biggest_sub, $left_smallest_sub, $result_out_sub) = @_;
	my ($FrontS, $FrontE, $AfterS, $AfterE);
	open IN, $file or die "can't open pslFile : $file";
	print "Reading the pslFile for A3A5 and RI\n";
	my $pre_gene = "";
	my $pre_end = 0;
	my $pre_line = "";
	while(<IN>){
		chomp;
		my @temp = split('\t', $_);
		my $chr = $temp[13];
		my $gene_id = $temp[9];
		my $strand=$temp[8];
		my @sizes=split(/,/,$temp[18]);
		my @starts=split(/,/,$temp[20]);
		next if(!exists $tr_gene{$gene_id});
		
		if($pre_gene eq ""){
			$pre_gene = $gene_id;
			$pre_end = $starts[-1] + $sizes[-1];
			$pre_line = $_;
		}else{	
			if(($starts[0]+1) <= $pre_end){
				if($tr_gene{$gene_id} ne $tr_gene{$pre_gene}){
					$overlap_genes{$tr_gene{$gene_id}} = 1;
					$overlap_genes{$tr_gene{$pre_gene}} = 1;
				}
			}
			$pre_gene = $gene_id;
			$pre_end = $starts[-1] + $sizes[-1] if(($starts[-1] + $sizes[-1]) > $pre_end);
			$pre_line = $_;
		}
		
		my @exonlist = ();
		next if(@starts < 2);
		for (my $i = 0; $i < @starts; $i++){
			push @exonlist, [ $starts[$i]+1, $starts[$i]+$sizes[$i] ];
		}
		
		for(my $j = 0; $j <= @exonlist - 2; $j++){
			($FrontS, $FrontE) = ($exonlist[$j]->[0], $exonlist[$j]->[1]);
			($AfterS, $AfterE) = ($exonlist[$j+1]->[0], $exonlist[$j+1]->[1]);
			# exon_edges
			if(!exists $Exon_Edge_Left{$chr}{$FrontE}){
				$Exon_Edge_Left{$chr}{$FrontE} = $FrontS;
			}else{
				$Exon_Edge_Left{$chr}{$FrontE} .= ";".$FrontS;
			}
			if(!exists $Exon_Edge_Right{$chr}{$FrontS}){
				$Exon_Edge_Right{$chr}{$FrontS} = $FrontE;
			}else{
				$Exon_Edge_Right{$chr}{$FrontS} .= ";".$FrontE;
			}
			if( $j == @exonlist-2 ){
				if(!exists $Exon_Edge_Left{$AfterE}){
					$Exon_Edge_Left{$chr}{$AfterE} = $AfterS;
				}else{
					$Exon_Edge_Left{$chr}{$AfterE} .= ";".$AfterS;
				}
				if(!exists $Exon_Edge_Right{$AfterS}){
					$Exon_Edge_Right{$chr}{$AfterS} = $AfterE;
				}else{
					$Exon_Edge_Right{$chr}{$AfterS} .= ";".$AfterE;
				}
			}

			if(!exists $result_out_sub->{$chr}{"$FrontE-$AfterS"}){
				$result_out_sub->{$chr}{"$FrontE-$AfterS"} = "$FrontS-$FrontE\t$AfterS-$AfterE";
			}else{
				my($first, $second) = split('\t', $result_out_sub->{$chr}{"$FrontE-$AfterS"});
				my $tmp_first = "$FrontS-$FrontE";
				my $tmp_second = "$AfterS-$AfterE";
				if($first =~ $tmp_first && $second =~ $tmp_second){
					next;
				}
				$first .= ";".$tmp_first;
				$second .= ";".$tmp_second;
				$result_out_sub->{$chr}{"$FrontE-$AfterS"} = $first."\t".$second;
			}
		 	
			if(!exists $right_biggest_sub->{$chr}{$AfterS}){
				$right_biggest_sub->{$chr}{$AfterS} = $FrontE;
			}else{
				$right_biggest_sub->{$chr}{$AfterS} = $FrontE
				if($right_biggest_sub->{$chr}{$AfterS} < $FrontE);
			}
			
			if(!exists $left_smallest_sub->{$chr}{$FrontE}){
				$left_smallest_sub->{$chr}{$FrontE} = $AfterS;
			}else{
				$left_smallest_sub->{$chr}{$FrontE} = $AfterS
				if($left_smallest_sub->{$chr}{$FrontE} > $AfterS);
			}
		}
	}
	print "Finish reading the pslFile for A3A5 and RI\n";
	close(IN);
}

sub read_JUNC_A3A5{
	my ($file, $junction_sub) = @_;
	my ($pos, @arr);
	open IN, $file or die "can't open junctionFile : $file";
	<IN>;
	print "Reading the junction File for A3A5\n";
	while(<IN>){
		chomp;
		@arr = split('\t', $_);
		next if(($arr[2]-$arr[1]) >= $intron_size_max);
		$junction_sub->{$arr[0]}{$arr[1]}{$arr[2]} = 1;
		$junction_sub->{$arr[0]}{$arr[2]}{$arr[1]} = 1;
	}
	close(IN);
	print "Finish reading the junction File for A3A5\n";
}

sub RetainedIntron{
	my ( $chr, $strand, $result_sub, @exons ) = @_;
	my ( $FrontS, $FrontE, $AfterS, $AfterE );
	for ( my $i = 0 ; $i <= @exons - 2 ; $i++ ) {
		( $FrontS, $FrontE ) = ( $exons[$i]->[0], $exons[$i]->[1] );
		( $AfterS, $AfterE ) = ( $exons[ $i + 1 ]->[0], $exons[ $i + 1 ]->[1] );
		# Intron Size Filter
		if($AfterS - $FrontE < 50 || $AfterS - $FrontE > $intron_size_max){
			next;
		}
		if ( exists $junctions{$chr}{"$FrontE-$AfterS"}
		 && !exists $uniq_RI_exon{"$chr-$FrontS-$FrontE-$AfterS-$AfterE"} 
		 && !exists $uniq_RI_junction{"$chr-$FrontE-$AfterS"} 
		 && exists $exons_depthes{$chr}{"$FrontS-$FrontE"}
		 && exists $exons_depthes{$chr}{"$AfterS-$AfterE"} 
		)
		{
			if($I eq "T" && exists $overlap_introns{$chr}{"$FrontE-$AfterS"}){
				next;
			}
			my $FrontE_I = $FrontE + 1;
			my $AfterS_I = $AfterS - 1;
			if ( exists $exons_covs{$chr}{"$FrontE_I-$AfterS_I"}
			  && exists $exons_depthes{$chr}{"$FrontE_I-$AfterS_I"}) {
				my $per = $exons_covs{$chr}{"$FrontE_I-$AfterS_I"};
				if ( $per >= $intron_cov_min ) {
					my $depth_F = $exons_depthes{$chr}{"$FrontS-$FrontE"}; 
					my $depth_A = $exons_depthes{$chr}{"$AfterS-$AfterE"};
					my $depth_I = $exons_depthes{$chr}{"$FrontE_I-$AfterS_I"};
					## intron
					if($depth_I <= ($depth_F * $intron_depth_min) && $depth_I <= ($depth_A * $intron_depth_min)){
						next;
					}
					$result_sub->{$RI} .= "$FrontE_I-$AfterS_I;";
					## check for the next
					$uniq_RI_exon{"$chr-$FrontS-$FrontE-$AfterS-$AfterE"} = 1;
					$uniq_RI_junction{"$chr-$FrontE-$AfterS"} = 1;
				}
			}
		}
	}
}

sub read_depth_cov{
	my ($file, $exons_depth_sub, $exons_cov_sub) = @_;
	open IN, $file or die "can't open depth_cov_file : $file";
	print "Reading the depth_cov File for RI\n";
	while(<IN>){
		chomp;
		my @arr = split('\t', $_);
		## chr start-end depth cov
		$exons_depth_sub->{$arr[0]}{$arr[1]} = $arr[2];
		$exons_cov_sub->{$arr[0]}{$arr[1]} = $arr[3];
	}
	close(IN);
	print "Finish reading the depth_cov File for RI\n";
}

sub read_JUNC{
	my ($file, $junctions_sub) = @_;
	print "Reading the junction File\n";
	open JUNC, "<".$file or die "Can't open junctions File : $file";
	my $head = <JUNC>;
	# Format: #
	# chr	side1	side2	strand	reads_number	junctions_name
	while(<JUNC>){
		chomp;
		my @arr = split('\t', $_);
		my $chr = $arr[0];
		my $strand = $arr[3];
		next if(($arr[2]-$arr[1]) >= $intron_size_max);
		my $pos = $arr[1]."-".$arr[2];
		$junctions_sub->{$chr}{$pos} = 1;
	}
	close(JUNC);
	print "Finish reading the junction File\n";
}

sub make_Dir{
	my $dir = shift;
	if(! -e $dir){
		mkdir $dir;
	}
}

sub read_ID{
	print "Reading the idFile\n";
	my ($file, $tr_gene_sub) = @_;
	open IN, $file or die "can't open idFile : $idFile";
	while(<IN>){
		chomp;
		my @arr = split('\t', $_);
		$tr_gene_sub->{$arr[1]} = $arr[0];
	}
	close(IN);
	print "Finish reading the idFile\n";
}

sub pro_overlap_genes{
	my %overlap_exons = ();
	my %overlap_trs = ();
	my %check_trs = ();
	open IN, $pslFile or die $!;
	while(<IN>){
		chomp;
		my @temp = split('\t', $_);
		my $line = $_;
		my $chr = $temp[13];
		my $gene_id = $temp[9];
		next if(!exists $tr_gene{$gene_id});
		next if(!exists $overlap_genes{$tr_gene{$gene_id}});
		my $strand=$temp[8];
		my $start = $temp[15];
		my $end = $temp[16];
		$overlap_exons{$chr}{$start}{$end} = $gene_id;
		$overlap_trs{$gene_id} = $line;
	}
	close(IN);
	
	foreach my $tr ( keys %overlap_trs ){
		my @temp = split('\t', $overlap_trs{$tr});
		my $chr = $temp[13];
		my $gene_id = $temp[9];
		my @sizes      = split( /,/, $temp[18] );
		my @starts     = split( /,/, $temp[20] );
		next if(@sizes == 1);
		for(my $i = 0; $i < @sizes - 1; $i++){
			my $FrontE = $starts[$i] + $sizes[$i];
			my $AfterS = $starts[$i+1] + 1;
			my $flag = 0;
			foreach my $start ( sort { $a <=> $b } keys %{$overlap_exons{$chr}} ){
				foreach my $end ( keys %{$overlap_exons{$chr}{$start}}){
					my $tmp_gene_id = $overlap_exons{$chr}{$start}{$end};
					if($tr_gene{$tmp_gene_id} eq $tr_gene{$gene_id}){
						next;
					}else{
						if( $FrontE >= $start && $FrontE <= $end ){
							$flag = 1;
							last;
						}elsif( $AfterS >= $start && $AfterS <= $end ){
							$flag = 1;
							last;
						}elsif( $AfterS < $start ){
							$flag = -1;
							last;
						}elsif( $FrontE > $end ){
							$flag = 0;
							next;
						}
					}
					last if($flag != 0);
				}
				last if($flag != 0);
			}
			if($flag == 1){
				$check_trs{$gene_id} = 1;
				$overlap_introns{$chr}{"$FrontE-$AfterS"} = 1;
			}
		}
	}
	open CHE, ">".$keyName.".out_psl" or die $!;
	foreach my $trID ( keys %check_trs ){
		print CHE $trID."\n";
	}
	close(CHE);
	
	#open TEMP, ">"."out.psl" or die $!;
	#for my $gene ( keys %overlap_trs ){
	#	print TEMP $overlap_trs{$gene}."\n";
	#}
	#close(TEMP);
}

