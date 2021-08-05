#########################################################################
#      File Name: Gene_Pan_genome.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Mon 25 Mar 2019 05:00:00 PM CST
#########################################################################

#!/usr/bin/perl -w
use strict;


my $spe_list=shift;


my @spe;
&read_spe_list();

my %hash_result;
my %hash_del;
my @tmp_spe;
my $loc=0;
my %hash_loc;

&construct_pan();

open(OUT,">Pan_result.new");

#output header:
print OUT "LOC";
foreach my $spe (@tmp_spe){
	print OUT "\t$spe";
}
print OUT "\n";

#output pan:
foreach my $loc (sort {$a<=>$b} keys %hash_result){
	print OUT "$loc";
	foreach my $spe (@tmp_spe){
		my $out="-";
		if(exists $hash_result{$loc}{$spe}){
			$out=join "|",@{$hash_result{$loc}{$spe}};
		}
		print OUT "\t$out";
	}
	print OUT "\n";
}
close OUT;


sub read_gene_list{
	my $spe_tmp=shift;
	my @genes_tmp;
	open(IN,"$spe_tmp.bed.list");

	while(<IN>){
	    chomp;
	    push(@genes_tmp,$_);
	}
	close IN;
	return \@genes_tmp;
	
}


sub read_spe_list{
	open(IN,"$spe_list");
	while(<IN>){
		chomp;
		push(@spe,$_);
	}
	close IN;
}

sub read_blocks{
	my $spe1=shift;
	my $spe2=shift;
	my $block="$spe1.$spe2.i1.blocks";
	#print "$block\n";
	my %hash_syn_tmp;
	open(IN,"$block");
	while(<IN>){
		chomp;
		my @line = split /\s+/,$_;
		next if ($line[1] eq ".");
		$hash_syn_tmp{$line[0]}{$spe2}=$line[1];
	}
	close IN;
	return \%hash_syn_tmp;
}


sub construct_pan{
	print "Construct Pan genome...\n";
	my $spe = shift @spe;
	push(@tmp_spe,$spe);
	print "\t$spe\t";
	my $genes=&read_gene_list($spe);
	my @genes = @{$genes};
	
	foreach my $gene (@genes){
		$loc++;
		push(@{$hash_result{$loc}{$spe}},$gene);
		$hash_loc{$gene}=$loc;
	}
	print "last locus idfss : $loc\n";

	for(my $i=0;$i<@spe;$i++){

		$spe=$spe[$i];
		print "\t$spe\t";
		$genes=&read_gene_list($spe);
		@genes=@{$genes};
		
		foreach my $spe_pre (@tmp_spe){
			my $syn=&read_blocks($spe,$spe_pre);
			my %hash_syn=%{$syn};
			foreach my $gene (@genes){
				next if (exists $hash_loc{$gene});
				if(exists $hash_syn{$gene}{$spe_pre}){

					my $target=$hash_syn{$gene}{$spe_pre};
					my $loc_tmp=$hash_loc{$target};
					if(!defined $loc_tmp){
						print "Wrong!!!\t$gene\t$spe_pre\t$target no locID\n";
					}
					$hash_loc{$gene}=$loc_tmp;
					push(@{$hash_result{$loc_tmp}{$spe}},$gene);
				}
			}
		}
		foreach my $gene (@genes){
			next if (exists $hash_loc{$gene});
			#print "new loc : $gene\n";
			$loc++;
			$hash_loc{$gene}=$loc;
			push(@{$hash_result{$loc}{$spe}},$gene);
		}
		push(@tmp_spe,$spe);
		print "$#tmp_spe\tlast locus is : $loc\n";
	}
}
