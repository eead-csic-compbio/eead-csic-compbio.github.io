#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;

# takes i) a list of sample names and ii) a VCF file
# returns a matrix of maerker SNPs in those samples
# by Anna and Bruno at EBI hackathon 2019

my $missing = '-/-';
my (%opts,$INP_sample_list,$INP_VCF_file);
my ($INP_only_polymorphic,$INP_min_freq,$INP_identity);
my (@samples,%sample2col,@wanted_cols);
my ($sample,$sample_vcf,$line,$col,$gt);
my $n_of_markers = 0;

getopts('hpm:s:v:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print   "\nusage: $0 [options]\n\n";
  print   "-h this message\n";
  print   "-v VCF file                             (required)\n";
  print   "-s sample list file, one per line       (optional, all taken by default)\n";
  print   "-p get only polymorphic markers         (optional)\n";
  print   "-m min genotype frequency               (optional, requires -p)\n";
  #print   "-i compute identity among samples       (optional)\n";
  exit(0);
}

if(defined($opts{'v'})){  $INP_VCF_file = $opts{'v'}; }
else{ die "# ERROR: neead an input VCF file\n" }

if(defined($opts{'p'})){  
	$INP_only_polymorphic = 1;
	
	if(defined($opts{'m'})){  $INP_min_freq = $opts{'m'} }
	else{ $INP_min_freq = 0 }
	
} else{ $INP_only_polymorphic = 0 }


## read in samples in list, optional
if(defined($opts{'s'})){  
	$INP_sample_list = $opts{'s'}; 

	open(SAMPLES,"<",$INP_sample_list) ||
		die "# ERROR: cannot read $INP_sample_list\n";
	while($line = <SAMPLES>) {
		next if($line =~ /^#/); #comment
		if($line =~ /^(\S+)/){ #regex for first non-blank fragment
			$sample = $1;
			push(@samples, $sample);
		}
	}
	close(SAMPLES);

	#printf("# total samples in list: %d\n",scalar(@samples));
} else {
	$INP_sample_list = '';
	#print "# take all samples\n";
}

## parse VCF file and extract samples in the list
# http://www.internationalgenome.org/wiki/Analysis/vcf4.0/
open(VCF,"<",$INP_VCF_file) || die "# ERROR: cannot read VCF file\n";
while($line = <VCF>) {

	my @vcfdata = split(/\t/,$line);

	if($vcfdata[0] eq '#CHROM') { #header, 1st line that matters in VCF file
		
		# parse samples and match their corresponding column index
		my $sample_list_string = $vcfdata[8];
		my @sample_list = split(/\s+/,$sample_list_string);

		for($col=1;$col<scalar(@sample_list);$col++) {
			$sample_vcf = $sample_list[$col];
			$sample2col{$sample_vcf} = $col; # index to the right column	
		}
		
		#printf("# total GT columns: %d\n",scalar(@sample_list)); 

		# find indexes of samples in input list
		if(scalar(@samples) > 0) {
			foreach $sample (@samples) {
				foreach $sample_vcf (@sample_list){
					if($sample_vcf =~ m/\Q$sample\E/) {
						push(@wanted_cols, $sample2col{$sample_vcf});
						#print "# $sample_vcf $sample2col{$sample_vcf}\n";	
					}
				}
			}

			#make sure at least one column was picked up
			if(scalar(@wanted_cols) == 0){
				die "# ERROR: cannot match your sample names in VCF file\n";
			} 

		} else { 
			for($col=1;$col<scalar(@sample_list);$col++) {
                        	$sample_vcf = $sample_list[$col];
				push(@wanted_cols, $sample2col{$sample_vcf});
				#print "# $sample_vcf $sample2col{$sample_vcf}\n";
			}
		} #print scalar(@wanted_cols); exit;

		# print metadata
		print "##fileformat=VCFv4.0\n";
		print "##reference=hopefullyIWGCCRefseq1.1\n";
		print "##source: $0 -s $INP_sample_list -v $INP_VCF_file ".
			"-p $INP_only_polymorphic -m $INP_min_freq\n";

		# print header
		foreach $col (0 .. 7) {
                        print "$vcfdata[$col]\t";
                }
		foreach $col (@wanted_cols){
                        print "$sample_list[$col]\t";
                }
		print "\n";

	} else { #actually read genotypes/SNPs

		my $subsetVCFline = '';

		# unchanged first 8 columns from VCF file
		foreach $col (0 .. 7) { 
			$subsetVCFline .= "$vcfdata[$col]\t";
		}
		
		# split genotype calls/markers and print only those wanted
		my $sample_list_string = $vcfdata[8];
                my @sample_list = split(/\s+/,$sample_list_string);
		my ($total_gt, $min_freq) = (0, 1.0);
		my %genotypes;
		foreach $col (@wanted_cols){
			$gt = $sample_list[$col];
			$subsetVCFline .= "$gt\t";
			if($gt ne $missing){
				$genotypes{$gt}++;
				$total_gt++;
			}
		}

		if($INP_min_freq) {
			my $freq;
			foreach $gt (keys(%genotypes)){
				$freq = sprintf("%1.3f",$genotypes{$gt} / $total_gt);
				#print "$gt $freq $genotypes{$gt} $total_gt\n";
				if($freq < $min_freq){ $min_freq = $freq }	
			}
		}

		# check the marker is polymorphic if required
		if($INP_only_polymorphic == 0) {
			print "$subsetVCFline\n";
			$n_of_markers++;
		} elsif(scalar(keys(%genotypes)) > 1) {
			if($min_freq >= $INP_min_freq){
				print "$subsetVCFline\n";
                        	$n_of_markers++;
			}
		}

	}	
}
close(VCF);

#print "# total markers = $n_of_markers\n";

