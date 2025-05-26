#!/usr/bin/perl
use strict;

# Utility to create HOR_14121.h.bed from PHG's vcf_dbs/hvcf_files/HOR_14121.h.vcf.gz,
# see align2graph.py
# 
# J Sarria, B Contreras-Moreira
# Copyright [2024-25] Estacion Experimental de Aula Dei-CSIC

my ($chr,$start,$end,$strand,$multi1toN,$genome,$sum); # indiv genome
my ($gchr,$gstart,$gend,$gstrand,$gsum,$intvl); # graph

# 1toN regions matches examples:
#
#zcat *vcf.gz |perl -lne 'if(/Regions=\"/){ print }' 
#1to2
##ALT=<ID=dee16bf3422778591b76cc0f92b3e389,Description="haplotype data for line: HOR_12184",Source="vcf_dbs//assemblies.agc",SampleName=HOR_12184,Regions="chr4H_OX460225.1:596550945-596550945,chr4H_OX460225.1:596550944-596536672",Checksum=dee16bf3422778591b76cc0f92b3e389,RefChecksum=2aae972083735768391317e5d25e1dd0,RefRange=chr4H_LR890099.1:604188691-604201641>
#1to3
##ALT=<ID=b7687097e22e772e0dadc66ac82f5427,Description="haplotype data for line: HOR_13942",Source="vcf_dbs//assemblies.agc",SampleName=HOR_13942,Regions="chr1H:6803400-6755916,chr1H:6755915-6756737,chr1H:6755092-6615635",Checksum=b7687097e22e772e0dadc66ac82f5427,RefChecksum=4c00d9dc7ec544cb84ee8341afb2e98a,RefRange=chr1H:5160593-5341120>

while(<>) {

	# single segment
	if(/SampleName=([^,]+),Regions=([^:]+):(\d+)-(\d+),Checksum=([^,]+),RefChecksum=([^,]+),RefRange=([^:]+):(\d+)-(\d+)/) { 

		($genome, $chr,$start,$end,$strand,$sum) = ($1, $2, $3, $4, "+", $5); 
		if($start > $end) {
			$strand="-";
			($start,$end)=($end,$start);
		}

		print "$chr\t$start\t$end\t$strand\t$sum\t$genome\t$7\t$8\t$9\t$6\n";

	} elsif(/SampleName=([^,]+),Regions="([^"]+)",Checksum=([^,]+),RefChecksum=([^,]+),RefRange=([^:]+):(\d+)-(\d+)/) {

		# 2+ segments, each with corresponding BED interval
		($genome,$multi1toN,$gsum,$sum,$gchr,$gstart,$gend) = ($1,$2,$3,$4,$5,$6,$7);
		
		foreach $intvl (split(/,/,$multi1toN)) {
			if($intvl =~ m/([^:]+):(\d+)-(\d+)/) {	
				($chr,$start,$end,$strand) = ($1, $2, $3, "+");
				if($start > $end) {
                 			$strand="-";
	                		($start,$end)=($end,$start);
        		        }
			}

			print "$chr\t$start\t$end\t$strand\t$sum\t$genome\t$gchr\t$gstart\t$gend\t$gsum\n";
		}
	}		
}
