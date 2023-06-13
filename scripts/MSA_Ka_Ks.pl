#!/usr/bin/env perl

# This script takes two FASTA files with peptide and nucleotide CDS sequences 
# and computes median Ka and Ks values. Can be called with GNU parallel.
# Note: requires perl, clustal-omega, pal2nal.pl and subopt-kaks, set paths below.
# Note2: get these packages from:
# http://www.clustal.org/omega
# http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz
# https://github.com/hyphaltip/subopt-kaks

# example call:
# ./MSA_Ka_Ks.pl 5825_Brdisv1ABR21034417m.faa 5825_Brdisv1ABR21034417m.fna

# 2023 Bruno Contreras-Moreira (1) and Pablo Vinuesa (2):
# 1: http://www.eead.csic.es/compbio (Laboratory of Computational Biology, EEAD-CSIC, Spain)
# 2: http://www.ccg.unam.mx/~vinuesa (Center for Genomic Sciences, UNAM, Mexico)

use strict;
use File::Basename;

my $PERLEXE      = 'perl';
my $CLUSTALOEXE  = '~/soft/clustal-omega-1.2.4/src/clustalo'; 
my $PAL2NALEXE   = '~/github/get_phylomarkers/pal2nal.pl'; 
my $CDSPREALNEXE = '~/github/subopt-kaks/src/yn00_cds_prealigned';

my $KEEPTMPFILES = 0; # set to 1 for debugging or for conserving MSA files




my ($n_of_seqs, $pepfile, $nuclfile, @trash) = (0);

if(!defined($ARGV[1])) {
	die "# usage: $0 <peptide CDS FASTA filename> <nucleotide CDS FASTA filename>\n";
} else {

	($pepfile, $nuclfile) = @ARGV;

	if(!-e $pepfile){
		die "# ERROR: cannot find $pepfile\n";

	} if(!-e $nuclfile){
                die "# ERROR: cannot find $nuclfile\n";

        } else {
	
		# check input file contains enough sequences
		open(FAA, $pepfile);
		while(<FAA>) {
			if(/^>/) {
				$n_of_seqs++;
			}
		}
		close(FAA);

		if($n_of_seqs < 4) {
			print "#pepfile\tnuclfile\tnseqs\tmedian_Ka\tmedian_Ks\tmedian_omega\n";
			print "$pepfile\t$nuclfile\t$n_of_seqs\tNA\tNA\tNA\n";
			exit(0);
		}
   	}
}

my ($pepfile, $nuclfile) = @ARGV;
my (@trash);


## 1) get filename prefix
my $filename_prefix = basename($pepfile);


## 2) remove stop codons from peptide FASTA
my $nostop_filename = $filename_prefix . '.nostop.faa';
system("$PERLEXE -lne 's/\\*//; print' $pepfile > $nostop_filename");
if(!-e $nostop_filename) {
	die "# ERROR: cannot remove stop codons in $pepfile\n";
} else {
	push(@trash, $nostop_filename);
}


## 3) compute multiple alignment of proteins (MSA)
my $MSA_filename = $filename_prefix . '.aln.faa';
system("$CLUSTALOEXE -i $nostop_filename -o $MSA_filename");
if(!-e $MSA_filename) {
        die "# ERROR: cannot align $nostop_filename\n";
} else {
        push(@trash, $MSA_filename);
}


## 4) make codon alignment based on previous MSA
my $codon_filename = $filename_prefix . '.aln.fna';
system("$PAL2NALEXE -output fasta $MSA_filename $nuclfile > $codon_filename");
if(!-e $codon_filename) {
        die "# ERROR: cannot align $codon_filename\n";
} else {
        push(@trash, $codon_filename);
}


## 5) compute median Ka, Ks and omega
my (@data, @Ka, @Ks);
open(KAKS,"$CDSPREALNEXE $codon_filename |") ||
	die "# ERROR: cannot run $CDSPREALNEXE $codon_filename\n";
while(<KAKS>) {
	#SEQ1	SEQ2	dN	dS	OMEGA	N	S	kappa	t	LENGTH
	#Brdisv1ABR21034417m	Brdisv1ABR31027070m	0.000000	0.000000	99.000000	...
	next if(/^SEQ/);
	@data = split(/\t/,$_);

	next if($data[4] eq '-nan');

	push(@Ka, $data[2]);
	push(@Ks, $data[3]);	
}
close(KAKS);


## 6) print summary output
my $median_Ka = calc_median( \@Ka );
my $median_Ks = calc_median( \@Ks );
my $median_omega = 'NA';
if($median_Ks > 0) {
	$median_omega = sprintf("%1.4f", 
		$median_Ka / $median_Ks );
}

print "#pepfile\tnuclfile\tnseqs\tmedian_Ka\tmedian_Ks\tmedian_omega\n";
print "$pepfile\t$nuclfile\t$n_of_seqs\t$median_Ka\t$median_Ks\t$median_omega\n";


## 7) clean tmp files
unlink(@trash) if(!$KEEPTMPFILES);




# Takes ref to list of numbers and returns the median
sub calc_median {

  my ($dataref) = @_;

  my $mid = int(scalar(@$dataref)/2);
  my @sorted = sort {$a<=>$b} (@$dataref);

  if(scalar(@sorted) % 2) {
    return $sorted[ $mid ]
  }
  else {
    return sprintf("%1.4f",($sorted[$mid-1] + $sorted[$mid])/2)
  }
}
