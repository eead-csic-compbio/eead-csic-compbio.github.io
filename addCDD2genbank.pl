#!/usr/bin/env perl

# This script adds domain annotations from CDD to protein sequences
# contained in CDS features within input GenBank file 

# 2017 Bruno Contreras-Moreira (1) and Pablo Vinuesa (2):
# 1: http://www.eead.csic.es/compbio (Laboratory of Computational Biology, EEAD/CSIC/Fundacion ARAID, Spain)
# 2: http://www.ccg.unam.mx/~vinuesa (Center for Genomic Sciences, UNAM, Mexico)

# Dependencies: bioperl

use strict;
use Getopt::Std;
use File::Temp qw/ tempfile /;
use Bio::SeqIO;

my $CDDERROR = 'need valid -c path to cdd_delta; ftp://ftp.ncbi.nlm.nih.gov/blast/db/cdd_delta.tar.gz';
my $BLASTERROR = 'need valid -p path to deltablast; ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/';

my ($MAKEBLASTDB,$DELTABLAST) = ( 'makeblastdb', 'deltablast' );
my ($INP_threads,$INP_evalue,$INP_CDD,$INP_blastpath,$INP_infile,$INP_outfile) = (2,0.00001,'cdd_delta');
my ($start,$end,$strand,$pepseq,$domstart,$domend,$inference,$note,$cddid);
my (%opts,%annotations,@trash);

getopts('hc:o:i:p:t:e:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print "\n[options]: \n";
  print "-h this message\n";
  print "-i input GenBank file                        (required)\n";
  print "-o output GenBank file with CDD annotations  (required)\n";
  print "-p path to deltablast & makeblastdb          \n";
  print "-c path to CDD formatted for deltablast      \n";
  print "-e max E-value of CDD domain matches         \n";
  print "-t number of threads/CPU cores               \n";
  exit(-1);
}

## 1) check params
if(defined($opts{'i'})){ $INP_infile = $opts{'i'} }
else{ die "\n# $0 : need -i infile, exit\n"; }

if(defined($opts{'o'})){ $INP_outfile = $opts{'o'} }
else{ die "\n# $0 : need -o outfile, exit\n"; }

if(defined($opts{'c'}))
{ 
  $INP_CDD = $opts{'c'};
  if(!-s $INP_CDD.'.pin'){ die "# $0 : $CDDERROR\n" }
}
else{ die "# $0 : $CDDERROR\n" }

if(defined($opts{'t'}) && $opts{'t'} >= 1){ $INP_threads = int($opts{'t'}) }

if(defined($opts{'e'}) && $opts{'e'} > 0){ $INP_evalue = $opts{'e'} }

if(defined($opts{'p'}))
{
  $INP_blastpath = $opts{'p'};
  $MAKEBLASTDB = $INP_blastpath.'/makeblastdb';
  $DELTABLAST  = $INP_blastpath.'/deltablast';
}

warn "\n# $0 -i $INP_infile -c $INP_CDD -o $INP_outfile -t $INP_threads -e $INP_evalue-p $INP_blastpath\n";

## 2) read input GenBank and extract peptide sequences
my $extracted_peptides = 0;
my ($tmpfh, $tmpfilename) = tempfile( unlink=>1 );

my $ingbk = new Bio::SeqIO(-file => $INP_infile, -format => 'genbank' );
while( my $seq = $ingbk->next_seq()) 
{    
  foreach my $f ($seq->get_SeqFeatures) 
  {
    next if($f->primary_tag() ne 'CDS' || !$f->has_tag('translation'));

    ($start,$end,$strand,$pepseq) = 
      ( $f->start(), $f->end(), $f->location()->strand(), ($f->get_tag_values('translation'))[0] );

    # write to tmp file using coords as identifier
    print $tmpfh ">$strand:$start:$end\n$pepseq\n";
    $extracted_peptides++;
  }
}
warn "# $0 : extracted $extracted_peptides peptide sequences to $tmpfilename\n";

## 3) annotate CDD domains in extracted peptides
open(MAKEDB,"$MAKEBLASTDB -in $tmpfilename -dbtype prot |") || 
  die "# $0 : need a valid \$MAKEBLASTDB path (please use -p path)\n";
while(<MAKEDB>){} # format sequence db
close(MAKEDB);
push(@trash,"$tmpfilename.phr","$tmpfilename.pin","$tmpfilename.psq");

my $annotated_domains = 0;
open(DELTA,"$DELTABLAST -num_threads $INP_threads -query $tmpfilename -db $tmpfilename ".
  "-rpsdb $INP_CDD -evalue $INP_evalue -show_domain_hits -outfmt \"6 std qlen stitle\" |") ||
  die "# $0 : need a valid \$DELTABLAST path (please use -p path)\n";
while(my $line = <DELTA>)
{
  #qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen stitle
  #1:46199:46990 gnl|CDD|214794  28.57 112 63  5 10  118 1 98  1e-04 40.0  263 smart00731, SprT...
  #-1:32308:35544 gnl|CDD|279750  38.00 50  27  1 854 899 178 227 0.001 40.8  1078  pfam01443, Viral...
  chomp;
  my @data = split(/\t/,$line);
  next if($data[1] !~ /CDD/); # skip sequence hits, it's domains we are after
  
  # compute domain boundaries in the genome
  ($strand,$start,$end) = split(/:/,$data[0],3); 
  $cddid = (split(/\|/,$data[1],3))[-1];

  if($strand == 1)
  {
    $domstart = $start + (($data[6]-1)*3);
    $domend   = $end   - (($data[12]-$data[7])*3);
  }
  else
  {
    $domstart = $end   - (($data[6]-1)*3);
    $domend   = $start + (($data[12]-$data[7])*3)
  }

  # format annotation for GenBank and save
  $inference = "protein motif:CDD:$cddid";
  $note = $data[13];
  push( @{$annotations{$strand}{$start}{$end}}, [ $inference, $note ] );
  $annotated_domains++;
}
close(DELTA);

warn "# $0 : annotated $annotated_domains CDD domains\n";

## 4) add domains as misc_features linked to original CDS features
$ingbk = Bio::SeqIO->new( ); '-format' => 'embl', '-file' => "$ARGV[0]");

my $seqout = new Bio::SeqIO('-format' => 'genbank', '-file' => ">$ARGV[1]");
while( my $seq = $seqio->next_seq) {
    $seqout->write_seq($seq)
}


$ingbk = new Bio::SeqIO(-file => $INP_infile, -format => 'genbank' );
while( my $seq = $ingbk->next_seq())
{
  foreach my $f ($seq->get_SeqFeatures)
  {
    next if($f->primary_tag() ne 'CDS' || !$f->has_tag('translation'));

    ($start,$end,$strand) =
      ( $f->start(), $f->end(), $f->location()->strand() );

    print "$start,$end,$strand\n";

    my $misc_feat = Bio::SeqFeature::Generic->new(
      -start =>$start,
      -end   =>$end,
      -primary_tag => 'misc_feature',
      -tag => { 
        inference => 'inference',
        note      => 'nota'
    });

    $ingbk->add_SeqFeature($misc_feat);
  }
}




## 5) clean tmp files
unlink(@trash);
