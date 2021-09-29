#!/usr/bin/env perl 
use strict;
use Getopt::Std;

# Reads in XMFA file produced by progressive MAUVE and produces a multi-FASTA file 
# containing a multiple sequence alignment (MSA) in which each input genome takes 
# a single line with syntenic blocks separated by $BLOCKSEPARATOR chars. 
# Columns of the resulting MSA take the corresponding coordinates of the first input genome.
# Blocks where first/reference genomes is absent have coordinates set to zero.
# 
# See http://darlinglab.org/mauve/user-guide/files.html for XMFA format definition
#
# B Contreras-Moreira, U Alonso EEAD-CSIC, Zaragoza, Spain

my $BLOCKSEPARATOR = 'n';
my $GAPCHAR = '-';

my ($input_XMFA_file,$output_FASTA_file,$output_LOG_file) = ('','','');
my $INP_polymorphic = 0;
my $INP_noindels = 0;
my $INP_biallelic = 0;
my %opts;

getopts('hpnbi:o:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print "\nusage: $0 [options]\n\n";
  print " -h this message\n";
  print " -i input XMFA file from progressive MAUVE\n";
  print " -o output multi-FASTA file\n";
  print " -p only polymorphic loci kept in the alignment [optional]\n";
  print " -n no indels are allowed, requires -p          [optional, gaps are valid alleles by default]\n";
  print " -b only biallelic, requires -p                 [optional]\n\n";
  exit(0);
}

if(defined($opts{'i'})){
  $input_XMFA_file = $opts{'i'};
  die "# EXIT : need a valid input XMFA file\n" if(!-e $input_XMFA_file);
}
else{ die "# EXIT : need a valid input XMFA file\n" }

if(defined($opts{'o'})){
    $output_FASTA_file = $opts{'o'};
    $output_LOG_file = $output_FASTA_file.'.log';
}
else{ die "# EXIT : need a valid output FASTA file\n" }

if(defined($opts{'p'})){ 
  $INP_polymorphic = 1;

  if(defined($opts{'n'})){ $INP_noindels = 1 }
  if(defined($opts{'b'})){ $INP_biallelic = 1 }
}

print "# $0 -i $input_XMFA_file -o $output_FASTA_file -p $INP_polymorphic -n $INP_noindels -b $INP_biallelic\n";
print "# \$BLOCKSEPARATOR='$BLOCKSEPARATOR'\n\n";

###################################################

my ($genome_index,$genome_filename,$line,$first_block_index,$block_coord);
my ($block_length,$raw_MSA_length,$final_MSA_length) = (0,0,0);
my ($is_polymorphic,$is_biallelic,$is_gapped) = (0,0,0);
my ($block_start,$block_end,$current_coord) = (0,0,0);
my (%finalMSA,%MSA,%block,%fullname);
my (@indexes,@MSAmatrix,@coords);

open(XMFA,$input_XMFA_file) || die "# ERROR: cannot read $input_XMFA_file\n";
while($line = <XMFA>){
  if($line =~ m/^#/){
    if($line =~ m/#Sequence(\d+)File\s+(\S+)/){
      ($genome_index,$genome_filename) = ($1,$2);
      $fullname{$genome_index} = $genome_filename;

      # save indexes as they appear in input file
      push(@indexes, $genome_index);
    } else{ next }
  }
  elsif($line =~ m/>\s*(\d+):(\d+)-(\d+)/){
    $genome_index = $1;

    # store original coordinates of this block
    if($genome_index == 1){
      ($block_start,$block_end) = ($2,$3);
    }
  }
  elsif($line =~ m/^=/){
    # process the last parsed block of sequences which might not include all genomes

    # check block length and increase total
    $first_block_index = (keys(%block))[0];
    $block_length = length($block{$first_block_index});
    $raw_MSA_length += $block_length;

    # keep track of block coordinates
    my @block_ref_seq = split(//,$block{'1'});
    for $block_coord (0 .. $block_length-1){
      if($block_start>0){
        # update current_coord only if first/reference genome is aligned in block
        $current_coord = $block_start;
        if($block_ref_seq[$block_coord] ne $GAPCHAR){
          $current_coord = $block_start++;
        }  
      }
      else{ $current_coord = '0' }
      push(@coords,$current_coord);
    }

    foreach $genome_index (@indexes){
      # for genomes present in block 
      if(defined($block{$genome_index})){
        $MSA{$genome_index} .= $block{$genome_index}.$BLOCKSEPARATOR;
      }
      else{
        $block{$genome_index} = '-' x $block_length;
        $MSA{$genome_index} .= $block{$genome_index}.$BLOCKSEPARATOR;
      }
    }  

    %block = ();
    $block_length = 0;
    ($block_start,$block_end) = (0,0);
  }
  else{
    chomp($line);
    $block{$genome_index} .= $line;
  }
}
close(XMFA);

printf("# input genomes=%d\n",scalar(@indexes));
printf("# input loci=%d input coordinates=%d\n\n",$raw_MSA_length,scalar(@coords)); 

## filter columns/loci from the raw multiple alignment (MSA)
# split MSA rows into single bases
foreach $genome_index (@indexes){
  $MSAmatrix[$genome_index] = [ split(//,$MSA{$genome_index}) ];
}

# loop through and filter out desired columns
open(LOG,">",$output_LOG_file) || die "# ERROR: cannot create $output_LOG_file\n";
print LOG "#column\tgenomic_coordinate\n";

for my $col (0 .. $raw_MSA_length-1){
  ($is_polymorphic,$is_biallelic,$is_gapped) = (0,0,0);
  my %alleles;

  foreach $genome_index (@indexes){
    $alleles{ $MSAmatrix[$genome_index][$col] }++;
  }

  if($INP_polymorphic == 1){
    if(scalar(keys(%alleles)) == 1){ next } # non-polymorphic column
    if($INP_noindels == 1 && defined($alleles{$GAPCHAR})){ next } # indel column
    if($INP_biallelic == 1 && scalar(keys(%alleles)) > 2){ next } # non-biallelic column
  }

  # concat this column to final MSA if not skipped by previous tests
  foreach $genome_index (@indexes){
    $finalMSA{$genome_index} .= $MSAmatrix[$genome_index][$col];
  }

  # update final MSA length
  $final_MSA_length++;

  # save coordinate of this column in reference/first genome
  print LOG "$final_MSA_length\t$coords[$col]\n";
}

close(LOG);

printf("# output loci=%d\n",$final_MSA_length);

## print output
open(OUT,">",$output_FASTA_file) || die "# ERROR: cannot create $output_FASTA_file\n";

foreach $genome_index (@indexes){
    print OUT ">$fullname{$genome_index}\n$finalMSA{$genome_index}\n";
}

close(OUT);

print "# output file: $output_FASTA_file\n";
print "# log file: $output_LOG_file\n";
