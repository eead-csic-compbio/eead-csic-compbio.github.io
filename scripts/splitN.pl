#!/usr/bin/env perl
use strict;
use warnings;

# Reads i) FASTA file with DNA scaffolds and ii) integer with max chr size (bp) to output a multi-FASTA string
# where sequences longer than max size GB are split into subsequences as needed, adding chunk number to headers.
# Note: N runs used to split are lost

# 2026 Bruno Contreras-Moreira

if(!$ARGV[1]){ 
  die "# usage: $0 <FASTA file> <max chr size (bp)>\n";
}

my ($infasta, $maxchrsize) = @ARGV;
my ($seq, $header, $chunk, $tot, $len, $chunkn) = ('','','',0);

open(FA,'<',$infasta) || 
  die "# ERROR: cannot read $infasta\n";
while(<FA>) {

  if (/^>(\S+)/) {

    $len = length($seq);
    if($len > 0) {
      printf(">%s_%d %d\n%s",
        $header,$chunkn,$len,wrap($seq));
      $tot+=$len;
    }

    $header = $1; 
    ($chunkn,$seq) = (1,'');

  } else { 
    chomp; 
    $seq .= $_;
    $len = length($seq);

    if($len > $maxchrsize) {
      ($chunk,$seq) = split(/N+([^N]*)$/, $seq);
      $len = length($chunk);
      printf(">%s_%d %d\n%s",
        $header,$chunkn,$len,wrap($chunk));
      $tot+=$len;
      $chunkn++;
    }
  }  
}
close(FA);

$len = length($seq);
if($len > 0) {
  printf(">%s_%d %d\n%s",
    $header,$chunkn,$len,wrap($seq));
  $tot+=$len;
}

warn "\n# INFO: sequence length = $tot\n";


sub wrap { 
  my ($seq) = @_;
    
  my $wrapped = '';
  while ($seq =~ m/(\S{1,80})/g) {
    $wrapped .= "$1\n";
  }

  return $wrapped
}
