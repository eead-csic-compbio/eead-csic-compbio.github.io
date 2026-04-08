#!/usr/bin/env perl

# Retrieves text taxonomies for one or more NCBI taxonIDs. Produces 2-column TSV: 
# 1) NCBI taxonID
# 2) full NCBI taxonomy in CSV string
# Requires https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz extracted in TAXONOMYPATH,
# particularly names.dmp and nodes.dmp (see below)
# 2026 Bruno Contreras-Moreira

use strict;
use warnings;

if(!$ARGV[0]){ 
  die "# usage: $0 <NCBI taxonID OR file with taxon IDs, one per file>\n" 
} 

## Edit as needed
my $TAXONOMYPATH = '/path/to/'; # path where taxdump.tar.gz was extracted
my $NAMESFILE = $TAXONOMYPATH . 'names.dmp';
my $NODESFILE = $TAXONOMYPATH . 'nodes.dmp';

my ($taxonID,$parentID,$taxonomy,%parent,%name,@taxonIDs,@lineage);

if($ARGV[0] =~ /^(\d+)$/) {
  push(@taxonIDs, $1)
} elsif(-e $ARGV[0]) {
  open(INFILE,"<",$ARGV[0]) || die "# $0 : cannot read $ARGV[0]\n";
  while(<INFILE>) {
    if(/^(\d+)$/){ 
      push(@taxonIDs, $1)
    }
  }
  close(INFILE);
} else {
  die "# wrong input, please use a valid NCBI taxonID or a file with several taxonIDs, one per line (see https://www.ncbi.nlm.nih.gov/taxonomy)\n"
}


## 1) retrieve parent graph
open(NODES,'<',$NODESFILE) || die "# $0 : cannot read $NODESFILE, please edit TAXONOMYPATH\n";
while(<NODES>) {
  if(/^(\d+)\t\|\t(\d+)\t/) {
    $parent{$1} = $2;
  }
}
close(NODES);

## 2) retrieve taxon names matching IDs from graph
open(NAMES,'<',$NAMESFILE) || die "# $0 : cannot read $NAMESFILE, please edit TAXONOMYPATH\n";
while(<NAMES>) {
  if(/^(\d+)\t\|\s+([^\t]+).*?scientific name/) {	  
    $name{$1} = $2;
  }
}
close(NAMES);

## 3) get full taxonomy (names) for input taxon IDs
foreach $taxonID (@taxonIDs) {

  if(!defined($parent{$taxonID})) {
    warn "# ERROR: cannot find $taxonID in $NODESFILE , skip it\n";
    print "$taxonID\tNA\n";
    next
  }  

  ($taxonomy,@lineage) = ('',());

  $parentID = $parent{$taxonID};
  push(@lineage,$taxonID);
  while($parentID ne $taxonID && $parentID > 1) {
    push(@lineage,$parentID);
    $parentID = $parent{$parentID};
  }

  foreach $parentID (reverse(@lineage)) {
    $taxonomy .= "$name{$parentID},"  
  }

  print "$taxonID\t$taxonomy\n";
}
