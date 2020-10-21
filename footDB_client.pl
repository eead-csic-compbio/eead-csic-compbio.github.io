#!/usr/bin/env perl

# This script allows you to query http://floresta.eead.csic.es/footprintdb
# from the command-line. It uses Perl module SOAP::Lite, which you'll find at 
# https://metacpan.org/pod/release/PHRED/SOAP-Lite-1.27/lib/SOAP/Lite.pm
# and install with $ sudo cpan -i SOAP::Lite

# Note: please space your queries if you plan to submit several jobs, 
# such as 1 every 10s for instance

# Read more at
# https://bioinfoperl.blogspot.pt/2017/10/soap-interface-of-footprintdb.html

use strict;
use warnings;
use Getopt::Std; 
use SOAP::Lite;

my %opts;
my ($username,$input,$server) = ('','');

getopts('hu:m:p:k:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
	print "\nusage: $0 [options]\n\n";
	print "-h this message\n";
	print "-m motif file in TRANSFAC format, example: -m matrix.tf\n";
	print "-p peptide sequence in FASTA format, example: -p sequence.faa\n";
	print "-k keyword for text query, example: -k myb\n";
	print "-u registered footprintDB username, optional\n\n";
	exit(0);
} 

if(!defined($opts{'m'}) && !defined($opts{'p'}) && !defined($opts{'k'})){ 
	print "# ERROR: choose one of -m, -p -or -k\n";
	exit(1);
}
else { # prepare Web Services connection
	$server = SOAP::Lite
		-> uri('footprintdb')
		-> proxy('http://floresta.eead.csic.es/footprintdb/ws.cgi');
}

if(defined($opts{'u'})){ $username = $opts{'u'} }





my ($result,$sequence,$sequence_name,$datatype,$keyword) = ('','','','','');


## sample protein sequence
$sequence_name = 'test';
$sequence = 'IYNLSRRFAQRGFSPREFRLTMTRGDIGNYLGLTVETISRLLGRFQKSGMLAVKGKYITIEN';

$result = $server->protein_query($sequence_name,$sequence,$footprintDBusername);
unless($result->fault()){
	print $result->result(); 
}else{
	print 'error: ' . join(', ',$result->faultcode(),$result->faultstring());
}

## sample regulatory motif sequence
#$sequence = 'TGTGANNN'; # possible format
#$sequence = "TGTGA\nTGTGG\nTGTAG"; # another format
# transfac format for position weight matrices

$sequence= <<EOM;
DE 1a0a_AB
01 1 93 0 2
02 0 96 0 0
03 58 33 3 2
04 8 78 6 4
05 8 5 75 8
06 1 2 47 46
07 1 2 84 9
XX
EOM

$result = $server->DNA_motif_query($sequence_name,$sequence,$footprintDBusername);
unless($result->fault()){
	print $result->result();
}else{
	print 'error: ' . join(', ',$result->faultcode(),$result->faultstring());
} 

$keyword = "myb";
$datatype = "site";
$result = $server->text_query($keyword,$datatype,$footprintDBusername);
unless($result->fault()){
	print $result->result();
}else{
	print 'error: ' . join(', ',$result->faultcode(),$result->faultstring());
}
