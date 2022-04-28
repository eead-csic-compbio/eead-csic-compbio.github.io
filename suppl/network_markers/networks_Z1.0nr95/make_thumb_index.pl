#!/usr/bin/perl -w

use strict;
use CGI;

my $CONVERTEXE = "/usr/bin/convert ";
my $THUMBSIZE = "200x200";

my $DIR = './';
my $TITLE = "networks_Z1.0nr95";
my $BUNDLE = "all.tgz";
my $REFERENCEINFO = "<pre>This is the set of transcriptional subnetworks for each non-redundant condition. You can download them all <a href=\"$BUNDLE\">here</a>. Please refer to ... for details.</pre>";

################################################

my ($n_of_thumbs,$pngthumb,$report) = (0);

# start HTML headers
$report .= "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n";
$report .= "<html>\n<head>\n<meta content=\"text/html\" http-equiv=\"content-type\">\n";
$report .= "<title>$TITLE</title></head><body>";
$report .= "<h2>$TITLE</h2>$REFERENCEINFO\n";
$report .= "<table border=\"0\" cellspacing=\"3\" cellpadding=\"3\"><tr>\n";

# read pdfs
opendir(DIR,$DIR) || die "# $0 : cannot list $DIR\n";
my @files = grep { /\.pdf/ && -f "$DIR/$_" } readdir(DIR);
closedir(DIR); 

# convert thumbs
foreach my $file (sort (@files))
{
	
	if($file =~ /(\S+)\.pdf/){ $pngthumb = $1 . '.png'; }
	system("$CONVERTEXE $file -resize $THUMBSIZE $pngthumb");
	
	$report .= "<td>$file<br><a href=\"$file\"><img src=\"$pngthumb\" alt=\"thumb\"></a></td>\n";
	
	$n_of_thumbs++;
	
	if($n_of_thumbs==5){ $report .= "</tr>\n<tr>"; $n_of_thumbs=0; }
}

if($n_of_thumbs && $n_of_thumbs < 5){ $report .= "</tr>"; }

$report .= "</table></body></html>";

print $report;
