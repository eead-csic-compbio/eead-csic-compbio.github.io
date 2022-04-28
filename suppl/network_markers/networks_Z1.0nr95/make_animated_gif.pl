#!/usr/bin/perl -w

use strict;
use CGI;

my $CONVERTEXE = "/usr/bin/convert ";
my $THUMBSIZE = "400x400";
my $DELAY     = 5; # half a second

my $DIR = './';
my $GIFdir  = $DIR . 'animation/'; 
my $GIFfile = "TRN_animation.gif";

################################################

my ($n_of_gifs,$command,$gif) = (0,'');

# read pdfs
opendir(DIR,$DIR) || die "# $0 : cannot list $DIR\n";
my @files = grep { /\.pdf/ && -f "$DIR/$_" } readdir(DIR);
closedir(DIR); 

# convert to GIFs
foreach my $file (sort (@files))
{
	
	if($file =~ /(\S+)\.pdf/){ $gif = $GIFdir . $1 . '.gif'; }
	system("$CONVERTEXE $file -resize $THUMBSIZE $gif");
	
	$n_of_gifs++;
}

# create animated GIF
system("$CONVERTEXE -delay $DELAY -loop 0 $GIFdir/*.gif $GIFfile");
