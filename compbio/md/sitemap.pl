#!/usr/bin/perl
# based on https://gist.github.com/pdp7/9784888

use strict;
use warnings;
use POSIX qw(strftime);

my $BASEURL = 'https://www.eead.csic.es/compbio/';
my @htmlfiles = glob("*.html");

print <<EOF;
<?xml version="1.0" encoding="UTF-8"?>
<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">
EOF

foreach my $file (@htmlfiles) {
  
    my $url = $BASEURL . $file; 
    
    # check modification time
    my $last_mod_time = strftime("%Y-%m-%d",
      localtime( (stat($file))[9] ));
    
    # take care of special chars
    $url =~ s/\s+$//;
    $url =~ s/^\s+//;
    $url =~ s/&/&amp;/;
    
    print <<EOF;
    <url>
        <loc>$url</loc>
        <lastmod>$last_mod_time</lastmod>
    </url>
EOF
}

print <<EOF;
</urlset>
EOF