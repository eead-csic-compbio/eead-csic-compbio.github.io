#!/usr/bin/perl -w

# extract CDS sequences from full GenBank records downloaded from the
# NCBI site with queries such as:
# https://www.ncbi.nlm.nih.gov/nuccore/?term=%22Hordeum+vulgare%22%5BOrganism%5D+AND+ddbj_embl_genbank%5Bfilter%5D+AND+gene%5BTitle%5D

use Bio::SeqIO;

my $TAKEORGANELLES   = 0;
my $MAXLENGTHPRODUCT = 30;

die "usage: $0 <sequences.gb>\n" unless(@ARGV);
 
my ($gbfile) = @ARGV;

my ($dna_fasta_file,$n_of_CDSs) = extract_CDSs_from_genbank($gbfile);


sub extract_CDSs_from_genbank
{
  my ($infile, $max_records) = @_;
	
  my ($out_dna_file,$n_of_CDS) = ($infile,0);
  $out_dna_file =~ s/\.gb/_cds.fna/;
  
  my ($taxon,$cv,$gene,$gbaccession,$protid,$product);
  my ($CDScoords,$CDSseq,$element,$isOrganelle,$chr);

  my $in = new Bio::SeqIO(-file => $infile, -format => 'genbank' );
	
  open(DNA,">$out_dna_file") || die "# extract_CDSs_from_genbank : cannot create $out_dna_file\n";
	
  SEQ: while( my $seq = $in->next_seq) {    
    
    $gbaccession = $seq->accession(); #print "$gbaccession\n";
    ($taxon,$cv,$chr,$isOrganelle) = ('','','',0);

    FEAT: foreach my $f ($seq->get_SeqFeatures) {
    
      if($f->primary_tag() =~ /source/) {
        if($f->has_tag('organelle')){ $isOrganelle = 1 }
        
        if($f->has_tag('organism')) {
          foreach $element ($f->each_tag_value('organism')){ $taxon = $element; last } 
        }
      
        if($f->has_tag('cultivar')) {
          foreach $element ($f->each_tag_value('cultivar')){ $cv = " cv $element"; last }
        }

        if($f->has_tag('chromosome')) {
          foreach $element ($f->each_tag_value('chromosome')){ $chr = $element; last }
        }
      }
      elsif($f->primary_tag() =~ /CDS/) {
      
        $CDScoords = $f->spliced_seq(); 
        $CDSseq = '';
       
        if($f->location->isa('Bio::Location::SplitLocationI')) {

          for my $location ( $f->location->sub_Location() ) {
            $CDSseq .= $seq->subseq($location->start,$location->end);
            if($location->strand == -1){ $rev = 1 }
          }

          if($rev){
            $CDSseq =~ tr/acgtnACGTN/tgcanTGCAN/;
            $CDSseq = reverse($CDSseq);
          } 
        } else {
          $CDSseq = $CDScoords->{'seq'};
        }
                
        if($CDSseq eq ''){
          warn "# extract_CDSs_from_genbank: cannot get CDS sequence for $gbaccession\n";
          next;
        }
        
        ($gene,$protid,$product) = ('','','');
        if($f->has_tag('gene')) {
          $gene = join(',',sort $f->each_tag_value('gene'));
        } 
        if ($f->has_tag('protein_id')) {
          $protid = join(',',sort $f->each_tag_value('protein_id'));
        }
        if ($f->has_tag('product')) {
          $product = join(',',sort $f->each_tag_value('product'));
          if(length($product)>$MAXLENGTHPRODUCT){ $product = substr($product,0,$MAXLENGTHPRODUCT) }
        }
				
        if($isOrganelle && $TAKEORGANELLES==0) {
          print "# extract_CDSs_from_genbank: skip $gbaccession $gene (TAKEORGANELLES=$TAKEORGANELLES)\n";
        }
        else {
				  print DNA ">$protid|$gene|$product|$gbaccession|$chr|[$taxon$cv]\n$CDSseq\n";
          $n_of_CDS++;	
        
          last SEQ if(defined($max_records) && $n_of_CDS == $max_records);
        }
      } 
    }
  }

  close(DNA);
	
	return ($out_dna_file,$n_of_CDS);
}
