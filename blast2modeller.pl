#!/usr/bin/env perl
use strict;
use warnings;
use File::Temp qw/ tempfile /;

# search for a PDB template online and build a MODELLER model with it
# Bruno Contreras Moreira, Pablo Vinuesa 2019

my $DEFBLASTPEVAL  = 0.00001; # PDB contains ~100K structures
my $BLASTPRESTURL  = "http://www.rcsb.org/pdb/rest/getBlastPDB1?eCutOff=$DEFBLASTPEVAL&matrix=BLOSUM62&outputFormat=XML";
my $PDBDOWNLOADURL = "https://files.rcsb.org/download/WXYZ.pdb.gz";
my $BLASTPEXE      = "blastp";  # add path if needed
my $MODELLEREXE    = 'mod9.21'; # add path if needed
my $PDBTEMPLATEDIR = './'; 
my $DEFALIGNMETHOD = 'blastp';
my $CLEANTEMPFILES = 1;
my $MAXMODELS      = 1;
my %key_residues   = ();

if(!$ARGV[0]){
	print "# usage: $0 <input.faa> [comma-separated residue numbers]\n";
	print "# example: $0 my_sequence.faa 12,34,56\n";
	exit(0);
}
elsif($ARGV[1]){
	foreach my $res (split(/,/,$ARGV[1])){
		$key_residues{$res} = 1;
	}
}

## read first sequence in infile
my ($inFASTAfile,$query_seq_full,$query_name,$query_blast_file);
my $num_of_seqs = 0;
open(FAA,"<",$ARGV[0]) || die "# ERROR: cannot read $ARGV[0]\n";
while(<FAA>){
	chomp;
	if(/^>(\S+)/) {
		chomp;
		$query_name = $1;
		$query_blast_file = $query_name .'.blast';
		$query_blast_file =~ s/[\|;&><\[\]]/_/g;
		$num_of_seqs++;
		if($num_of_seqs > 1){
			warn "# WARNING: only first sequence in FASTA file will be processed\n";
		       	last;	
		}
	}
	else{ $query_seq_full .= $_; }
}
close(FAA);	

## run BLASTP remotely
if(!-e $query_blast_file){
	system("wget -q '$BLASTPRESTURL&sequence=$query_seq_full' -O $query_blast_file");
}

if(!-s $query_blast_file){
	die "# ERROR: cannot perform remote BLASTP to scan PDB templates, ".
		"please check network connection\n";
}

## parse BLASTP results
my ($PDBcode,$PDBchain,$raw_model_filename,$model_filename);
my ($bitscore,$Evalue,$from,$to,$identical,$perc_ident,$query_cover,$query_length);
my ($query_seq,$template_seq); # actually locally aligned sequence fragments 
my ($num_of_models,$num_of_highlighted,$DOPE) = 0;

# print column headers
print "PDB\tperc_ident\tperc_qcover\tbitscore\tEvalue\tqstart\tqend\thighlighted\tDOPE\tmodel_filename\n";

open(BLASTP,"<",$query_blast_file) || 
	die "# ERROR: could not run BLASTP remotely, please check your network connection\n";
while(<BLASTP>){
	if(/<BlastOutput_query-len>(\d+)</){ $query_length = $1 }
	elsif(/<Hit_def>(\w{4}):\d+:(\w+)/){ 
		($PDBcode,$PDBchain) = ($1,$2); 
		($query_seq,$template_seq) = ('',''); # init vars for this template
		($bitscore,$Evalue,$from,$to,$identical,$perc_ident,$query_cover) = 
			(0,$DEFBLASTPEVAL,0,0,0,0,0);
	}
	elsif(/<Hsp_bit-score>(\S+?)</){ $bitscore = $1 }
	elsif(/<Hsp_evalue>(\S+?)</){ $Evalue = $1 }
	elsif(/<Hsp_query-from>(\S+?)</){ $from = $1 }
	elsif(/<Hsp_query-to>(\S+?)</){ 
		$to = $1;
		$query_cover = sprintf("%1.1f",100*($to-$from+1)/$query_length);
	}
	elsif(/<Hsp_identity>(\S+?)</){ 
		$identical = $1;
		$perc_ident = sprintf("%1.1f",100*$identical/(($to-$from)+1));
	}
	elsif(/<Hsp_qseq>(\S+?)</){ $query_seq = $1 } # might contains gaps-
       	elsif(/<Hsp_hseq>(\S+?)</){ 
		$template_seq = $1; # last relevant piece of info for this template

		# process this template	
		$model_filename = 'NA';
		$DOPE = 9999.99;
		$num_of_highlighted = 0;
		$num_of_models++;

		# print details of aligned template, model file is added below
		print "$PDBcode\_$PDBchain\t$perc_ident\t$query_cover\t$bitscore\t$Evalue\t$from\t$to";

		if($num_of_models > $MAXMODELS){
			print "\t$num_of_highlighted\t$DOPE\t$model_filename\n";
			next;
		}

		my $pdbfilename = download_PDB_chain($PDBDOWNLOADURL,$PDBTEMPLATEDIR,$PDBcode,$PDBchain);
		if($pdbfilename eq 'download_error'){
			print "\t$num_of_highlighted\t$DOPE\tERROR: cannot download $PDBcode\n";
			next;
		}
		elsif($pdbfilename eq 'no_chain'){
			print "\t$num_of_highlighted\t$DOPE\tERROR: cannot find chain $PDBchain\n";
			next;
		}

		# get actual sequence contained in ATOM records of PDB file
		my ($atomic_protein_seq, $ref_residue_number2id) = extract_protein_sequence($pdbfilename);
	
		my ($align_file,$script_file) = make_infiles4modeller(
			$PDBTEMPLATEDIR,$PDBcode,$PDBchain,
			$atomic_protein_seq,$ref_residue_number2id,
			$query_name, $query_seq_full, $DEFALIGNMETHOD );
		
		($DOPE,$raw_model_filename) = run_modeller( $query_name, $script_file, $CLEANTEMPFILES);

		if($raw_model_filename =~ m/ERROR/){
			# terminate TSV string
			print "\t$num_of_highlighted\t$DOPE\t$raw_model_filename\n";
		}
		else{
			$model_filename = "$query_name.$DEFALIGNMETHOD.pdb";
			rename($raw_model_filename, $model_filename);

			if(keys(%key_residues)){
				$num_of_highlighted = 
					highlight_residues($model_filename, \%key_residues);
			}

			# add model name to TSV string
                        print "\t$num_of_highlighted\t$DOPE\t$model_filename\n";
		}

	} # elsif(/<Hsp_hseq>(\S+?)</){
}
close(BLASTP);

if($CLEANTEMPFILES){
	unlink( $query_name.'.blast' );
}

################################################################

# takes a PDB file and reference to a hash of residues ordinals
# which will be set with 99.99 B-factor compared to otherwise 
# values of 00.00 
sub highlight_residues {
	my ($PDB_filename, $ref_residues) = @_;

	my (@PDBcontents,$resid,$Bfactor,$atom);
	my %key_residues_found;

	open(PDBMODEL,"<",$PDB_filename) || die "# cannot read $PDB_filename\n";
	while(<PDBMODEL>){
		if(/^ATOM/){ 
                	$resid = substr($_,22,4);
			$resid =~ s/\s+//g;
			if($ref_residues->{$resid}){
				$Bfactor = ' 99.99';
				$key_residues_found{$resid}++;
			}
			else{
				$Bfactor = ' 00.00';
			}

			$atom = substr($_,0,60) . $Bfactor . substr($_,72);
			push(@PDBcontents,$atom);
		}
		else{ push(@PDBcontents,$_) }
	}
	close(PDBMODEL);
	
	if(keys(%key_residues_found)){
		my($fh, $filename) = tempfile( SUFFIX => '.pdb');
        	print $fh @PDBcontents;
        	close($fh);

		rename($filename, $PDB_filename);
	}

	return scalar(keys(%key_residues_found));
}	

# call modeller and clean temp files if requested
# returns a i) scalar with DOPE score and a ii) string which can 
# be a PDB filename or an ERROR message
# uses globals: MODELLEREXE
sub run_modeller {
	my ($query_name, $script_file, $clean_temp_files) = @_;

	my ($DOPE,$model_filename);
        open(MODELLER,"$MODELLEREXE $script_file 2>&1 |") ||
        	die "# ERROR: failed running $MODELLEREXE, check path\n";
        while(<MODELLER>){
        	if(/Top model: (\S+) DOPE: (\S+)/){
                	$model_filename = $1;
			$DOPE = $2;
                }
                elsif(/_modeller.ModellerError:/){
			$DOPE = 9999.99;
			$model_filename = "ERROR: $_";
                        last;
                }
	}
	close(MODELLER);

	# clean MODELLER tmp files
        if($clean_temp_files){
        	my @modeller_temp_files = glob( "$query_name.D* $query_name.V* $query_name*.log");
                unlink(@modeller_temp_files);
                unlink($query_name.'.ini', $query_name.'.rsr', $query_name.'.sch' );
        }

	return ($DOPE,$model_filename);
}

# computes pairwise alignment with selected method 
# returns: $qstart,$qend,$sstart,$send,$qseq,$sseq
# uses globals: $BLASTPEXE, 
sub pairwise_align {
	my ($queryseq,$templateseq,$method) = @_;

	my($fhq, $filenameq) = tempfile( SUFFIX => '.faa');
        print $fhq ">query\n$queryseq\n";
        close($fhq);

        my($fht, $filenamet) = tempfile( SUFFIX => '.faa');
        print $fht ">template\n$templateseq\n";
        close($fht);

        my ($qstart,$qend,$sstart,$send,$qseq,$sseq);

	if(!$method || $method eq 'blastp') {
		
        	open(BLAST,"$BLASTPEXE -query $filenameq -subject $filenamet -max_hsps 1 ".
                        "-use_sw_tback -outfmt \"6 qstart qend sstart send qseq sseq\" |") ||
                	die "# ERROR: cannot run $BLASTPEXE\n";
        	while(<BLAST>){
                	chomp;
                	($qstart,$qend,$sstart,$send,$qseq,$sseq) = split(/\h/,$_);
        	}
        	close(BLAST);
	}	

	return ($qstart,$qend,$sstart,$send,$qseq,$sseq);
}



# computes pairwise alignment and returns 
# i) filename of alignment formatted for modeller
# ii) filename of script for modeller
sub make_infiles4modeller {
	my ($tdir,$PDBcode,$chain,$PDBseq,$ref_tnum2id,$qname,$queryseq,$method) = @_;

	my $tname = "$PDBcode\_$chain";	

	my ($qstart,$qend,$sstart,$send,$qseq,$sseq) = pairwise_align($queryseq,$PDBseq,$method);

	# convert alignment coords of template to actual residue ids in ATOM records
	$sstart = $ref_tnum2id->{$sstart};
       	$send   = $ref_tnum2id->{$send};	

	# write alignment to file
	my $align_filename = "$qname.$tname.ali";
	open(ALI,">",$align_filename) || die "# ERROR: cannot create $align_filename\n";
	print ALI ">P1;$qname\nsequence:$qname:$qstart:A:$qend:A:query:::\n$qseq*\n";
	print ALI ">P1;$tname\nstructureX:$tname:$sstart:$chain:$send:$chain:template:::\n$sseq*\n"; 
	close(ALI);

	# write script file
	my $script_filename = "$qname.$tname.py";
        open(SCRIPT,">",$script_filename) || die "# ERROR: cannot create $script_filename\n";
	my $script = <<~MODSTRING;
		import sys
		from modeller.automodel import *   
		class MyModel(automodel):
			def special_patches(self, aln):
				self.rename_segments(segment_ids='A', renumber_residues=$qstart)

		log.verbose()    
		env = environ() 
		env.io.atom_files_directory = '$PDBTEMPLATEDIR'
		a = automodel(env,alnfile='$align_filename',
			knowns='$tname',sequence='$qname',
			assess_methods=(assess.DOPE))      
		a.starting_model= 1 
		a.ending_model  = 1                 
		a.make()
		ok_models = [x for x in a.outputs if x['failure'] is None]
		m = ok_models[0]
		sys.stderr.write("Top model: %s DOPE: %.3f" % (m['name'], m['DOPE score']))

		MODSTRING
	print SCRIPT $script;
        close(SCRIPT);

	return ($align_filename,$script_filename);
}

# takes input PDB plain filename 
# returns i) string with one-letter amino acid sequence actually contained in ATOM records
# and ii) ref to hash with residue ordinal to actual residue ID in PDB 
sub extract_protein_sequence {
	my ($PDB_chain_file) = @_;

	my %aa3to1;
        $aa3to1{"ALA"} = "A"; $aa3to1{"CYS"} = "C"; $aa3to1{"ASP"} = "D"; $aa3to1{"GLU"} = "E";
        $aa3to1{"PHE"} = "F"; $aa3to1{"GLY"} = "G"; $aa3to1{"HIS"} = "H"; $aa3to1{"ILE"} = "I";
        $aa3to1{"LYS"} = "K"; $aa3to1{"LEU"} = "L"; $aa3to1{"MET"} = "M"; $aa3to1{"ASN"} = "N";
        $aa3to1{"PRO"} = "P"; $aa3to1{"GLN"} = "Q"; $aa3to1{"ARG"} = "R"; $aa3to1{"SER"} = "S";
        $aa3to1{"THR"} = "T"; $aa3to1{"VAL"} = "V"; $aa3to1{"TRP"} = "W"; $aa3to1{"TYR"} = "Y";
        $aa3to1{"MSE"} = "M";

	my ($atom_record_sequence,%number2id,$atom_name,$res_name,$resid);
	my $resnumber = 0;
        open(PDB,"<",$PDB_chain_file) || die "# ERROR: cannot read $PDB_chain_file\n";
        while(<PDB>){
		$atom_name = substr($_,12,4);
                if($atom_name eq ' CA '){
			$resnumber++;
	                $res_name = substr($_,17,3);
			$resid    = substr($_,22,4);
                        $atom_record_sequence .= $aa3to1{$res_name};
			$number2id{$resnumber} = $resid;
                }
        }
        close(PDB);

	return ($atom_record_sequence, \%number2id);
}

# returns name of uncompressed PDB file in $dir or 
# empty string in case of failure
sub download_PDB_chain {
	my ($URL,$dir,$PDB_code,$PDB_chain) = @_;

	my $fullURL = $URL;
	my $PDB_code_uc = uc($PDB_code);
	$fullURL =~ s/WXYZ/$PDB_code_uc/;

	my $tmp_pdb_file_gz = $PDB_code . '.pdb.gz';
	my $pdb_chain_file_name = "$dir/$PDB_code\_$PDB_chain.pdb";

	if(-s $pdb_chain_file_name){
		return $pdb_chain_file_name;
	}

	# download compressed, full PDB entry
	system("wget -q $fullURL -O $tmp_pdb_file_gz");

	# parse selected chain
	my (@chain_contents,$alt);
	open(PDBGZ,"zcat $tmp_pdb_file_gz |") || return "download_error";
	while(<PDBGZ>){

		last if(/^ENDMDL/); # NMR structures
                if(/^HETATM/ && /MSE/){ $_ =~ s/^HETATM/ATOM  / }
                next if not (/^ATOM/ || length($_) < 21);

		# skip alternative atom positions
                $alt = substr($_,16,1);
                next if($alt ne ' ' && $alt ne 'A');

		if(substr($_,21,1) eq $PDB_chain) {
			push(@chain_contents,$_);
		}
	}
	close(PDBGZ);
	unlink($tmp_pdb_file_gz);

	# save coords in plain file
	if(@chain_contents){
		open(CHAINFILE,">",$pdb_chain_file_name); 
		foreach my $atom (@chain_contents){ 
			print CHAINFILE $atom; 
		}
		close(CHAINFILE);
	}
	else{ return "no_chain" }

	return ($pdb_chain_file_name);
}
