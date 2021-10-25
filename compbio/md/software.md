
# Computational & structural biology group


<table align="center" width=100%>
  <tr>
    <td align="center"><b>[Home](index.html)</b></td>
    <td align="center"><b>[Members](staff.html)</b></td>
    <td align="center"><b>[Publications](publications.html)</b></td>
    <td align="center"><b>[Software](software.html)</b></td>
    <td align="center"><b>[Material educativo](matdidactico.html)</b></td>
    <td align="center"><a href="http://bioinfoperl.blogspot.com"><b>Blog</b></a></td>
    <td align="center"><a href="http://www.eead.csic.es"><img src="pics/logoEEAD.jpeg"></a></td>
  </tr>
</table>


### Web resources 

-   [RSAT::Plants](http://plants.rsat.eu): A plant-dedicated server for
    the analysis of regulatory sequences
-   [BARLEYMAP](http://floresta.eead.csic.es/barleymap/): a tool to
    search the position of barley genetic markers on the Barley Physical
    Map and the POPSEQ map
-   [footprintDB](http://floresta.eead.csic.es/footprintdb/): a database
    of transcription factors with annotated cis elements and binding
    interfaces
-   [TFcompare](http://floresta.eead.csic.es/tfcompare/): a tool for
    structural alignment of DNA motifs and protein domains from DNA
    binding protein complexes
-   [3D-footprint](http://floresta.eead.csic.es/3dfootprint/): a
    database for the structural analysis of protein-DNA complexes
<!---   [TFmodeller](http://www.ccg.unam.mx/tfmodeller/): a server to build and analyse comparative models of DNA-binding proteins -->
-   [primers4clades](http://floresta.eead.csic.es/primers4clades/):
    server for the design of primers for phylogenetic analyses based on
    PCR amplifications, created in collaboration with [Pablo
    Vinuesa](http://www.ccg.unam.mx/%7Evinuesa/) (UNAM,México)

### Source code and binaries 

-   [GET\_HOMOLOGUES](https://github.com/eead-csic-compbio/get_homologues):
    a versatile software package for pan-genome analysis (Linux/MacOSX
    Perl, R scripts, binaries, [Docker
    image](https://hub.docker.com/r/csicunam/get_homologues),
    [manual](http://eead-csic-compbio.github.io/get_homologues/manual/),
    [microbia](http://aem.asm.org/cgi/pmidlookup?view=long&pmid=24096415),
    [plants](http://journal.frontiersin.org/article/10.3389/fpls.2017.00184/full))
-   [GET\_PHYLOMARKERS](https://github.com/vinuesa/get_phylomarkers): software package designed to identify optimal genomic markers for phylogenomics, population genetics and genomic taxonomy.
-   [coexpression_motif_discovery](https://eead-csic-compbio.github.io/coexpression_motif_discovery): recipes to discover cis-regulatory motifs within proximal promoters of plants in RSAT::Plants.
-   [Multigenomic Entropy-Based
    Score (MEBS)](https://github.com/eead-csic-compbio/metagenome_Pfam_score):
    protocol for finding informative protein families and then using
    them to score metagenomic sets
    ([PDF](https://academic.oup.com/gigascience/advance-article/doi/10.1093/gigascience/gix096/4561660))
-   [Chloroplast assembly
    protocol](https://github.com/eead-csic-compbio/chloroplast_assembly_protocol):
    set of scripts for the assembly of chloroplast genomes out of
    whole-genome sequencing reads
-   [DNAPROT](./soft/dnaprot.php): takes protein-DNA complex in PDB
    format and calculates structure-based position weight matrices
    (Linux 64bit Perl scripts, binaries,
    [manual](http://www.eead.csic.es/compbio/soft/manual_dnaprot.pdf),
    [PDF](http://www.biomedcentral.com/1471-2105/9/436))

### Assorted utilities 

-   [split\_pairs](https://github.com/eead-csic-compbio/split_pairs):
    efficient
    [kseq](http://lh3lh3.users.sourceforge.net/kseq.shtml)-based program
    to sort and find paired reads within FASTQ/FASTA files, with the
    ability to edit headers with the power of Perl-style regular
    expressions
-   [split\_blast](http://bioinfoperl.blogspot.com.es/2013/04/splitblastpl-real-multicore-blast.html): Perl script to take advantage of multi-core CPUs for doing BLAST searches that fit in RAM (also part of GET\_HOMOLOGUES)
-   [addCDD2genbank.pl](https://github.com/eead-csic-compbio/eead-csic-compbio.github.io/blob/master/scripts/addCDD2genbank.pl): adds domain annotations from CDD to protein sequences contained in CDS features within input GenBank file
-   [xmfa2fasta.pl](https://github.com/eead-csic-compbio/eead-csic-compbio.github.io/blob/master/scripts/xmfa2fasta.pl):    reads in XMFA file produced by progressive MAUVE and produces a multi-FASTA file containing a multiple sequence alignment
-   [barleyGO](https://github.com/eead-csic-compbio/barleyGO): annotates barley sequences (Perl scripts, C++ source, Linux 64bit binaries)

### GitHub repositories

Main: https://github.com/eead-csic-compbio

Other repos with code contributed by members of the lab:

* https://github.com/Cantalapiedra 
* https://github.com/vinuesa/get_phylomarkers 
* https://github.com/Ensembl/plant-scripts 

