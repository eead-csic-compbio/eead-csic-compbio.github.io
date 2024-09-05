---
title: "Computational & structural biology software"
---

<table align="center" width=100%>
  <tr>
    <td align="center"><b>[Home](index.html)</b>&nbsp;</td>
    <td align="center"><b>[Members](staff.html)</b>&nbsp;</td>
    <td align="center"><b>[Publications](publications.html)</b>&nbsp;</td>
    <td align="center"><b>[Software](software.html)</b>&nbsp;</td>
    <td align="center"><b>[Material educativo](matdidactico.html)</b>&nbsp;</td>
    <td align="center"><a href="https://bioinfoperl.blogspot.com"><b>Blog</b></a>&nbsp;</td>
    <td align="center"><a href="https://www.eead.csic.es"><img src="pics/logoEEAD.png"></a></td>
  </tr>
</table>


### Web resources 

* [RSAT::Plants](http://plants.rsat.eu): server for the analysis of regulatory sequences in plant genomes (available also as Docker container).
* [BARLEYMAP](https://barleymap.eead.csic.es): a tool to
    search the position of genetic markers on barley genomic, physical and POPSEQ maps.
* [barley_pangenes](https://eead-csic-compbio.github.io/barley_pangenes): clusters of gene models/alleles in a similar genomic location, linked from [BARLEYMAP](https://barleymap.eead.csic.es).
* [PRUNUSMAP](https://prunusmap.eead.csic.es): a tool to
    search the position of genetic markers and protein on annotated *Prunus* genomes.
* [footprintDB](https://footprintdb.eead.csic.es): a database
    of transcription factors with annotated cis elements and binding interfaces. Includes our own databases
  * [3D-footprint](https://3dfootprint.eead.csic.es), which runs structural analysis of protein-DNA complexes from the [Protein Data Bank](https://www.rcsb.org). 
  * [EEADannot](https://github.com/eead-csic-compbio/EEADannot), with manually curated DNA motifs and cis regulatory sites, mostly from plants.

    
### Source code, binaries and Docker containers

-   [plant-scripts](https://github.com/Ensembl/plant-scripts): code examples for interrogating 
[Ensembl Plants](https://plants.ensembl.org) from your own scripts, masking & annotating repeats and [calling pangenes](https://github.com/Ensembl/plant-scripts/tree/master/pangenes) in plant genomes ([GitHub downloads](https://tooomm.github.io/github-release-stats/?username=ensembl&repository=plant-scripts)).
-   [GET\_HOMOLOGUES](https://github.com/eead-csic-compbio/get_homologues):
    a versatile software package for gene-based pangenome analysis of microbes and plants (Linux/MacOSX
    Perl, R scripts, binaries, 
    [bioconda](https://anaconda.org/bioconda/get_homologues),
    [**[Docker]**](https://hub.docker.com/r/csicunam/get_homologues),
    [manuals and tutorials](https://github.com/eead-csic-compbio/get_homologues?tab=readme-ov-file#documentation),
    [GitHub downloads](https://tooomm.github.io/github-release-stats/?username=eead-csic-compbio&repository=get_homologues)).
-   [GET\_PHYLOMARKERS](https://github.com/vinuesa/get_phylomarkers): software package designed to identify optimal genomic markers for phylogenomics, population genetics and genomic taxonomy.

-   [**[RSAT Docker]**](https://hub.docker.com/r/biocontainers/rsat/tags): ready-to-use, requires downloading/installing genomes (see [docs](https://rsa-tools.github.io/installing-RSAT/RSAT-Docker/RSAT-Docker-tuto.html) and [protocol for coexpression motif discovery in plants](https://eead-csic-compbio.github.io/coexpression_motif_discovery)).
-   [Multigenomic Entropy-Based
    Score (MEBS)](https://github.com/eead-csic-compbio/metagenome_Pfam_score):
    protocol for finding informative protein families and then using
    them to score metagenomic sets.
-   [Chloroplast assembly
    protocol](https://github.com/eead-csic-compbio/chloroplast_assembly_protocol):
    set of scripts for the assembly of chloroplast genomes out of
    whole-genome sequencing reads.
-   [DNAPROT](https://hub.docker.com/r/eeadcsiccompbio/dnaprot): takes protein-DNA complex in PDB
    format and calculates structure-based position weight matrices ([manual](suppl/manual_dnaprot.pdf)), **[legacy Docker]**.
-   [primers4clades](https://hub.docker.com/r/csicunam/primers4clades): PCR primers for cross-species amplification of sequences from metagenomic DNA or selected lineages **[legacy Docker]**.
-   [TFcompare](https://hub.docker.com/r/eeadcsiccompbio/tfcompare): a tool for
    structural alignment of DNA motifs and protein domains from DNA-binding protein complexes **[legacy Docker]**.
-   [TFmodeller](https://hub.docker.com/r/eeadcsiccompbio/tfmodeller): comparative modelling of protein-DNA complexes **[legacy Docker]**.
    

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
-   see [all scripts](https://github.com/eead-csic-compbio/eead-csic-compbio.github.io/tree/master/scripts)

### Data files

Available [here](https://github.com/eead-csic-compbio/eead-csic-compbio.github.io/tree/master/data)


### GitHub repositories

Main: https://github.com/eead-csic-compbio

Other repos with code contributed by members of the lab:

* https://github.com/brunocontrerasmoreira
* https://github.com/Cantalapiedra 
* https://github.com/valdeanda/mebs
* https://github.com/vinuesa/get_phylomarkers 
* https://github.com/rsa-tools
* https://github.com/Ensembl/plant-scripts

