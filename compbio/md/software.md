\

laboratory of computational & structural biology
------------------------------------------------

<div id="header">

  ----------------------------------- ------------------------------------------------------------------------ ---------------------------------------- ------------------------------------
  ![logo CSIC](./pics/logoCSIC.png)   [Estación Experimental Aula Dei](http://www.eead.csic.es/)\              [Fundación ARAID](http://www.araid.es)   ![logo ARID](./pics/logoARAID.gif)
                                      [Consejo Superior de Investigaciones Científicas](http://www.csic.es/)                                            
  ----------------------------------- ------------------------------------------------------------------------ ---------------------------------------- ------------------------------------

</div>

<div id="content">

-   [Home](index.html)
-   [Lab members](staff.html)
-   [Our collaborators](collaborators.html)
-   [Publications](publications_computational_biology_bioinformatics.html)
-   [Software](software_computational_biology_bioinformatics.html)
-   [Material
    educativo](material_didactico_biologia_computacional_bioinformatica.html)

\
### Web resources {#web-resources style="text-align:left"}

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
-   [TFmodeller](http://www.ccg.unam.mx/tfmodeller/): a server to build
    and analyse comparative models of DNA-binding proteins
-   [primers4clades](http://floresta.eead.csic.es/primers4clades/):
    server for the design of primers for phylogenetic analyses based on
    PCR amplifications, created in collaboration with [Pablo
    Vinuesa](http://www.ccg.unam.mx/%7Evinuesa/) (UNAM,México)

### Code and binaries {#code-and-binaries style="text-align:left"}

-   [GET\_HOMOLOGUES](https://github.com/eead-csic-compbio/get_homologues):
    a versatile software package for pan-genome analysis (Linux/MacOSX
    Perl, R scripts, binaries, [Docker
    image](https://hub.docker.com/r/csicunam/get_homologues),
    [manual](http://eead-csic-compbio.github.io/get_homologues/manual/),
    [microbia](http://aem.asm.org/cgi/pmidlookup?view=long&pmid=24096415),
    [plants](http://journal.frontiersin.org/article/10.3389/fpls.2017.00184/full))
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

### Software utilities {#software-utilities style="text-align:left"}

-   [split\_pairs](https://github.com/eead-csic-compbio/split_pairs):
    efficient
    [kseq](http://lh3lh3.users.sourceforge.net/kseq.shtml)-based program
    to sort and find paired reads within FASTQ/FASTA files, with the
    ability to edit headers with the power of Perl-style regular
    expressions
-   [split\_blast](http://bioinfoperl.blogspot.com.es/2013/04/splitblastpl-real-multicore-blast.html):
    Perl script to take advantage of multi-core CPUs for doing BLAST
    searches that fit in RAM (also part of GET\_HOMOLOGUES)
-   [addCDD2genbank.pl](https://github.com/eead-csic-compbio/eead-csic-compbio.github.io/blob/master/addCDD2genbank.pl)
    adds domain annotations from CDD to protein sequences contained in
    CDS features within input GenBank file
-   [xmfa2fasta.pl](https://github.com/eead-csic-compbio/eead-csic-compbio.github.io/blob/master/xmfa2fasta.pl)
    reads in XMFA file produced by progressive MAUVE and produces a
    multi-FASTA file containing a multiple sequence alignment
-   [barleyGO](https://github.com/eead-csic-compbio/barleyGO): annotates
    barley sequences (Perl scripts, C++ source, Linux 64bit binaries)
-   [FRAGBENCH](./soft/fragbench.tgz): estimates the accuracy limits for
    template-based predictions of protein structure and was tested with
    CASP5 data (Perl scripts, [PDF](./papers/fragbench_2005.pdf))
-   [MSUPER](./soft/msuper.tgz): calculates rigid-body progressive
    multiple alignments of protein structures (Linux 64bit binary,
    algorithm published as an appendix in the [PhD
    thesis](./papers/PhD_thesis_Bruno_Contreras-Moreira2003.pdf)
    of B.Contreras-Moreira)

GitHub repository of the lab: <https://github.com/eead-csic-compbio>\
\
Other repos with code contributed by members of the lab:\
<https://github.com/Cantalapiedra>\
<https://github.com/vinuesa/get_phylomarkers>\
\
### Other web projects and curated-databases in which we have been involved: {#other-web-projects-and-curated-databases-in-which-we-have-been-involved style="text-align:left"}

-   DomainFishing: splits protein sequences in domains and looks for
    structural templates to model them
-   3D-JIGSAW: a comparative modelling web server that integrated
    DomainFishing
-   <span style="font-style:italic">in silico</span> Protein
    Recombination: web server that combines protein models with a
    genetic algorithm, benchmarked in the lab and during
    [CASP5](http://predictioncenter.org/casp5/Casp5.html)
-   [RegulonDB](http://regulondb.ccg.unam.mx) and
    [EcoCyc](http://ecocyc.org): curated databases of *Escherichia coli*
    regulation and metabolism

</div>
