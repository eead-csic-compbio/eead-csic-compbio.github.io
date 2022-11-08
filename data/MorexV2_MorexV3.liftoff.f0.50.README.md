
TSV file [MorexV2_MorexV3.liftoff.f0.50.tsv](./MorexV2_MorexV3.liftoff.f0.50.tsv) contains mappings of gene IDs from barley MorexV2/TRITEX  
to MorexV3 gene models. These were produced by:
+ Lifting over v2 genes with [Liftoff](https://github.com/agshumate/Liftoff) using the 
v2 and v3 toplevel genomic sequences at [Ensembl Plants](https://plants.ensembl.org/index.html) release 52. 


      liftoff -g Barley_Morex_V2_gene_annotation_PGSB.all.gff3 -o Barley_Morex_V2.PGSB.all.liftoff.gff GCA_904849725.1_MorexV3_pseudomolecules.chrnames.fna Barley_Morex_V2_pseudomolecules.fasta

+ Intersecting the lifted over v2 genes and the v3 genes with [bedtools](https://bedtools.readthedocs.io) requiring 50% overlap.

      bedtools intersect -sorted -b Barley_Morex_V2.PGSB.allgenes.liftoff.gff -a ../../gene_annotation/Hv_Morex.pgsb.Jul2020.allgenes.gff -wo -f 0.5

The column names of the TSV file are:

IBSCv2 gene, MorexV3 gene, v2 chr, v2 start v2 end, v3 chr, v3 start, v3 end

Each row is a gene model pair, both HC and LC models are included.
