
TSV file [IBSCv2_MorexV3.liftoff.f0.50.tsv](./IBSCv2_MorexV3.liftoff.f0.50.tsv) contains mappings of gene IDs from barley IBSCv2  
to MorexV3 gene models. These were produced by:
+ Lifting over v2 genes with [Liftoff](https://github.com/agshumate/Liftoff) using the 
v2 and v3 toplevel genomic sequences at [Ensembl Plants](https://plants.ensembl.org/index.html) release 52. 

    liftoff -g Hordeum_vulgare.IBSCv2.gff3 -o Hordeum_vulgare.IBSCv2.liftoof.gff3 hordeum_vulgare.dna.fasta hordeum_vulgare.old.fasta

+ Intersecting the lifted over v2 genes and the v3 genes with [bedtools](https://bedtools.readthedocs.io) requiring 50% overlap.

    bedtools intersect -sorted -b Hordeum_vulgare.IBSC_v2.liftoof.gff -a Hv_Morex.pgsb.Jul2020.HC.gene.gff -wo -f 0.5
    bedtools intersect -sorted -b Hordeum_vulgare.IBSC_v2.liftoof.gff -a Hv_Morex.pgsb.Jul2020.LC.gene.gff -wo -f 0.5


