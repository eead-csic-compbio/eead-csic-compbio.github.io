
TSV file [IBSCv2_MorexV3.barleymap.tsv](./IBSCv2_MorexV3.barleymap.tsv) contains mappings of gene IDs from barley IBSCv2  
to MorexV3 gene models. These were produced by:
+ Mapping v2 longest transcript isoforms with [Barleymap](https://floresta.eead.csic.es/barleymap/) to the MorexV3 
sequence at [Ensembl Plants](https://plants.ensembl.org/index.html) release 52.

    bmap_align --maps=MorexV3 -g HC2017.fna -u > HC2017.bmap 
    bmap_align --maps=MorexV3 -g LC2017.fna -u > HC2017.bmap
