
TSV file [IBSCv1_MorexV3.barleymap.tsv](./IBSCv1_MorexV3.barleymap.tsv) contains mappings of gene IDs from barley IBSCv1 (2012)
to MorexV3 gene models. These were produced by mapping v1 transcript isoforms with [Barleymap](https://floresta.eead.csic.es/barleymap/) to the MorexV3 sequence at [Ensembl Plants](https://plants.ensembl.org/index.html) release 52.

    bmap_align --maps=MorexV3 -g HC_complete_genes.short_header.fna -b > HC2012.bmap 
    bmap_align --maps=MorexV3 -g LC_complete_genes.short_header.fna -b > LC2012.bmap

The column names of the TSV file are:

IBSCv1 gene, MorexV3 gene, v2 chr, v2 start v2 end, v3 chr, v3 start, v3 end

Each row is a gene model pair, with HC models first and LC models towards the end of the file.
