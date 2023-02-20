
The pangene matrices in BED and TSV format describe clusters of gene models 
found to be collinear across barley assemblies and annotations in the following 
order: 

    MorexV3 Morex Barke HOR10350 HOR3081 HOR9043 OUN333 HOR7552 Igri HOR21599 HOR13942 Akashinriki HOR8148 RGT_Planet HOR13821 B1K-04-12 HOR3365 Hockett ZDM02064 ZDM01467 Golden_Promise

The gene annotations come from https://doi.org/10.5447/ipk/2020/24, 
https://doi.org/10.1111/tpj.15871 (BaRTv2)
and https://plants.ensembl.org (MorexV3).
 
The **TSV** file contains one cluster per row, the 1st column being the cluster name.

The **BED** file describes the same clusters with the first four columns indicating the
position of the relevant gene model in the MorexV3 reference and the cluster name.
Clusters that do not include MorexV3 genes do not have an exact position, they are placed
within an interval defined by two MorexV3 genes. Those rows start with '#'.

These were produced by clustering collinear genes with 
[get_pangenes.pl](https://github.com/Ensembl/plant-scripts/tree/master/pangenes) with the following arguments:

    plant-scripts/pangenes/get_pangenes.pl -d barley -s '^chr\d+H' -m cluster -r MorexV3 -H -t 0 &> log.barley.H.t0.MorexV3.txt

The log file is available [here](./MorexV3_highrep_0taxa_5neigh_algMmap_split_/log.barley.H.t0).

This protocol has been preprinted at https://www.biorxiv.org/content/10.1101/2023.01.03.520531v1.
