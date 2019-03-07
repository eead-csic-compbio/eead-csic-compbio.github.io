
import logging
import requests, sys
import argparse


FORMAT = '%(asctime)s - %(name)s(%(lineno)d) - %(levelname)s - %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger('ensembl_position_from_snp')
logger.setLevel('DEBUG')
#logger.setLevel('INFO')
logger.info('Starting')

# ------------------------------------------------------
# INPUT :  A file with a list of variant of interest (first 9 columns of a VCF, example: "chr7A	636901870	BA00796146;AX-94381170	G	C	.	.	.	GT")
# OUTPUT : List of affected gene IDS plus a summary of the data pulled from the ensemble genome API
# ------------------------------------------------------


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', help="A file with a list of variant of interest (first 9 columns of a VCF, example: 'chr7A	636901870	BA00796146;AX-94381170	G	C	.	.	.	GT')", type=str, required=True)
args = parser.parse_args()
_INPUT_SNP_FILE = args.input_file



# parse input snps
snps_regions = []
with open(_INPUT_SNP_FILE, 'r') as f:
    for line in f:
        if not line.startswith("#"):
            snps_regions.append(line.replace("\t"," ").replace("\n"," "))

logger.debug(str(snps_regions))

server = "http://rest.ensemblgenomes.org"
ext = "/vep/triticum_aestivum/region"
headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

fh = open("affectedGenes.tsv", "w")



for variant in snps_regions:
    data_string='{ "variants" : ["'+variant+'" ] }'
    logger.debug(data_string)
    r = requests.post(server+ext, headers=headers, data=data_string)

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    print(repr(decoded))

    if "transcript_consequences" in decoded[0] :
        colocated_variants = ""
        if 'colocated_variants' in decoded[0]:
            for col_var in decoded[0]['colocated_variants']:
                colocated_variants += col_var['id']+":"+decoded[0]['seq_region_name']+":"+col_var['start']+":"+col_var['end']
        for my_consequence in decoded[0]['transcript_consequences']:
            affected_gene= my_consequence['gene_id'] + \
                           "\tcolocated_variants:" + colocated_variants+ \
                           "\tconsequence_terms:" + ",".join(my_consequence['consequence_terms'])+\
                           "\ttranscript_id:"+ my_consequence['transcript_id']+ \
                           "\timpact:"+ my_consequence['impact']
            if "distance" in my_consequence:
                affected_gene+="\tdistance:" + str(my_consequence['distance'])+"\n"
            else:
                affected_gene+="\n"
            fh.write (str(affected_gene))
            logger.debug(affected_gene)
    # No Gene ID with intergenic consequence , TODO: check this
    # if 'intergenic_consequences' in decoded[0]:
    #     for my_consequence in decoded[0]['intergenic_consequences']:
    #         affected_gene= my_consequence['gene_id'] + \
    #                        "\tcolocated_variants:" + colocated_variants+ \
    #                        "\tconsequence_terms:" + ",".join(my_consequence['consequence_terms'])+ \
    #                        "\timpact:"+ my_consequence['impact']+"\n"
    #         fh.write (str(affected_gene))
    #         logger.debug(affected_gene)

fh.close()

exit(0)