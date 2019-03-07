#!/usr/bin/env/python

#bring python3 print() into python2
from __future__ import print_function
import os, math, traceback, datetime, urllib, requests, json, argparse, sys

pheno = ["yield", "grain"]

def getgs(genedesign):
    '''Open gene list and use them to search Knetminer along with keywords'''
    with open(genedesign, "r") as gk:
        with open("genome.json", "w") as af:
            genes = []
            for line in gk:
                col = line.split("\t")
                genes.append(col[0])
            genelist = (",").join(genes) #join all iterative elements by ,
            #print(genelist)
            #pheno = ["yield", "grain size"]
            #use str.join() to convert multiple elments in a list into one string.
            keyw = "%20OR%20".join(pheno)
            url = "http://knetminer.rothamsted.ac.uk/wheatknet/genome?keyword={}&list={}".format(keyw, genelist)
            #print(url)
            r = requests.get(url)
            r.json()
            r.status_code #check if request is successful.
            print(r.text, file=af)
        af.close()
    gk.close()
    return

def parsejs(genetable):
    ''' deserialise json into dictionary and extract the genetable which hopefully provide right genes and score given right url'''
    with open("genome.json", "r") as jf:
        content = json.load(jf) #deserialise content of json, which will be dictionary object.
        #print(type(content))
        with open(genetable, "w") as g:
            print(content[u'geneTable'], file=g) #r.json will put a u infront of the keys of json dictionary
        g.close()
    jf.close()
    return

def gene_score(genetable, scores):
    '''Extract the scores only.'''
    with open(genetable, "r") as f:
        next(f)
        with open(scores, "w") as sf:
            for line in f:
                if(line == "\n"): continue 
                col = line.split("\t")
                #print (str(col))
                score=str(col[6])
                gene_id=col[1]
                gene_name=col[2]
                #pheno = ['yield','grain']
                keyw = "%20OR%20".join(pheno)
                parameters = {"keyword":keyw, "list":gene_id}
                link="http://knetminer.rothamsted.ac.uk/wheatknet/genepage?"
                r=requests.get(link, params=parameters)
                print("{}\t{}\t{}\t{}".format(gene_id, gene_name, score, r.url), file=sf)
        sf.close()
    f.close()
#End of block 7 to print genes, scores and url into scores.txt

if __name__ == "__main__":
    print("KnetMiner API retrieval script")
    #genedesign = "input.txt"
    input_data = sys.argv[1]
    try:
        getgs(input_data)
        print("Query KnetMiner genome API...")
    except Exception:
        traceback.print_exc()
    
    
    genetable="genetable.tab"
    try:
        parsejs(genetable)
        print("Processing KnetMiner output...")
    except Exception:
        traceback.print_exc()

    # Compare input gene list with KnetMiner output
    genetable="genetable.tab"
    scores="scores.tab"
    try:
        gene_score(genetable, scores)
        print("Create output...")
    except Exception:
        traceback.print_exc()


print(datetime.datetime.now())