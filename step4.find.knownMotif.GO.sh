##############校正后的P值小于0.05的motif对应的TF用来做GO KEGG; TF已经由step3转换为Gene ID
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/01.homer_motif_analysis/Important_files/"

/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/software/anaconda3/bin/findGO.pl ${wd}know_motif_gene_ID.txt human ${wd}Pathway.Results/ -human -ontology

# findGO.pl assesses the enrichment of various categories of gene function, biological pathways, domain structure, chromosome location, etc., in your gene list relative to a set of background gene IDs.  Enrichment is calculated assuming the cumulative hypergeometric distribution, much in the same way that HOMER scores motif enrichment.
#网址：http://homer.ucsd.edu/homer/microarray/go.html
