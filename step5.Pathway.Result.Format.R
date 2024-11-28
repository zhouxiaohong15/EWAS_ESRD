setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/01.homer_motif_analysis/Important_files/Pathway.Results")

kegg=read.csv("kegg.txt",header=T,sep="\t")

library(dplyr)
 kegg_unique <- kegg %>% distinct(Term, .keep_all = TRUE)
 kegg_unique$P.adj=p.adjust(kegg_unique$Enrichment,method = "fdr")

kegg_unique[1:10,c(1:3,12)]

write.table(kegg_unique,"kegg.final.txt",row.names=F)
