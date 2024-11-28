setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Important.files")
info=read.csv("Blood.Smaller.P.tissue.txt",header=T,sep="\t")

info=info[,1:4]

frame=read.table("/hwfssz5/ST_HEALTH/P20Z10200N0041/zhouxiaohong/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/All.Enrichment.Fisher.Pval.txt",header=T)
file=merge(frame,info,by="Sample.ID")

write.table(file,"/hwfssz5/ST_HEALTH/P20Z10200N0041/zhouxiaohong/PartG.Functional_analysis/04.chromatin_accessibility/Enrichment.Fisher.Pval.result.txt",row.names=F)
