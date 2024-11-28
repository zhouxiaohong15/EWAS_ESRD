######################This cript was aimed to transfer TFs format to Gene ID, facilitating to do the KEGG and GO analysis. 主要是对knownResults数据库中的motif结合的TF进行转换。

setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/01.homer_motif_analysis/Important_files/")

result=read.csv("knownResults.txt",header=T,sep="\t")

######提取显著的TF的名称
result$Motif <- str_extract(result$Motif.Name, ".+(?=\\()")
 result$TF <- str_extract(result$Motif, ".+(?=\\()")
names(result)[5]="Pvalue.adj"
sig=subset(result,Pvalue.adj < 0.05)

############将显著的TF的GENE SYMBOL转换为Entrez Gene IDs, 网站https://metascape.org/gp/index.html#/main/step1上进行转换
gene=sig$TF
sig[,10:11]
########转换后的结果保存
