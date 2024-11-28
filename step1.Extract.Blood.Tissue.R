#####################本脚本是基于ATAC数据，分析591个位点在血液组织中的染色质开放区域的富集程度

###########################本富集分析需要两个文件：1）所用的染色质开放区域的ATAC bed文件下载于https://bio.liclab.net/ATACdb/；文件名：Accessible_chromatin_region_all.bed2） 各个组织和样本的表型文件： https://bio.liclab.net/ATACdb/Download.php，若要获取血液组织，则在该页面输入blood，然后选择显示100行，将显示结果粘贴到execl表格中，最终得到500多个sample对应的blood组织和细胞的表型文件：Blood.tissue.txt

setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Important.files")


##########提取所有血液组织和细胞的ATAC bed文件
library(data.table)
bed=fread("/hwfssz1/pub/database/bio.liclab.net/ATACdb/download/packages/Accessible_chromatin_region_all.bed")

blood=read.csv("Blood.tissue.txt",header=T,sep="\t")
head(blood)

blood$ID=gsub("Sample_", "", blood$Sample.ID)

data=subset(bed, sample_ID %in% blood$ID)
head(data)

write.table(data,"Blood.Tissue.ATAC.bed",row.names=F)

######共508个sample对应的bed文件
