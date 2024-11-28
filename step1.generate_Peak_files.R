##########################This script was aimed to generated peak files to do the motif analysis.
anno=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/01.homer_motif_analysis/Important_files/EPIC.Annotation.BED.files.txt",header=T)
DMP=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/01.homer_motif_analysis/Important_files/Bacon.Annotate.DMP.txt",header=T)

anno_sub=subset(anno,CpG %in% DMP$CpG)
anno_sub=anno_sub[,c(1,3,4,5,2)]
 names(anno_sub)=c("PeakID","Chr","Start","End","Strand")
anno_sub$Strand[which(anno_sub$Strand=="F")] <- "+"
anno_sub$Strand[which(anno_sub$Strand=="R")] <- "-"


write.table(anno_sub,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/01.homer_motif_analysis/Important_files/DMPpeaks.txt",row.names=F,sep="\t")
