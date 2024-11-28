group="adipose	blood_t_cell	bone	brain	digestive	endocrine	endothelial	epithelial	esc	es	eye	heart	hsc_b_cell	ipsc	kidney	liver	lung	lymphoblastoid	mesench	muscle	myosat	neurosph	pancreas	placenta	pns	reproductive	smmuscle	spleen	stromal	thymus	urinary";

for i in ${group} ; 
do
echo " #!/usr/bin/env Rscript
library(data.table)
setwd(\"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Important.files/\")

sentinels=read.table(\"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.Annotate.DMP.txt\",header=T)
bed=fread(\"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartH.Target_gene_annotation/epimap_link/Epimap.by.Pergroup/links_by_group.$i.tsv\")
all= read.table(\"All.Anno.CpG.Mean.SD.txt\",header=T)
load(\"matches.Rdata\")
matches=as.data.frame(matches)

#############################Background CpG sites
#result=c()
#file=c()
#length=c()

#	for (j in 1:ncol(matches)) { 
#	result=c(); 
#	print (paste0(\"Background.CpG.set\",j));
#	anno=subset(all, CpG %in% matches[,j] )[,1:3];  ########第一个循环是读取每一个background CpG set, 并注释位置信息;
#		for ( h in 1:nrow(sentinels)) {
#		CpG.open=subset(bed, chr %in% anno\$CHR[h] & start < anno\$BP[h] & end > anno\$BP[h]);if ( nrow(CpG.open) > 1) {result=as.data.frame(rbind(result,CpG.open[1,]))} else {result=as.data.frame(rbind(result,CpG.open))}};		
#	length=c(length,nrow(result))
#	};  #############第二个循环是对每一个background CpG set中的nrow(sentinels)个CpG进行 open chromatin的注释，并得到与open chromatin overlap的数量;
#
#       file=as.data.frame(t(length))
#	sum=sum(length)
#	colnames(file)=paste0(\"Open.Num\",1:ncol(file)) 
#	file\$Probability=sum/(nrow(sentinels)*1000)
#       rownames(file)=\"$i\"
#write.table(file,\"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/Background.${i}.Number.txt\")

#############################Target CpG sites

anno=sentinels[,1:3]
num=c()
result=c()
	for ( h in 1:nrow(sentinels)) {
                CpG.open=subset(bed, chr %in% anno\$Chromosome[h] & start < anno\$Start[h] & end > anno\$Start[h]);if ( nrow(CpG.open) > 1) {result=as.data.frame(rbind(result,CpG.open[1,]))} else {result=as.data.frame(rbind(result,CpG.open))}
	}
num=as.data.frame(nrow(result))
colnames(num)[1]=\"Open.num\"
num\$Open.num=as.numeric(num\$Open.num)
num\$Close.num=nrow(sentinels)- num\$Open.num
num\$Sample.ID=\"$i\"
num=num[,c(3,1,2)]  ####Sample.ID Open.num Close.num

write.table(num,\"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/DMP.${i}.Number.txt\")

num=read.table(\"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/DMP.${i}.Number.txt\",header=T)
file=read.table(\"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/Background.${i}.Number.txt\",header=T)

file=as.data.frame(file)
binomial <- binom.test(num\$Open.num, n=nrow(sentinels), file\$Probability)
frame=data.frame(Sample.ID=\"$i\",Obv.Freq=num\$Open.num/nrow(sentinels),Exp.Freq=file\$Probability,FC=(num\$Open.num/nrow(sentinels))/file\$Probability,Pval=binomial\$p.value)


write.table(frame, \"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/Final.Enrichment.Pval.${i}.txt\",row.names=F)  " > /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/permutation.Enh.${i}.Bino.R;
done
