Ref:The use of two-sample methods for Mendelian randomization analyses on single large datasets 



#####ESRD mendelian randomization


###########################################Part1: Format exposure data


# Step1. Generate input files for exporsure (meQTL)

index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
Important="/share/home/liujianjun/ESRD_Project/Important_files/"
Comp="/share/home/liujianjun/ESRD_Project/Part5.Complications/Complication_Results/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
geno="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.hg38.PLINK/"
ASA_dir="/share/home/liujianjun/BGI_data/BGItrans/ESRD_ASA/PhD.paper/Impotant.files/"
QTLme="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/Common.Sample.Methylation/"
QTLme2="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.Methylation/"
ChinaMap="/share/home/liujianjun/BGI_data/BGItrans/ESRD_ASA/PhD.paper/01.Aligment.hg37tohg38/Aligment.Data.hg38.PLINK/"
obtain="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/"
Imp="/share/home/liujianjun/ESRD_Project/Important_files/"
geno_out="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Imputation.ChinaMAP/"
wd7="/share/home/liujianjun/GWAS.OF.ESRD/Raw.Data/3.PLINK.QCed.after.UpdatedID/"
var="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.hg38.PLINK/Genotype/"
split="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.hg38.PLINK/Split.Chromosome/"
class="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/meQTL.by.CpG/"


mr_ex="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Exposure_meQTL/"
mr_out="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Outcome_GWAS/"
mr_result="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Results/"
mr_code="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Code"



#Extract signifcant meQTL pairs, and add the information of position of meQTL-SNP and meQTL-CpG
library(data.table)
full=fread(paste0(obtain,"Full.cis-meQTL.Results.txt"))


CpG_loc=fread(paste0(QTLme2, "shared.DMP.hg38.Coloc")) #CpG位点的位置
 names(CpG_loc)[1]="gene"
names(CpG_loc)[2]="CpG.chr"
names(CpG_loc)[3]="CpG_pos"
CpG=CpG_loc[,1:3]

# Load the required package
library(stringr)
full$SNP_pos <- str_extract(full$SNP, "(?<=:)(\\d+)(?=:)") #SNP的位置

##Merge
anno=merge(full,CpG,by="gene")
 library(dplyr)
 anno = anno %>% rename(SNP.chr = CHR) %>% rename(P = `p-value`)
 sig=subset(anno, P < 5e-8)

##Format
library(dplyr)
sig$se=abs(sig$beta)/qnorm(1 - sig$P / 2)
sig=sig %>% select (SNP.chr,SNP_pos,SNP,A1,A2,MAF,beta,se,P,Group,gene,NCHROBS)
data = sig %>% rename(CpG = gene, CHR = SNP.chr, BP=SNP_pos)
data$N=sig$NCHROBS/2
data$NCHROBS=c()
data$CHR=paste0("chr",data$CHR)

#save
write.table(data,paste0(mr_ex,"cis-meQTL.raw.txt"),row.names=F,quote=F)



###step2 Convert SNP ID to rsID
#!usr/bin/bash

#!usr/bin/bash
mr_ex="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Exposure_meQTL/"
wd0="/share/home/liujianjun/BGI_data/BGItrans/ESRD_ASA/PhD.paper/DBSNP/"
mr_code="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Code/"


awk '
    NR==FNR {
        key = $1 " " $2
        array[key] = $0
        next
    }
    {
        key = $1 " " $2
        if (key in array) {
            print array[key], $3
        }
    }
' ${mr_ex}cis-meQTL.raw.txt  ${wd0}rsID.hg38 > ${mr_ex}cis-meQTL.rsID.Results.txt


mr_code="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Code/"
bsub -e error.log -o output.log <  ${mr_code}Convert.exposure.cis-SNP.to.rsID.sh  #SNP ID 转换为rsID的脚本

sed -i "1iCHR BP SNP A1 A2 MAF beta se P Group CpG N rsid"  ${mr_ex}cis-meQTL.rsID.Results.txt




###step3 Pruned the meQTL
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
Important="/share/home/liujianjun/ESRD_Project/Important_files/"
Comp="/share/home/liujianjun/ESRD_Project/Part5.Complications/Complication_Results/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
geno="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.hg38.PLINK/"
ASA_dir="/share/home/liujianjun/BGI_data/BGItrans/ESRD_ASA/PhD.paper/Impotant.files/"
QTLme="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/Common.Sample.Methylation/"
QTLme2="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.Methylation/"
ChinaMap="/share/home/liujianjun/BGI_data/BGItrans/ESRD_ASA/PhD.paper/01.Aligment.hg37tohg38/Aligment.Data.hg38.PLINK/"
obtain="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/"
Imp="/share/home/liujianjun/ESRD_Project/Important_files/"
geno_out="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Imputation.ChinaMAP/"
wd7="/share/home/liujianjun/GWAS.OF.ESRD/Raw.Data/3.PLINK.QCed.after.UpdatedID/"
var="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.hg38.PLINK/Genotype/"
split="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.hg38.PLINK/Split.Chromosome/"
class="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/meQTL.by.CpG/"


mr_ex="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Exposure_meQTL/"
mr_out="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Outcome_GWAS/"
mr_result="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Results/"
mr_code="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Results/Code"
panel="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Causal_Effect/Reference_Panel/"




library(data.table)
library(ieugwasr)
library(TwoSampleMR)

exp <- read_exposure_data(filename=paste0(mr_ex,"cis-meQTL.rsID.Results.txt"),sep = " ",snp_col ="rsid",beta_col = "beta", se_col = "se",effect_allele_col = "A1",other_allele_col = "A2", pval_col = "P",phenotype_col="CpG",eaf="MAF")
exp=subset(exp, pval.exposure < 5e-4)
exp=subset(exp,mr_keep.exposure=="TRUE")

a=dplyr::tibble(rsid=exp$SNP, pval=exp$pval.exposure, id=exp$exposure)


##########Pruning (HPC)

list=unique(a$id)
result=c()
for ( i in 1:length(list)) { data=subset(a,id==list[i]);data_clump=ld_clump(dat=data,clump_kb = 1000,clump_r2 = 0.2,plink_bin = "/share/home/liujianjun/software/plink",bfile =paste0(panel,"EAS"));result=rbind(result,data_clump)}

#可尝试将clump_r2 = 0.2 改为clump_r2 = 0.4


names(result)=c("SNP","pval","exposure")
LD=result[,c(1,3)]
final=merge(LD,exp,by=c("SNP","exposure"))
#save
write.table(final,paste0(mr_ex,"Pruned.cis-meQTL.exposure.txt"),row.names=FALSE,quote=FALSE)








###########################################Part2: Format outcome data

####Step1: Format the GWAS results
#input :"CHR"   "SNP"   "BP"    "A1"    "TEST"  "NMISS" "OR"    "STAT"  "P"   
#output:  "CHR" "BP"  "SNP" "A1"  "A2"  "MAF" "OR"  "b"   "se"  "P"   "N"  

#!usr/bin/Rscirpt
library(data.table)
library(dplyr)
wd1="/share/home/liujianjun/GWAS.OF.ESRD/Important.file/"


gwas=fread(paste0(wd1,"ESRD.GWAS.Results.txt"))
freq=fread(paste0(wd1,"snp_freq.frq"))
head(freq)

freq=freq %>% select(SNP,A2,MAF,NCHROBS)
head(freq)
freq$N=freq$NCHROBS/2
freq$NCHROBS=c()
head(freq)
all=merge(gwas,freq,by="SNP")

all$b=log(all$OR)  #β
all$se=abs(all$b)/qnorm(1 - all$P / 2)  #se
all$CHR=paste0("chr",all$CHR)
all = all %>% select (CHR,BP,SNP,A1,A2,MAF,OR,b,se,P,N)


write.table(all,paste0(wd1,"ESRD.GWAS.Formated.Results.txt"),row.names=F,quote=F)




#!usr/bin/bash
wd1="/share/home/liujianjun/GWAS.OF.ESRD/Important.file/"
wd0="/share/home/liujianjun/BGI_data/BGItrans/ESRD_ASA/PhD.paper/DBSNP/"

awk '
    NR==FNR {
        key = $1 " " $2
        array[key] = $0
        next
    }
    {
        key = $1 " " $2
        if (key in array) {
            print array[key], $3
        }
    }
' ${wd1}ESRD.GWAS.Formated.Results.txt  ${wd0}rsID.hg38 > ${wd1}ESRD.GWAS.rsID.Results.txt

bsub -e error.log -o output.log <  ${wd1}Convert.SNP.to.rsID.sh  #SNP ID 转换为rsID的脚本


mr_ex="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Exposure_meQTL/"
mr_out="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Outcome_GWAS/"
mr_result="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Results/"
mr_code="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Code/"
wd1="/share/home/liujianjun/GWAS.OF.ESRD/Important.file/"

cp ${wd1}ESRD.GWAS.rsID.Results.txt  ${mr_out}ESRD.GWAS.rsID.Results.txt
sed -i '1iCHR BP SNP A1 A2 MAF OR b se P N rsid' ${mr_out}ESRD.GWAS.rsID.Results.txt
cp  ${wd1}Convert.SNP.to.rsID.sh   ${mr_code}Convert.Outcome.SNP.to.rsID.sh 



####################################Harmonize data :ESRD
library(data.table)
library(ieugwasr)
library(TwoSampleMR)

out<- read_outcome_data(filename=paste0(mr_out,"ESRD.GWAS.rsID.Results.txt"),snps=NULL,sep = " ",snp_col ="rsid",beta_col = "b", se_col = "se",effect_allele_col = "A1",other_allele_col = "A2",eaf_col = "MAF", pval_col = "P")
out$outcome="ESRD"

final=fread(paste0(mr_ex,"Pruned.cis-meQTL.exposure.txt"))


dat <- harmonise_data(exposure_dat = final,outcome_dat = out)
df=as.data.frame(unlist(table(dat$exposure)))
df$Method <- ifelse(df$Freq >= 2, "ESRD_and_eGFR", "only_eGFR")
names(df)[1]="CpG"

write.table(dat,paste0(mr_out,"cis-Harmonised.Data.txt"),row.names=FALSE,quote=FALSE)
write.table(df,paste0(mr_result,"cis-ESRD_Instruments_counts.txt"),row.names=FALSE,quote=FALSE)







##########################MR analysis
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"


hm=fread(paste0(mr_out,"cis-Harmonised.Data.txt"))
hm=subset(hm,mr_keep.outcome == "TRUE")

#F-statistics
dat_ESRD=hm
dat_ESRD$samplesize.exposure=823
dat_ESRD$samplesize.outcome <- as.numeric(dat_ESRD$samplesize.outcome)
dat_ESRD$samplesize.outcome=3882
dat_ESRD=steiger_filtering(dat_ESRD)
dat_ESRD=subset(dat_ESRD,steiger_dir== "TRUE")

dat=dat_ESRD

list=unique(dat$exposure)
F=c()
length(list)

for ( i in 1:length(list)) {
    data=subset(dat,exposure==list[i]);n=1848;k=nrow(data); 
R2=data$rsq.exposure;left=(n-k-1)/k; 
right=R2/(1-R2);
data$F=left*right; 
data$Fmean=(sum(data$F))/k;
F=rbind(F,data)
}

hm=subset(F,F >=10)

##ivw
ESRD_ivw <- mr(hm , method_list = "mr_ivw")


##Wald ratio
dat_ESRD=hm
freq_table <- table(dat_ESRD$exposure) 
data<- data.frame(
  CpG = names(freq_table),
  Freq = as.numeric(freq_table)
)

sub=subset(data,Freq ==1)  #工具变量只有1个的CpG位点

list=unique(dat_ESRD$exposure)  ##Wald ratio
file=c()
for ( i in 1:length(list)) { data=subset(dat_ESRD,exposure ==list[i]);new=data[which(data$pval.exposure==min(data$pval.exposure)),];file=as.data.frame(rbind(file,new))}
file=subset(file, exposure %in% sub$CpG)
ESRD_wald<- mr(file , method_list = "mr_wald_ratio")


sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")
names(sentinels)[1]="CpG"
assoc=sentinels[,c("CpG","Bacon_b_DM","UCSC_RefGene_Name","UCSC_RefGene_Group","CHR","MAPINFO")]
names(ESRD_ivw)[4]="CpG"
ivw_result=merge(ESRD_ivw,assoc,by="CpG",All.x=T)
ivw_result$dir=ifelse(ivw_result$Bacon_b_DM * ivw_result$b >0,"yes","no")
subset(ivw_result,dir == "yes" & pval < 0.05)


sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")
names(sentinels)[1]="CpG"
assoc=sentinels[,c("CpG","Bacon_b_DM","UCSC_RefGene_Name","UCSC_RefGene_Group","CHR","MAPINFO")]
names(ESRD_wald)[4]="CpG"
wald_result=merge(ESRD_wald,assoc,by="CpG",All.x=T)
wald_result$dir=ifelse(wald_result$Bacon_b_DM * wald_result$b >0,"yes","no")
subset(wald_result,dir == "yes"& pval < 0.05)

all_ESRD=rbind(ivw_result,wald_result)



write.table(ivw_result,paste0(mr_result,"ESRD_ivw_MR_Result.txt"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(wald_result,paste0(mr_result,"ESRD_wald_MR_Result.txt"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(all_ESRD,paste0(mr_result,"ESRD_All_MR_Result.txt"),quote=FALSE,row.names=FALSE,sep="\t")


########################################Harmonized data: eGFR GWAS summary were downloaded from https://www.nature.com/articles/s41467-024-53516-7#data-availability
library(ieugwasr)
library(TwoSampleMR)
 library(data.table)

path_1="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Causal_Effect/Published.GWAS.Datasets/"

final=fread(paste0(mr_ex,"Pruned.cis-meQTL.exposure.txt")) #exposure

#eGFR Asian
out_eGFR<- read_outcome_data(filename=paste0(path_1,"TWB2_BBJ_eGFR_hg19_METAL_FUMA_noNA"),snps=NULL,sep = " ",snp_col ="SNP",beta_col = "Effect", se_col = "standard_error",effect_allele_col = "Allele1",other_allele_col = "Allele2", pval_col = "P")
out_eGFR$outcome="eGFR"

dat_eGFR <- harmonise_data(exposure_dat = final,outcome_dat = out_eGFR)
 dat_eGFR=subset(dat_eGFR,mr_keep == "TRUE")


#eGFR Gabriel 

out_eGFR1<- read_outcome_data(filename=paste0(path_1,"UCSF.Gabriel.eGFR.GWAS.tsv"),snps=NULL,sep = "\t",snp_col ="rsid",beta_col="beta", se_col = "standard_error",effect_allele_col = "effect_allele",other_allele_col = "other_allele",eaf_col="effect_allele_frequency", pval_col = "p_value")
out_eGFR1$outcome="eGFR1"

dat_eGFR1 <- harmonise_data(exposure_dat = final,outcome_dat = out_eGFR1)
 dat_eGFR1=subset(dat_eGFR1,mr_keep == "TRUE")

 #F-statistics  （https://www.sciencedirect.com/science/article/pii/S1389945723004148； Selection of genetic instruments in Mendelian randomisation studies of sleep traits）


###该文献报导有两种方法：1) (F=(n−k−1/k) (R2/1−R2))  2)F = β2/SE2;    我们选择方法(F = β2/SE2）

 #When the R2 is unknown, the ‘t-statistic’ summary-level method can be used (F = β2/SE2). In this case, the F-statistic will be an approximation because it uses the sample size for the discovery GWA study, not the one from the data under analysis. Finally, the “MendelianRandomization” R package allows the calculation of the F-statistic [13].


dat_eGFR1$samplesize.exposure=923
dat_eGFR1$samplesize.outcome=406504
dat_eGFR1=steiger_filtering(dat_eGFR1)
dat_eGFR1=subset(dat_eGFR1,steiger_dir== "TRUE")

dat=dat_eGFR1

list=unique(dat$exposure)
length(list)

Fdata=c()
for ( i in 1:length(list)) {
data=subset(dat,exposure==list[i]);
k=nrow(data);
beta2=(data$beta.exposure)*(data$beta.exposure);
se2=(data$se.exposure)*(data$se.exposure);
Fvalue=beta2/se2
data$Fvalue=Fvalue; 
data$Fmean=(sum(data$Fvalue))/k;
Fdata=rbind(Fdata,data)
}

dat_eGFR1=subset(Fdata, F >=10)


write.table(dat_eGFR1,paste0(table,"Table.eGFR.Instruments.txt",sep="\t",row.names=F,quote=F)



#BUN
out_BUN<- read_outcome_data(filename=paste0(path_1,"BUN.mean.GWAS.tsv"),samplesize_col="n",snps=NULL,sep = "\t",snp_col ="rsid",beta_col = "beta", se_col = "standard_error",effect_allele_col = "effect_allele",other_allele_col = "other_allele", pval_col = "p_value",eaf_col="effect_allele_frequency")
out_BUN$outcome="BUN"

dat_BUN<- harmonise_data(exposure_dat = final,outcome_dat = out_BUN)
dat_BUN=subset(dat_BUN,mr_keep == "TRUE")




#F-statistics
dat_BUN=steiger_filtering(dat_BUN)
dat_BUN=subset(dat_BUN,steiger_dir== "TRUE")

dat=dat_BUN

list=unique(dat$exposure)
length(list)

Fdata=c()
for ( i in 1:length(list)) {
data=subset(dat,exposure==list[i]);
k=nrow(data);
beta2=(data$beta.exposure)*(data$beta.exposure);
se2=(data$se.exposure)*(data$se.exposure);
Fvalue=beta2/se2
data$Fvalue=Fvalue; 
data$Fmean=(sum(data$Fvalue))/k;
Fdata=rbind(Fdata,data)
}


dat_BUN=subset(Fdata, Fvalue >=10)


##########################MR analysis
#eGFR
frame <- mr(dat_eGFR , method_list = "mr_ivw")  #IVW

list=unique(dat_eGFR$exposure)  ##Wald ratio
file=c()
for ( i in 1:length(list)) { data=subset(dat_eGFR,exposure ==list[i]);new=data[which(data$pval.exposure==min(data$pval.exposure)),];file=as.data.frame(rbind(file,new))}
wald <- mr(file , method_list = "mr_wald_ratio")


#eGFR1

eGFR1_ivw <- mr(dat_eGFR1 , method_list = "mr_ivw")  #IVW



freq_table <- table(dat_eGFR1$exposure)
data<- data.frame(
  CpG = names(freq_table),
  Freq = as.numeric(freq_table)
)

sub=subset(data,Freq ==1)  #工具变量只有1个的CpG位点


list=unique(dat_eGFR1$exposure)  ##Wald ratio
file=c()
for ( i in 1:length(list)) { data=subset(dat_eGFR1,exposure ==list[i]);new=data[which(data$pval.exposure==min(data$pval.exposure)),];file=as.data.frame(rbind(file,new))}
file=subset(file, exposure %in% sub$CpG)
eGFR1_wald<- mr(file , method_list = "mr_wald_ratio")



sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")
names(sentinels)[1]="CpG"
assoc=sentinels[,c("CpG","Bacon_b_DM","UCSC_RefGene_Name","UCSC_RefGene_Group","CHR","MAPINFO")]
names(eGFR1_ivw)[4]="CpG"
ivw_result=merge(eGFR1_ivw,assoc,by="CpG",All.x=T)
ivw_result$dir=ifelse(ivw_result$Bacon_b_DM * ivw_result$b <0,"yes","no")
subset(ivw_result,dir == "yes" & pval < 0.05)

names(eGFR1_wald)[4]="CpG"
wald_result=merge(eGFR1_wald,assoc,by="CpG",All.x=T)
wald_result$dir=ifelse(wald_result$Bacon_b_DM * wald_result$b <0,"yes","no")
subset(wald_result,dir == "yes"& pval < 0.05)

all_eGFR1=rbind(ivw_result,wald_result)

write.table(ivw_result,paste0(mr_result,"eGFR1_ivw_MR_Result.txt"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(wald_result,paste0(mr_result,"eGFR1_wald_MR_Result.txt"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(all_eGFR1,paste0(mr_result,"eGFR1_All_MR_Result.txt"),quote=FALSE,row.names=FALSE,sep="\t")





#BUN
BUN_ivw <- mr(dat_BUN , method_list = "mr_ivw")  #IVW

freq_table <- table(dat_BUN$exposure)
data<- data.frame(
  CpG = names(freq_table),
  Freq = as.numeric(freq_table)
)

sub=subset(data,Freq ==1)  #工具变量只有1个的CpG位点


list=unique(dat_BUN$exposure)  ##Wald ratio
file=c()
for ( i in 1:length(list)) { data=subset(dat_BUN,exposure ==list[i]);new=data[which(data$pval.exposure==min(data$pval.exposure)),];file=as.data.frame(rbind(file,new))}
file=subset(file, exposure %in% sub$CpG)
BUN_wald<- mr(file , method_list = "mr_wald_ratio")

sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")
names(sentinels)[1]="CpG"
assoc=sentinels[,c("CpG","Bacon_b_DM","UCSC_RefGene_Name","UCSC_RefGene_Group","CHR","MAPINFO")]
names(BUN_ivw)[4]="CpG"
ivw_result=merge(BUN_ivw,assoc,by="CpG",All.x=T)
ivw_result$dir=ifelse(ivw_result$Bacon_b_DM * ivw_result$b >0,"yes","no")
subset(ivw_result,dir == "yes" & pval < 0.05)

names(BUN_wald)[4]="CpG"
wald_result=merge(BUN_wald,assoc,by="CpG",All.x=T)
wald_result$dir=ifelse(wald_result$Bacon_b_DM * wald_result$b >0,"yes","no")
subset(wald_result,dir == "yes"& pval < 0.05)

all_BUN=rbind(ivw_result,wald_result)

write.table(ivw_result,paste0(mr_result,"BUN_ivw_MR_Result.txt"),quote=F,row.names=F,sep="\t")
write.table(wald_result,paste0(mr_result,"BUN_wald_MR_Result.txt"),quote=F,row.names=F,sep="\t")
write.table(all_BUN,paste0(mr_result,"BUN_All_MR_Result.txt"),quote=F,row.names=F,sep="\t")









##合并MR 位点与andidate gene


path1="/share/home/liujianjun/ESRD_Project/Table/Input/"
path2="/share/home/liujianjun/ESRD_Project/Table/Output/"
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
mr_result="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Results/"

library(data.table)
library(dplyr)

setwd(mr_result)
#MR results
eGFR=fread("eGFR1_All_MR_Result.txt")



#result=subset(disease, pval < 0.05)
result=eGFR
result$id.outcome=c()
result$id.exposure=c()

#Prioritization gene
setwd(path1)
gene=fread("Final.ESRD_DMP_Gene.txt")


data1=gene
data1$eQTM_blood <- ifelse(rowSums(data1[, c("eqtm_2015", "eqtm_2018_GTP", "eqtm_2018_MESA", "eqtm_2023")]) >= 1, 1, 0)
data1$eQTM_kidney <- ifelse(rowSums(data1[, c("eqtm_2022_kidney")]) >= 1, 1, 0)
data1$Enhancer_linked_gene <- ifelse(rowSums(data1[, c("atac", "epimap", "hic","ABC")]) >= 1, 1, 0)

new_dt1=subset(data1, select = -c(atac, epic_nearest, epimap, hic, eqtm_2015, eqtm_2018_GTP, eqtm_2018_MESA, eqtm_2022_kidney, eqtm_2023, ABC, group, sum))

###选取eQTM 位点
new_dt2=subset(new_dt1,eQTM_blood > 0 )

new_dt3=setDT(new_dt2)[, .(gene = paste(gene, collapse = ";")), by = CpG]
names(new_dt3)[2]="eQTM_gene"

###合并
result1=merge(result,new_dt3,by="CpG",all.x=T)
 result1=result1[,c("CpG","eQTM_gene")]
write.table(result1,paste0(mr_result,"eGFR_MR_result_eQTMgene.txt"),row.names=F,quote=F)

###选取enhancer 位点
new_dt2=subset(new_dt1,Enhancer_linked_gene > 0 )


new_dt3=setDT(new_dt2)[, .(gene = paste(gene, collapse = ";")), by = CpG]
names(new_dt3)[2]="Enhancer_linked_gene"

##合并
result1=merge(result,new_dt3,by="CpG",all.x=T)
 result1=result1[,c("CpG","Enhancer_linked_gene")]
write.table(result1,paste0(mr_result,"eGFR_MR_result_Enhancer_linked_gene.txt"),row.names=F,quote=F)




####CAD Complivations: stroke MI HF

#exposure data

library(ieugwasr)
library(TwoSampleMR)
library(data.table)

path_1="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Causal_Effect/Published.GWAS.Datasets/"

final=fread(paste0(mr_ex,"Pruned.cis-meQTL.exposure.txt")) #exposure



##Outcome data

comp="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/ESRD_ASA/PartD2.ESRD_MR_causal_analysis/ESRD.Complications/"

#stroke
stroke<- read_outcome_data(filename=paste0(comp,"Stroke.GWAS.tsv"),sep = "\t",snps=NULL,snp_col ="rsid",beta_col = "beta", se_col = "standard_error",effect_allele_col = "effect_allele",other_allele_col = "other_allele", pval_col = "p_value",eaf_col="effect_allele_frequency")
stroke$outcome = "Stroke"

##harmonised data
dat_stroke <- harmonise_data(exposure_dat = final,outcome_dat = stroke)


###MR 
stroke_ivw <- mr(dat_stroke , method_list = "mr_ivw")  #IVW



freq_table <- table(dat_stroke$exposure)  #Wald ratio
data<- data.frame(
  CpG = names(freq_table),
  Freq = as.numeric(freq_table)
)

sub=subset(data,Freq ==1)  #工具变量只有1个的CpG位点


list=unique(dat_stroke$exposure)  ##Wald ratio
file=c()
for ( i in 1:length(list)) { data=subset(dat_stroke,exposure ==list[i]);new=data[which(data$pval.exposure==min(data$pval.exposure)),];file=as.data.frame(rbind(file,new))}
file=subset(file, exposure %in% sub$CpG)
stroke_wald<- mr(file , method_list = "mr_wald_ratio")



sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")
names(sentinels)[1]="CpG"
assoc=sentinels[,c("CpG","Bacon_b_DM","UCSC_RefGene_Name","UCSC_RefGene_Group","CHR","MAPINFO")]
names(stroke_ivw)[4]="CpG"
ivw_result=merge(stroke_ivw,assoc,by="CpG",All.x=T)
ivw_result$dir=ifelse(ivw_result$Bacon_b_DM * ivw_result$b >0,"yes","no")
subset(ivw_result,dir == "yes" & pval < 0.05)

names(stroke_wald)[4]="CpG"
wald_result=merge(stroke_wald,assoc,by="CpG",All.x=T)
wald_result$dir=ifelse(wald_result$Bacon_b_DM * wald_result$b >0,"yes","no")
subset(wald_result,dir == "yes"& pval < 0.05)

all_stroke=rbind(ivw_result,wald_result)

write.table(ivw_result,paste0(mr_result,"stroke_ivw_MR_Result.txt"),quote=F,row.names=F,sep="\t")
write.table(wald_result,paste0(mr_result,"stroke_wald_MR_Result.txt"),quote=F,row.names=F,sep="\t")
write.table(all_stroke,paste0(mr_result,"stroke_All_MR_Result.txt"),quote=F,row.names=F,sep="\t")




#HF
HF<- read_outcome_data(filename=paste0(comp,"HF.GWAS.tsv"),sep = "\t",snps=NULL,snp_col ="rsid",beta_col = "beta", se_col = "standard_error",effect_allele_col = "effect_allele",other_allele_col = "other_allele", pval_col = "p_value",eaf_col="effect_allele_frequency")
HF$outcome = "HF"



##harmonised data
dat_HF <- harmonise_data(exposure_dat = final,outcome_dat = HF)


###MR 
HF_ivw <- mr(dat_HF , method_list = "mr_ivw")  #IVW

freq_table <- table(dat_HF$exposure) ##Wald ratio
data<- data.frame(
  CpG = names(freq_table),
  Freq = as.numeric(freq_table)
)

sub=subset(data,Freq ==1)  #工具变量只有1个的CpG位点


list=unique(dat_HF$exposure)  ##Wald ratio
file=c()
for ( i in 1:length(list)) { data=subset(dat_HF,exposure ==list[i]);new=data[which(data$pval.exposure==min(data$pval.exposure)),];file=as.data.frame(rbind(file,new))}
file=subset(file, exposure %in% sub$CpG)
HF_wald<- mr(file , method_list = "mr_wald_ratio")


sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")
names(sentinels)[1]="CpG"
assoc=sentinels[,c("CpG","Bacon_b_DM","UCSC_RefGene_Name","UCSC_RefGene_Group","CHR","MAPINFO")]
names(HF_ivw)[4]="CpG"
ivw_result=merge(HF_ivw,assoc,by="CpG",All.x=T)
ivw_result$dir=ifelse(ivw_result$Bacon_b_DM * ivw_result$b >0,"yes","no")
subset(ivw_result,dir == "yes" & pval < 0.05)

names(HF_wald)[4]="CpG"
wald_result=merge(HF_wald,assoc,by="CpG",All.x=T)
wald_result$dir=ifelse(wald_result$Bacon_b_DM * wald_result$b >0,"yes","no")
subset(wald_result,dir == "yes"& pval < 0.05)

all_HF=rbind(ivw_result,wald_result)


write.table(ivw_result,paste0(mr_result,"HF_ivw_MR_Result.txt"),quote=F,row.names=F,sep="\t")
write.table(wald_result,paste0(mr_result,"HF_wald_MR_Result.txt"),quote=F,row.names=F,sep="\t")
write.table(all_HF,paste0(mr_result,"HF_All_MR_Result.txt"),quote=F,row.names=F,sep="\t")







#MI
MI<- read_outcome_data(filename=paste0(comp,"MI.GWAS.tsv"),sep = "\t",snps=NULL,snp_col ="rsid",beta_col = "beta", se_col = "standard_error",effect_allele_col = "effect_allele",other_allele_col = "other_allele", pval_col = "p_value",eaf_col="effect_allele_frequency")
MI$outcome = "MI"



##harmonised data
dat_MI <- harmonise_data(exposure_dat = final,outcome_dat = MI)


###MR 
MI_ivw <- mr(dat_MI , method_list = "mr_ivw")  #IVW

freq_table <- table(dat_MI$exposure) ##Wald ratio
data<- data.frame(
  CpG = names(freq_table),
  Freq = as.numeric(freq_table)
)

sub=subset(data,Freq ==1)  #工具变量只有1个的CpG位点


list=unique(dat_MI$exposure)  ##Wald ratio
file=c()
for ( i in 1:length(list)) { data=subset(dat_MI,exposure ==list[i]);new=data[which(data$pval.exposure==min(data$pval.exposure)),];file=as.data.frame(rbind(file,new))}
file=subset(file, exposure %in% sub$CpG)
MI_wald<- mr(file , method_list = "mr_wald_ratio")



sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")
names(sentinels)[1]="CpG"
assoc=sentinels[,c("CpG","Bacon_b_DM","UCSC_RefGene_Name","UCSC_RefGene_Group","CHR","MAPINFO")]
names(MI_ivw)[4]="CpG"
ivw_result=merge(MI_ivw,assoc,by="CpG",All.x=T)
ivw_result$dir=ifelse(ivw_result$Bacon_b_DM * ivw_result$b >0,"yes","no")
subset(ivw_result,dir == "yes" & pval < 0.05)

names(MI_wald)[4]="CpG"
wald_result=merge(MI_wald,assoc,by="CpG",All.x=T)
wald_result$dir=ifelse(wald_result$Bacon_b_DM * wald_result$b >0,"yes","no")
subset(wald_result,dir == "yes"& pval < 0.05)


all_MI=rbind(ivw_result,wald_result)


write.table(ivw_result,paste0(mr_result,"MI_ivw_MR_Result.txt"),quote=F,row.names=F,sep="\t")
write.table(wald_result,paste0(mr_result,"MI_wald_MR_Result.txt"),quote=F,row.names=F,sep="\t")
write.table(all_MI,paste0(mr_result,"MI_All_MR_Result.txt"),quote=F,row.names=F,sep="\t")

