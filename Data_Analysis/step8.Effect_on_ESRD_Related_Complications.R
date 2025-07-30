

#############################################Complications


################Extract Shared DMP.Residuals:
#!usr/bash/bin

dir="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"

awk -F ' ' '{print $1 }' ${index}Full.Rep.DMP.Results  > ${index}CpG.list.txt
grep -F -f ${index}CpG.list.txt ${dir}GDPH_998Samples.residual_mvals.txt  > ${index}Shared.DMP.Residuals.txt
head -1 ${dir}GDPH_998Samples.residual_mvals.txt > ${dir}header.txt
cp ${dir}header.txt ${index}header.txt
cat  ${index}Shared.DMP.Residuals.txt >> ${index}header.txt
mv ${index}header.txt ${index}Shared.DMP.Residuals.txt

 


##################################################Analysis

#!usr/bin/Rscript
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
Important="/share/home/liujianjun/ESRD_Project/Important_files/"
Comp="/share/home/liujianjun/ESRD_Project/Part5.Complications/Complication_Results/"

library(data.table)
data=read.csv(paste0(index,"Full.Rep.DMP.Results"),sep=" ")
head(data)
res=fread(paste0(index,"Shared.DMP.Residuals.txt"))
res=as.data.frame(res)
rownames(res)=res$V1
res$V1=c()

residual=as.data.frame(t(res))
residual$Sample_ID=rownames(residual)

################Pheno file
pheno=read.csv(paste0(Important,"GDPH.1166.Samples.Pheno.txt"),sep="\t")
names(pheno)[36]="PTH"
names(pheno)[32]="Ca"
names(pheno)[33]="P"
names(pheno)[35]="Hb"
names(pheno)[45]="CVD"
cov=pheno[,c(2:19,32,33,36,35,45)]
names(cov)
cov$Smoking=cov$PredictedSmokingStatus
cov$Smoking[which(cov$PredictedSmokingStatus %in% c("Never Smoker","Former Smoker"))] <-0
cov$Smoking[which(cov$PredictedSmokingStatus %in% c("Current Smoker"))] <- 1



#############################################################severe hyperparathyroidism (PTH concentrations ≥ 600 pg/mL) 根据KDIGO（肾脏病改善全球预后）指南，建议将ESRD患者的PTH水平维持在正常上限的2-9倍。以常见的实验室正常上限65 pg/mL为例，这一范围大约是130-585 pg/mL。
#####https://pmc.ncbi.nlm.nih.gov/articles/PMC4991918/Moderate to severe hyperparathyroidism (PTH concentrations ≥ 600 pg/mL) may increase risk of cardiovascular death,7 though the causality of this association is debatable.
###############Pheno file + methylation value
cov$PTH_Group=cov$PTH
cov$PTH_Group[which((cov$PTH > 600) )] <- 1
cov$PTH_Group[which(cov$PTH < 600 & cov$PTH > 150)] <- 0
cov$PTH_Group[which(cov$PTH < 150)] <- "NA"

cov$PTH_Group=as.numeric(cov$PTH_Group)

cov$PTH_Group=as.numeric(cov$PTH_Group)

all=c()
all=merge(cov,residual,by="Sample_ID")
all=subset(all,Sample_Group == 1)  #只在ESRD患者中进行分析
table(all$PTH_Group)

result=c()
for ( i in 26:ncol(all)) { fit = glm(all$PTH_Group~ all[,i] + all$Sample_Group+all$Age+all$Gender+all$Smoking,family=binomial); result=as.data.frame(rbind(result,c(colnames(all)[i],coef(summary(fit))[2,c(1,2,4)])));names(result)=c("CpG","Estimate" ,"SE" ,"Pvalue")}


head(result)
result$Pvalue
str(result)
result$Pvalue=as.numeric(result$Pvalue)
result$SE=as.numeric(result$SE)
result$Estimate=as.numeric(result$Estimate)
result$Trait="PTH"

result=result[order(result$Pvalue),]
head(result)

###############PTH reuslt + ESRD result
names(result)[2:5]=paste0("PTH_",colnames(result)[2:5])
final= merge(result,data,by="CpG")
summary(glm(final$PTH_Estimate ~ final$Bacon_b_meta,family=gaussian))

final_PTH=subset(final,PTH_Pvalue < 0.05/52)
names(final_PTH)[2:5]=c("Estimate","SE","Pvalue","Trait")

cor(final$PTH_Estimate, final$Bacon_b_meta, method = "spearman")

write.table(final,paste0(Comp,"PTH.Results.txt"),row.names=F,quote=F)
###slope: 0.20153  P: 2.97e-10




#####################Ca levels (low: 0-2.1 high: >2.1) : KDIGO（2017）慢性肾脏病-矿物质与骨异常（CKD-MBD）指南，针对终末期肾病（ESRD）患者，**低钙血症（hypocalcemia）**的定义为< 2.1 mmol/L

pheno=read.csv(paste0(Important,"GDPH.1166.Samples.Pheno.txt"),sep="\t")
names(pheno)[36]="PTH"
names(pheno)[32]="Ca"
names(pheno)[33]="P"
names(pheno)[35]="Hb"
names(pheno)[45]="CVD"
cov=pheno[,c(2:19,32,33,36,35,45)]
names(cov)
cov$Smoking=cov$PredictedSmokingStatus
cov$Smoking[which(cov$PredictedSmokingStatus %in% c("Never Smoker","Former Smoker"))] <-0
cov$Smoking[which(cov$PredictedSmokingStatus %in% c("Current Smoker"))] <- 1

cov$Ca_Group[which(cov$Ca < 2.1 )] <- 1
cov$Ca_Group[which(cov$Ca > 2.1 )] <- 0


cov$Ca_Group=as.numeric(cov$Ca_Group)
table(cov$Ca_Group)
all=merge(cov,residual,by="Sample_ID")
all=subset(all,Sample_Group == 1)  #只在ESRD患者中进行分析
table(all$Ca_Group)


result=c()
for ( i in 26:ncol(all)) { fit = glm(all$Ca_Group~all[,i] + all$Sample_Group+all$Age+all$Gender+all$Smoking,family=binomial); result=as.data.frame(rbind(result,c(colnames(all)[i],coef(summary(fit))[2,c(1,2,4)])));names(result)=c("CpG","Estimate" ,"SE" ,"Pvalue")}

head(result)
result$Pvalue
str(result)
result$Pvalue=as.numeric(result$Pvalue)
result$SE=as.numeric(result$SE)
result$Estimate=as.numeric(result$Estimate)
result$Trait="Ca"

result=result[order(result$Pvalue),]
head(result)


###############Ca reuslt + ESRD result
names(result)[2:5]=paste0("Ca_",colnames(result)[2:5])
final= merge(result,data,by="CpG")
summary(glm(final$Ca_Estimate ~ final$Bacon_b_meta,family=gaussian))

final_Ca=subset(final,Ca_Pvalue < 0.05/52)
names(final_Ca)[2:5]=c("Estimate","SE","Pvalue","Trait")

 cor(final$Ca_Estimate, final$Bacon_b_meta, method = "spearman")

write.table(final,paste0(Comp,"Ca.Results.txt"),row.names=F,quote=F)


#####################P levels (low: 0-1.45 high: >1.45) :mmol/L

pheno=read.csv(paste0(Important,"GDPH.1166.Samples.Pheno.txt"),sep="\t")
names(pheno)[36]="PTH"
names(pheno)[32]="Ca"
names(pheno)[33]="P"
names(pheno)[35]="Hb"
names(pheno)[45]="CVD"
cov=pheno[,c(2:19,32,33,36,35,45)]

names(cov)
cov$Smoking=cov$PredictedSmokingStatus
cov$Smoking[which(cov$PredictedSmokingStatus %in% c("Never Smoker","Former Smoker"))] <-0
cov$Smoking[which(cov$PredictedSmokingStatus %in% c("Current Smoker"))] <- 1

cov$P_Group=cov$P
cov$P_Group[which(cov$P < 1.45)] <- 0
cov$P_Group[which(cov$P >= 1.45)] <- 1

table(cov$P_Group)
all=merge(cov,residual,by="Sample_ID")
all=subset(all,Sample_Group == 1)  #只在ESRD患者中进行分析
table(all$P_Group)


result=c()
for ( i in 26:ncol(all)) { fit = glm(all$P_Group~all[,i] + all$Sample_Group +all$Age+all$Gender+all$Smoking,family=binomial); result=as.data.frame(rbind(result,c(colnames(all)[i],coef(summary(fit))[2,c(1,2,4)])));names(result)=c("CpG","Estimate" ,"SE" ,"Pvalue")}

head(result)
result$Pvalue
str(result)
result$Pvalue=as.numeric(result$Pvalue)
result$SE=as.numeric(result$SE)
result$Estimate=as.numeric(result$Estimate)
result$Trait="P"

result=result[order(result$Pvalue),]
head(result)


###############P reuslt + ESRD result
names(result)[2:5]=paste0("P_",colnames(result)[2:5])
final= merge(result,data,by="CpG")
summary(glm(final$P_Estimate ~ final$Bacon_b_meta,family=gaussian))
final_P=subset(final,P_Pvalue < 0.05/52)
names(final_P)[2:5]=c("Estimate","SE","Pvalue","Trait")

cor(final$P_Estimate, final$Bacon_b_meta, method = "spearman")

write.table(final,paste0(Comp,"P.Results.txt"),row.names=F,quote=F)







###################################Hb：115 g/dL
#severe anemia定义:Update on Anemia in ESRD and Earlier Stages of CKD (https://www.ajkd.org/article/S0272-6386(17)31076-4/fulltext)
#severe anemia定义KDIGO Clinical Practice Guideline for Anemia in Chronic (https://kdigo.org/wp-content/uploads/2016/10/KDIGO-2012-Anemia-Guideline-English.pdf)
#Renal anaemia is a frequent complication in patients with chronic kidney disease (CKD). Severe anaemia (haemoglobin <90 g/l) is associated with increased risks of mortality and cardiac complicationshttps://pubmed.ncbi.nlm.nih.gov/23438972/
pheno=read.csv(paste0(Important,"GDPH.1166.Samples.Pheno.txt"),sep="\t")
names(pheno)[36]="PTH"
names(pheno)[32]="Ca"
names(pheno)[33]="P"
names(pheno)[35]="Hb"
names(pheno)[45]="CVD"
cov=pheno[,c(2:19,32,33,36,35,45)]
names(cov)
cov$Smoking=cov$PredictedSmokingStatus
cov$Smoking[which(cov$PredictedSmokingStatus %in% c("Never Smoker","Former Smoker"))] <-0
cov$Smoking[which(cov$PredictedSmokingStatus %in% c("Current Smoker"))] <- 1


cov$Hb_Group[which((cov$Hb < 115  & cov$Gender == 0) |  (cov$Hb < 115  & cov$Gender == 1))] <- 1
cov$Hb_Group[which((cov$Hb > 115  & cov$Gender == 0) |  (cov$Hb > 115  & cov$Gender == 1))] <- 0
table(cov$Hb_Group)
all=merge(cov,residual,by="Sample_ID")
all=subset(all,Sample_Group == 1)  #只在ESRD患者中进行分析
table(all$Hb_Group)


result=c()
for ( i in 26:ncol(all)) { fit = glm(all$Hb_Group~all[,i] +all$Sample_Group+all$Age+all$Gender+all$Smoking,family=binomial); result=as.data.frame(rbind(result,c(colnames(all)[i],coef(summary(fit))[2,c(1,2,4)])));names(result)=c("CpG","Estimate" ,"SE" ,"Pvalue")}

head(result)
result$Pvalue
str(result)
result$Pvalue=as.numeric(result$Pvalue)
result$SE=as.numeric(result$SE)
result$Estimate=as.numeric(result$Estimate)
result$Trait="anemia"

result=result[order(result$Pvalue),]
head(result)


###############Hb reuslt + ESRD result
names(result)[2:5]=paste0("Hb_",colnames(result)[2:5])
final= merge(result,data,by="CpG")
summary(glm(final$Hb_Estimate ~ final$Bacon_b_meta,family=gaussian))
final_Hb=subset(final,Hb_Pvalue < 0.05/52)
names(final_Hb)[2:5]=c("Estimate","SE","Pvalue","Trait")

cor(final$Hb_Estimate, final$Bacon_b_meta, method = "spearman")

write.table(final,paste0(Comp,"Hb.Results.txt"),row.names=F,quote=F)





######整合所有对并发症有影响的位点，P < 0.05/52 (fdr < 0.05)
out="/share/home/liujianjun/ESRD_Project/Table/Output/"
library(dplyr)
output=rbind(final_PTH,final_P,final_Ca,final_Hb)
output1=output %>% 
      select ( CpG,CHR,MAPINFO,Strand,
              Estimate,SE,Pvalue,Trait,
              Bacon_b_DM,Bacon_b_SLEN,Bacon_b_IgAN,Bacon_b_PKD,Bacon_b_HRD,
              Bacon_P_DM, Bacon_P_SLEN, Bacon_P_IgAN, Bacon_P_PKD, Bacon_P_HRD)



# Assuming 'output1' is the data frame provided
# Merge Bacon_b_ columns into association_b with semicolon separation
output1$association_b <- paste(output1$Bacon_b_DM, output1$Bacon_b_SLEN, 
                            output1$Bacon_b_IgAN, output1$Bacon_b_PKD, 
                            output1$Bacon_b_HRD, sep = ";")

# Merge Bacon_P_ columns into association_P with semicolon separation
output1$association_P <- paste(output1$Bacon_P_DM, output1$Bacon_P_SLEN, 
                            output1$Bacon_P_IgAN, output1$Bacon_P_PKD, 
                            output1$Bacon_P_HRD, sep = ";")

# Optionally, remove the original columns if needed
output1 <- output1[, !names(output1) %in% c("Bacon_b_DM", "Bacon_b_SLEN", 
                                       "Bacon_b_IgAN", "Bacon_b_PKD", 
                                       "Bacon_b_HRD", "Bacon_P_DM", 
                                       "Bacon_P_SLEN", "Bacon_P_IgAN", 
                                       "Bacon_P_PKD", "Bacon_P_HRD")]
                                       
write.table(output1,paste0(out,"Table.Association.with.Complication.txt"),quote=F,row.names=F)



###############最后计算各个Complication之间的相关性，保证各自独立
################Pheno file
pheno=read.csv(paste0(Important,"GDPH.1166.Samples.Pheno.txt"),sep="\t")
names(pheno)[36]="PTH"
names(pheno)[32]="Ca"
names(pheno)[33]="P"
names(pheno)[35]="Hb"
names(pheno)[45]="CVD"
cov=pheno[,c(2:19,32,33,36,35,45)]
names(cov)
cov$Smoking=cov$PredictedSmokingStatus
cov$Smoking[which(cov$PredictedSmokingStatus %in% c("Never Smoker","Former Smoker"))] <-0
cov$Smoking[which(cov$PredictedSmokingStatus %in% c("Current Smoker"))] <- 1

all=subset(cov,Sample_Group == 1)  #只在ESRD患者中进行分析
comp=all[,c(19:22)] #选择P PTH Ca Hb 四列并发症
names(comp)=c( "Hypocalcemia", "Hyperphosphatemia", "Hyperparathyroidism", "severe_Anemia")
df_ln <- log(comp) #取对数
df=cor(df_ln, use = "pairwise.complete.obs", method = "spearman")




