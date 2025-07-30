
#############################Perform fixed-effect meta-analysis and the inspection of results. 

bcm <- meta(bc)
save(bcm,file=paste0(index,"Bacon.Subtype.Meta.Results.Rdata"))
save(bc,file=paste0(index,"Bacon_Subtype.Rdata"))


merged_data$Bacon_b_meta=es(bcm)[,6]
merged_data$Bacon_se_meta=se(bcm)[,6]
merged_data$Bacon_P_meta=pval(bcm)[,6]
df <- merged_data 
all=df[order(df$Bacon_P_meta),]


##########亚组effect size 是否一致
all$dir <- ifelse(
  (all$Bacon_b_DM < 0 & all$Bacon_b_SLEN < 0 & all$Bacon_b_HRD < 0 & all$Bacon_b_PKD < 0 & all$Bacon_b_IgAN < 0) |
  (all$Bacon_b_DM > 0 & all$Bacon_b_SLEN > 0 & all$Bacon_b_HRD > 0 & all$Bacon_b_PKD > 0 & all$Bacon_b_IgAN > 0),
  1, 
  0
)


write.table(all,paste0(index,"All_Bacon_Result.txt"),row.names=F,quote=F)

################################计算每个亚型DMP数量

out <- NULL  # Initialize an empty variable to store the results
for (i in c("SLEN", "DM", "IgAN", "HRD", "PKD")) {
  # Dynamically construct the column name
  CpG_column <- paste0("Bacon_P_", i)
  
  # Subset the data where the corresponding column value is less than the threshold
  i_CpG <- subset(all, all[[CpG_column]] < 6.5e-8 / 5)
  
  # Get the number of rows in the subset
  num <- nrow(i_CpG)
  
  # Create the output string
  str <- paste0("The number of significant DMP for ", i, " EWAS is: ", num)
  
  # Append the result to the output
  out <- rbind(out, str)
}
 write.table(out,paste0(index,"Number.of.Subtype.DMPs.txt"),row.names=F)

cut.DM=0.05/4/nrow(subset(all,Bacon_P_DM < 6.528932e-08/5))
cut.SLEN=0.05/4/nrow(subset(all,Bacon_P_SLEN < 6.528932e-08/5))
cut.IgAN=0.05/4/nrow(subset(all,Bacon_P_IgAN < 6.528932e-08/5))
cut.PKD=0.05/4/nrow(subset(all,Bacon_P_PKD < 6.528932e-08/5))
 cut.HRD=0.05/4/nrow(subset(all,Bacon_P_HRD < 6.528932e-08/5))
 

 IgAN.CpG= subset(all,Bacon_P_IgAN < 6.528932e-08/5 & Bacon_P_PKD < cut.PKD & Bacon_P_HRD < cut.HRD & Bacon_P_SLEN < cut.SLEN & Bacon_P_DM < cut.DM)$CpG
 PKD.CpG= subset(all,Bacon_P_PKD < 6.528932e-08/5 & Bacon_P_IgAN < cut.IgAN & Bacon_P_HRD < cut.HRD & Bacon_P_SLEN < cut.SLEN & Bacon_P_DM < cut.DM)$CpG
 HRD.CpG= subset(all,Bacon_P_HRD < 6.528932e-08/5 & Bacon_P_PKD < cut.PKD  &  Bacon_P_IgAN < cut.IgAN & Bacon_P_SLEN < cut.SLEN & Bacon_P_DM < cut.DM)$CpG
 DM.CpG= subset(all,Bacon_P_DM < 6.528932e-08/5 & Bacon_P_PKD < cut.PKD  &  Bacon_P_IgAN < cut.IgAN & Bacon_P_SLEN < cut.SLEN & Bacon_P_HRD < cut.HRD)$CpG
 SLEN.CpG= subset(all,Bacon_P_SLEN < 6.528932e-08/5 & Bacon_P_PKD < cut.PKD & Bacon_P_HRD < cut.HRD & Bacon_P_IgAN < cut.IgAN & Bacon_P_DM < cut.DM)$CpG

can.CpG=unique(c(IgAN.CpG,PKD.CpG,HRD.CpG))
 DMP=subset(all,CpG %in% can.CpG &   dir==1)

data=DMP

write.table(DMP,paste0(index,"ESRD.Shared.DMPs.txt"),row.names=F,quote=F)


dat=data.frame(DM=length(DM.CpG),SLEN=length(SLEN.CpG),IgAN=length(IgAN.CpG),PKD=length(PKD.CpG),HRD=length(HRD.CpG),meta=nrow(DMP))
 rownames(dat)="count"


####################计算异质性
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
library(metafor)

data=fread(paste0(index,"ESRD.Shared.DMPs.txt"))

meta_results <- data.frame(
  CpG = character(0),
  I2 = numeric(0),
  Q = numeric(0),
  Phet= numeric(0),
  stringsAsFactors = FALSE
)

# 循环计算每个CpG的meta分析和异质性（I²）
for(i in 1:nrow(data)) {
  # 为每个CpG提取相应的beta和标准误
  yi <- c(data$Bacon_b_DM[i],data$Bacon_b_SLEN[i],data$Bacon_b_IgAN[i], data$Bacon_b_PKD[i], data$Bacon_b_HRD[i])
  sei <- c(data$Bacon_se_DM[i],data$Bacon_se_SLEN[i],data$Bacon_se_IgAN[i], data$Bacon_se_PKD[i], data$Bacon_se_HRD[i])
  
  # 执行固定效应meta分析
  meta_result <- rma(yi, sei = sei, method = "FE")
  
  # 将每个CpG的结果存储到meta_results数据框中
  meta_results <- rbind(meta_results, data.frame(
    CpG = data$CpG[i],
    I2 = meta_result$I2,
    Q = meta_result$QE,
    Phet=meta_result$QEp
  ))
}

result=merge(data,meta_results,by="CpG")
result=result[order(result$Bacon_P_meta),]
write.table(result,paste0(index,"DMP.Heterogeneity.txt"),quote=F,row.names=F)


########################读取EPIC芯片注释信息
Important="/share/home/liujianjun/ESRD_Project/Important_files/"
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"

library(data.table)
library(bacon)
anno=fread(paste0(Important,"Slim.MethylationEPIC_v-1-0_B4.txt"))
result=fread(paste0(index,"DMP.Heterogeneity.txt"))

sub=subset(anno,IlmnID %in% result$CpG)
names(sub)[1]="CpG"
#sub=sub[,c(1,6,7)]
file=merge(result,sub,by="CpG")

########################SG 队列进行验证
sg=fread("/share/home/liujianjun/ESRD_Project/Part2.Subtype_ESRD_EWAS/KPTH/step3.annotation.results.files/DM.linear_result_all.txt")
library(bacon)
es=sg$Estimate
se=sg$SE
bc=bacon(NULL,es,se)
sg_bacon=data.frame(CpG=sg$CpG,bacon_Estimate=es(bc),bacon_SE=se(bc),bacon_P=pval(bc))
names(sg_bacon)=c("CpG","Replication_b","Replication_se","Replication_P")
rep=subset(sg_bacon,CpG %in% file$CpG  & Replication_P < 0.05/60)
all=merge(file,rep,by="CpG")

###########发现组和验证组方向是否一致
all$dir_rep <- ifelse(
  (all$Replication_b <0 & all$Bacon_b_DM < 0 & all$Bacon_b_SLEN < 0 & all$Bacon_b_HRD < 0 & all$Bacon_b_PKD < 0 & all$Bacon_b_IgAN < 0) |
  (all$Replication_b >0 &all$Bacon_b_DM > 0 & all$Bacon_b_SLEN > 0 & all$Bacon_b_HRD > 0 & all$Bacon_b_PKD > 0 & all$Bacon_b_IgAN > 0),
  1, 
  0
)

write.table(all,paste0(index,"Full.Rep.DMP.Results"),row.names=F,quote=F)



