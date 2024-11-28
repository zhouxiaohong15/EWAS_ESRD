################This script was aimed to replicated our EWAS meta analysis result by integrating 4 reported studys.The replicated threshold P-value is 10*e-6.

################读取差异甲基化位点
DMP=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.Annotate.DMP.txt",header=T)
names(DMP)[1]="probeID"

###############Step 1. CKD 相关的位点来进行验证，共纳入2篇文献
######第一篇 1.NC_2017
CKD1=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/1.Kidney_EWAS_literatrue/1.NC_2017/CKD/META_EWAS_CKD.TBL",header=T)
CKD1sig=subset(CKD1,P.value < 0.000001)
CKD1_DMP=subset(CKD1sig, TargetID %in% DMP$probeID)

#######第二篇 3.NC_2021
CKD2=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/1.Kidney_EWAS_literatrue/3.NC_2021/CKD.csv",header=T,sep="\t")
CKD2_sig=subset(CKD2, P.value < 0.000001)
CKD2_DMP=subset(CKD2_sig, probeID  %in% DMP$probeID)

names(CKD2_DMP)[2]="probeID"
names(CKD1_DMP)[1]="probeID"

######### 合并所有CKD replication到的位点
CKD_DMP=merge(CKD1_DMP,CKD2_DMP,by="probeID",all=T)
DMP_CKD=subset(DMP,probeID %in% CKD_DMP$probeID)
CKD_final=merge(DMP_CKD,CKD_DMP,merge="probeID")
write.table(CKD_final,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/3.Replicated_Results/CKD_Replication_results.txt")


####################Step 1. eGFR相关的位点来进行验证，共纳入5篇文献,阈值为原文所设定的显著性阈值如FDR < 0.05, P < 5E-7 以及P < 5E-5; 注意：这都是文章所设定的P value.
#########第一篇：1.NC_2017;P meta<1e-5；450K
eGFR1=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/1.Kidney_EWAS_literatrue/1.NC_2017/eGFR/META_EWAS_eGFR.TBL",header=T)
names(eGFR1)[2]="probeID"
head(eGFR1)
eGFR1_sig=subset(eGFR1, P.value < 0.00001)
eGFR1_DMP=subset(eGFR1_sig, probeID %in% DMP$probeID)
overlap.eGFR1.DMP=nrow((eGFR1_DMP))
eGFR1.len=length(intersect(DMP$probeID,eGFR1$probeID))
eGFR1.Journal="NC_2017"

########第二篇：2.GM_2021; 原文显著定义：FDR < 0.05;450K,EPIC
eGFR2=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/1.Kidney_EWAS_literatrue/2.GM_2021/eGFR_DMP_FDR0.05.txt",header=T,sep="\t")
dim(eGFR2)
eGFR2$Paper="2021_GM"
names(eGFR2)[1] = "probeID"
eGFR2_sig=eGFR2
eGFR2_DMP=subset(eGFR2_sig, probeID %in% DMP$probeID)
overlap.eGFR2.DMP=nrow((eGFR2_DMP))

eGFR2.full=fread("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/1.Kidney_EWAS_literatrue/2.GM_2021/eGFR_discovery_transethnic.txt")
names(eGFR2.full)[2]="probID"
eGFR2.len=length(intersect(eGFR2.full$probID, DMP$probeID))
overlap.eGFR2.DMP=nrow((eGFR2_DMP)
eGFR2.Journal="GM_2021"

#######第三篇：3.NC_2021; combined p-value<1.1E-7; 450K,EPIC
eGFR3=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/1.Kidney_EWAS_literatrue/3.NC_2021/eGFR.csv",header=T,sep="\t")
head(eGFR3)
eGFR3_sig=subset(eGFR3, P.value < 0.00000011)
eGFR3_DMP=subset(eGFR3_sig, probeID %in% DMP$probeID)

overlap.eGFR3.DMP=nrow((eGFR3_DMP))
eGFR3_len=length(intersect(DMP$probeID,eGFR3$probeID))  ##可评估的ESRD DMP位点数量
eGFR3.Journal="NC_2021"

######第四篇：4.2020_PNAS; P < 0.00005; EPIC
eGFR4=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/1.Kidney_EWAS_literatrue/4.2020_PNAS/eGFR_EWAS.txt",header=T,sep="\t")
head(eGFR4)
names(eGFR4)[2]="probeID"
head(eGFR4)
eGFR4_sig=subset(eGFR4,Pvalue < 0.00005)
eGFR4_DMP=subset(eGFR4_sig,probeID %in% DMP$probeID)

overlap.eGFR4.DMP=nrow((eGFR4_DMP))
eGFR4_len=length(intersect(DMP$probeID,eGFR4$probeID))  ##可评估的ESRD DMP位点数量
eGFR4.Journal="PNAS_2020"

######第五篇：5.2023_NC；FDR = 0.05;  450K
eGFR5=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/1.Kidney_EWAS_literatrue/8.2023_NC/eGFR.baseline.EWAS.FDR0.05.txt",header=T,sep="\t")
names(eGFR5)[1]="probeID"
eGFR5_sig=eGFR5
eGFR5_DMP=subset(eGFR5_sig,probeID %in% DMP$probeID)

overlap.eGFR5.DMP=nrow((eGFR5_DMP))
eGFR5_len=length(intersect(DMP$probeID,eGFR5$probeID))  ##可评估的ESRD DMP位点数量
eGFR5.Journal="NC_2021"

########合并所有五篇验证到的eGFR结果
new_eGFR=merge(eGFR1_DMP,eGFR2_DMP,by="probeID",all=T)
new_eGFR1=merge(eGFR3_DMP,eGFR4_DMP,by="probeID",all=T)
new_eGFR2=merge(new_eGFR,eGFR5_DMP,by="probeID",all=T)
new_eGFR_final=merge(new_eGFR1,new_eGFR2,by="probeID",all=T)

#######合并: 验证到的eGFR结果和我们的meta DMP结果
head(new_eGFR_final)
DMP_eGFR_final=subset(DMP,probeID %in% new_eGFR_final$probeID)
dim(DMP_eGFR_final)
final_eGFR=merge(DMP_eGFR_final,new_eGFR_final,by="probeID",all=T)

write.table(final_eGFR,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/3.Replicated_Results/eGFR_Significant_Replication_results.txt",row.names=F)


#####################################################Step2.1 统计：有多少我们的DMP在过往使用450K芯片的eGFR研究中，同样可以达到全基因组显著。
DMP=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.Annotate.DMP.txt",header=T)
names(DMP)[1]="probeID"
#file=as.data.frame(fread("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartH.Target_gene_annotation/MethylationEPIC_v-1-0_B4.csv"))
#Array = file[,c(2,41)]
#write.table(Array,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartH.Target_gene_annotation/Array450K.txt",row.names=F)

file=as.data.frame(fread("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartH.Target_gene_annotation/Array450K.txt"))
names(file)[1]="probeID"
Array=subset(file,Methyl450_Loci == "TRUE")

###########在既往450K研究中有信息的DMP个数
DMP450=length(intersect(DMP$probeID,Array$probeID))

##########在既往450K研究中全基因组显著的DMP个数,共三篇使用450K：1）第一篇：1.NC_2017;P meta<1e-5；450;  2)3.NC_2021; combined p-value<1.1E-7; 450K,EPIC; 3)5.2023_NC；FDR = 0.05;  450K
DMP450.sig= length(unique( c(eGFR1_DMP$probeID,eGFR3_DMP$probeID,eGFR5_DMP$probeID)))

#######450K显著能验证的比例=在既往450K研究中全基因组显著的DMP个数/在既往450K研究中有信息的DMP个数
Ratio.450K.sig=DMP450.sig/DMP450

#####################################################Step2.2统计：有多少我们的DMP在过往使用850K芯片的eGFR研究中，达到全基因组显著
DMP=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.Annotate.DMP.txt",header=T)
names(DMP)[1]="probeID"

##########在既往850K研究中全基因组显著的DMP个数，共两篇使用850K: 1) 2.GM_2021; 原文显著定义：FDR < 0.05;450K,EPIC; 2) 4.2020_PNAS; P < 0.00005; EPIC
DMP850.sig= length(unique( c(eGFR2_DMP$probeID,eGFR4_DMP$probeID)))

row=unique(c(eGFR2$probeID,eGFR4$probeID))
DMP850=length(intersect(DMP$probeID,row))

#######850K显著能验证的比例=在既往850K研究中全基因组显著的DMP个数/在既往850K研究中有信息的DMP个数
Ratio.850K.sig=DMP850.sig/DMP850





###############################Step 3. eGFR相关的位点来进行验证,注意：与step2不同，现在阈值P value < 0.05.

DMP=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.Annotate.DMP.txt",header=T)
names(DMP)[1]="probeID"

#########第一篇：1.NC_2017;P meta<0.05；450K. Full summary statistics: 430136个位点
eGFR1=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/1.Kidney_EWAS_literatrue/1.NC_2017/eGFR/META_EWAS_eGFR.TBL",header=T)
names(eGFR1)[2]="probeID"
head(eGFR1)
eGFR1_sig=subset(eGFR1, P.value < 0.05)
nrow(eGFR1_sig)
eGFR1_DMP=subset(eGFR1_sig, probeID %in% DMP$probeID)
nrow(eGFR1_DMP)

########第二篇：2.GM_2021; 原文显著定义：FDR < 0.05;450K,EPIC. 共纳入10个研究队列。6个用EPIC，4个用450K.Full summarry statistics: 655797个位点
eGFR2=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/1.Kidney_EWAS_literatrue/2.GM_2021/eGFR_discovery_transethnic.txt",header=T,sep="\t")
dim(eGFR2)
names(eGFR2)[2] = "probeID"
eGFR2_sig=subset(eGFR2, P_value < 0.05)
nrow(eGFR2_sig)
eGFR2_DMP=subset(eGFR2_sig, probeID %in% DMP$probeID)
nrow(eGFR2_DMP)

#######第三篇：3.NC_2021; combined p-value<1.1E-7; 450K,EPIC; full summary:441871个位点
eGFR3=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/1.Kidney_EWAS_literatrue/3.NC_2021/eGFR.csv",header=T,sep="\t")
head(eGFR3)
eGFR3_sig=subset(eGFR3, P.value < 0.05)
nrow(eGFR3_sig)
eGFR3_DMP=subset(eGFR3_sig, probeID %in% DMP$probeID)
nrow(eGFR3_DMP)


######第四篇：4.2020_PNAS; P < 0.00005; EPIC; full summary: 778515个位点
eGFR4=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/1.Kidney_EWAS_literatrue/4.2020_PNAS/eGFR_EWAS.txt",header=T,sep="\t")
head(eGFR4)
names(eGFR4)[2]="probeID"
head(eGFR4)
eGFR4_sig=subset(eGFR4,Pvalue < 0.05)
nrow(eGFR4_sig)
eGFR4_DMP=subset(eGFR4_sig,probeID %in% DMP$probeID)
nrow(eGFR4_DMP)
########合并所有四篇验证到的eGFR结果（P < 0.05)
new_eGFR=merge(eGFR1_DMP,eGFR2_DMP,by="probeID",all=T)
new_eGFR1=merge(eGFR3_DMP,eGFR4_DMP,by="probeID",all=T)
new_eGFR_final=merge(new_eGFR,new_eGFR1,by="probeID",all=T)

#######合并: 验证到的eGFR结果和我们的meta DMP结果
head(new_eGFR_final)
DMP_eGFR_final=subset(DMP,probeID %in% new_eGFR_final$probeID)
dim(DMP_eGFR_final)
final_eGFR=merge(DMP_eGFR_final,new_eGFR_final,by="probeID",all=T)

write.table(final_eGFR,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/3.Replicated_Results/eGFR_Nominal_Replication_results.txt",row.names=F)


#####################################################Step3.1 统计：有多少我们的DMP在过往使用450K芯片的eGFR研究中，同样可以达到nominal显著。P<0.05
DMP=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.Annotate.DMP.txt",header=T)
names(DMP)[1]= "probeID"
file=as.data.frame(fread("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartH.Target_gene_annotation/Array450K.txt"))
names(file)[1]="probeID"
Array=subset(file,Methyl450_Loci == "TRUE")

###########在既往450K研究中有信息的DMP个数
DMP450= length(intersect(DMP$probeID,Array$probeID))

##########在既往450K研究中nominal显著的DMP个数,共两篇使用450K且有full summary statistics：1）第一篇：1.NC_2017  2)3.NC_2021
DMP450.nominal= length(unique( c(eGFR1_DMP$probeID,eGFR3_DMP$probeID)))

#######450K显著能验证的比例=在既往450K研究中nominal显著的DMP个数/在既往450K研究中有信息的DMP个数
Ratio.450K.nominal=DMP450.nominal/DMP450

#####################################################Step2.2统计：有多少我们的DMP在过往使用850K芯片的eGFR研究中，达到nominal显著
DMP=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Meta_analysis/01.Bacon_meta_ALL_results/Bacon.Annotate.DMP.txt",header=T)

##########在既往850K研究中全基因组显著的DMP个数，共两篇使用850K: 1) 2.GM_2021;  2) 4.2020_PNAS; 
DMP850.nominal= length(unique( c(eGFR2_DMP$probeID,eGFR4_DMP$probeID)))

row=unique(c(eGFR2$probeID,eGFR4$probeID))
DMP850=length(intersect(DMP$probeID,row))
#######850K显著能验证的比例=在既往850K研究中全基因组显著的DMP个数/在既往850K研究中有信息的DMP个数
Ratio.850K.nominal=DMP850.nominal/DMP850




##############################Step 4. eGFR slope 相关的位点来进行验证，共纳入3篇文献

######第一篇：2020_PNAS;P < 0.00005; EPIC
slope=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation//1.Kidney_EWAS_literatrue/4.2020_PNAS/eGFR_slope_EWAS.txt",header=T)
head(slope)
names(slope)[2]= "probeID"
head(slope)
slope_sig=subset(slope, Pvalue < 0.00005 & Pvalue != 0)
slope_len=length(intersect(DMP$probeID,slope$probeID))  ##可评估的ESRD DMP位点数量

nrow(slope_sig)
nrow(slope_len)
slope_DMP=subset(slope_sig,probeID %in% DMP$probeID)
nrow(slope_DMP)

######第二篇：2018_KI;FDR < 0.05;450K
slope1=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation//1.Kidney_EWAS_literatrue/6.2018_KI/eGFR.slope.DMP.FDR0.05.txt",header=T,sep="\t")
names(slope1)[1]= "probeID"
slope1_sig=slope1
length(intersect(DMP$probeID, slope1$probeID))  ###可评估的ESRD DMP位点数量
slope1_DMP=subset(slope1_sig,probeID %in% DMP$probeID)

#######第三篇：2019_NC;FDR < 0.05;450K
slope2=read.csv("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation//1.Kidney_EWAS_literatrue/7.2019_NC/eGFR.slope.EWAS.FDR0.05.txt",header=T,sep="\t")
names(slope2)[1]= "probeID"
slope2_sig=slope2
length(intersect(DMP$probeID, slope2$probeID)) 
slope2_DMP=subset(slope2_sig,probeID %in% DMP$probeID)
dim(slope2_DMP)

######第四篇：2023_NC;FDR < 0.05; 450K
slope3=read.csv("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation//1.Kidney_EWAS_literatrue/8.2023_NC/eGFR.slope.EWAS.FDR0.05.txt",header=T,sep="\t")
names(slope3)[1]= "probeID"
slope3_sig=slope3
length(intersect(DMP$probeID, slope3$probeID))
slope3_DMP=subset(slope3_sig,probeID %in% DMP$probeID)
dim(slope3_DMP)



#########合并: 验证到的eGFR slope结果和我们的meta DMP结果
DMP_slope=subset(DMP,probeID %in% slope_DMP$probeID)
slope_final=merge(DMP_slope,slope_DMP,slope1_DMP,slope1_DMP,slope3_DMP,by="probeID")
#未验证到结果
#write.table(slope_final,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/3.Replicated_Results/eGFR_Slope_Replication_results.txt",row.names=F)


###################Step5.  eGFR slope相关的位点来进行验证,注意：与step4不同，现在阈值P value < 0.05. 注意：只有一篇文章有full summary信息。

######第一篇：2020_PNAS;P < 0.00005; EPIC
slope=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation//1.Kidney_EWAS_literatrue/4.2020_PNAS/eGFR_slope_EWAS.txt",header=T)
head(slope)
names(slope)[2]= "probeID"
names(DMP)[1]="probeID"
head(slope)
slope_sig=subset(slope, Pvalue < 0.05)
slope_DMP=subset(slope_sig,probeID %in% DMP$probeID)


#########合并: 验证到的eGFR slope结果和我们的meta DMP结果

DMP_slope=subset(DMP,probeID %in% slope_DMP$probeID)
slope_final=merge(DMP_slope,slope_DMP,by="probeID")

write.table(slope_final,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/3.Replicated_Results/eGFR_Slope_Nominal_Replication_results.txt",row.names=F)

#################在EPIC研究中有信息的DMP位点个数
DMP.slope.850k=length(intersect(slope$probeID,DMP$probeID))

#################在EPIC研究中nominal显著的位点个数
DMP.slope.850k.nominal=nrow(slope_DMP)

##################验证到 的比例
Ratio.slope.850k.nominal=DMP.slope.850k.nominal/DMP.slope.850k


#########################################################整合所有eGFR以及eGFR Slope的统计结果
#####注意，对于eGFR slope在进行nomial验证时，只有一篇文献有full summary结果。
Name=c("eGFR.450k","eGFR.450k.Significant.Ratio","eGFR.850k","eGFR.850k.Significant.Ratio","DMP450K.nominal","eGFR.450K.nominal.Ratio","DMP850K.nominal","eGFR.850k.Nominal.Ratio","eGFR.slope.Significant","eGFR.slope.850k","eGFR.slope.850k.Nominal.Ratio")
Num=c(DMP450,Ratio.450K.sig,DMP850,Ratio.850K.sig,DMP450.nominal,Ratio.450K.nominal,DMP850.nominal,Ratio.850K.nominal,0,DMP.slope.850k,Ratio.slope.850k.nominal)

frame=data.frame(Name=Name,Num=Num)

write.table(frame,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/3.Replicated_Results/Summary.DMP.Replication.txt",row.names=F)




##############################Step 4. T1D ESRD 相关的位点来进行验证，共纳入1篇文献
DMP=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Meta_analysis/04.Bacon_meta_ALL_results/Bacon.meta.DMP.txt",header=T)

########################ESRD数据集是指case为107例肾移植和透析的ESRD患者，control为253例TIDM患者
ESRD=read.csv("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation//1.Kidney_EWAS_literatrue/5.2021_CE/T1D_ESRD_360samples.txt",header=T,sep="\t")
names(ESRD)[c(1,11,12,13)]=c("probeID","P.value","Fold.Change","Description") 
ESRD=ESRD[,c(1,11,12,13)]
ESRD$case <- "Transplant_Dialysis"
ESRD_DMP=subset(ESRD,probeID %in% DMP$probeID)


######################transplant数据集是指case为73例肾移植的ESRD患者，control为253例TIDM患者
transplant=read.csv("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation//1.Kidney_EWAS_literatrue/5.2021_CE/T1D_ESRD_transplant_326samples.txt",header=T,sep="\t")
names(transplant)[c(1,11,12,13)]=c("probeID","P.value","Fold.Change","Description") 
transplant=transplant[,c(1,11,12,13)]
transplant$case <- "Transplant"
transplant_DMP=subset(transplant,probeID %in% DMP$probeID)

#########################合并所有ESRD validation到的位点
all_DMP=merge(ESRD_DMP,transplant_DMP,by="probeID",all=T)
DMP_all=subset(DMP,probeID %in% all_DMP$probeID)
all_final=merge(DMP_all,all_DMP,merge="probeID")
write.table(all_final,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/3.Replicated_Results/ESRD_Replication_results.txt",row.names=F)



######################Step5 读取UACR的文章数据;P < 1.1e-7;
UACR=read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/1.Kidney_EWAS_literatrue/3.NC_2021/UACR.csv",header=T)

UACR_sig=subset(UACR, P.value < 0.00000011)
UACR_DMP=subset(UACR_sig, probeID %in% DMP$probeID)
names(UACR_DMP)[2]="probeID"
DMP_UACR=subset(DMP,probeID %in% UACR_DMP$probeID)
UACR_final=merge(DMP_UACR,UACR_DMP,merge="probeID")

write.table(UACR_final,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartF.External_Validation/3.Replicated_Results/UACR_Replication_results.txt",row.names=F)

