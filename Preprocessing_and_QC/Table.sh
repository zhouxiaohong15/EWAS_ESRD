#Table 1: EWAS results of subtype ESRD
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
Important="/share/home/liujianjun/ESRD_Project/Important_files/"
table="/share/home/liujianjun/ESRD_Project/EWAS_Manuscript_Table/"



library(data.table)
library(dplyr)
data=fread(paste0(index,"All_Bacon_Result.txt"))

anno=fread(paste0(Important,"Slim.MethylationEPIC_v-1-0_B4.txt"))
names(anno)[1]="CpG"
all=data

################################计算每个亚型DMP数量,并列出每个亚型的DMP信息整理为table
out <- NULL  # Initialize an empty variable to store the results
result= c()
for (i in c("SLEN", "DM", "IgAN", "HRD", "PKD")) {
  # Dynamically construct the column name
  P_column <- paste0("Bacon_P_", i)
  b_column <- paste0("Bacon_b_", i)
  se_column <- paste0("Bacon_se_", i)
  
  # Subset the data where the corresponding column value is less than the threshold
  i_CpG <- subset(all, all[[P_column]] < 6e-8 / 5)
  i_CpG = i_CpG %>% select(CpG, b_column,se_column,P_column,)
  names(i_CpG)=c("CpG","b","se","P")
  i_CpG$Group=i
  frame=merge(i_CpG,anno,by="CpG")
  frame <- bind_rows(data.frame(CpG = i), frame)
  
  # Get the number of rows in the subset
  num <- nrow(i_CpG)
  
  # Create the output string
  str <- paste0("The number of significant DMP for ", i, " EWAS is: ", num)
  
  # Append the result to the output
  out <- rbind(out, str)
  result=as.data.frame(rbind(result,frame))
}

write.table(result,paste0(table,"Table2.ALL.Subtype.DMP"),quote=F,row.names=F)



########Table 2: Discovery cohort and Replication cohort
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
Important="/share/home/liujianjun/ESRD_Project/Important_files/"
table="/share/home/liujianjun/ESRD_Project/EWAS_Manuscript_Table/"

library(data.table)
library(dplyr)

sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")
names(sentinels)[1]="CpG"

#Reported kidney trait related CpGs
DMP=fread(paste0(Important,"Kidney_Trait_CpGs.txt"))
 names(DMP)[1]="CpG"

intersect(sentinels$CpG,DMP$CpG)
dt=unique(subset(DMP,CpG %in% sentinels$CpG))
dt_agg <- dt[, .(
  beta = paste(beta, collapse = ";"),
  Trait = paste(Trait, collapse = ";"),
  PMID = paste(PMID, collapse = ";")
), by = CpG]


anno=fread(paste0(Important,"Slim.MethylationEPIC_v-1-0_B4.txt"))
names(anno)[1]="CpG"
anno=anno[,1:3]

all=merge(sentinels,dt_agg,by="CpG",all.x=TRUE)
final=merge(all,anno,by="CpG")

write.table(final,paste0(table,"Table3.Discovery_Replication.txt"),quote=F,row.names=F,sep=",")









##############Gene Prioritation : tro CpG --- Gene ---eGFR ----ESRD
# Load necessary libraries
library(tidyr)
library(data.table)
library(dplyr)

# Define file paths
base="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Database/"
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
out="/share/home/liujianjun/ESRD_Project/Table/Output/"
Imp="/share/home/liujianjun/ESRD_Project/Important_files/"
fi="/share/home/liujianjun/ESRD_Project/Part13.DMP_to_Protein/Final_Results/"
inp="/share/home/liujianjun/ESRD_Project/Table/Input/"
out="/share/home/liujianjun/ESRD_Project/Table/Output/"



gene <- fread(paste0(out, "Final.ESRD_DMP_Gene.txt"), sep = "\t")

# Load data files
eqtm2015 <- fread(paste0(Imp, "2015_09_02_cis_eQTMsFDR0.05-CpGLevel_blood.txt"))

# Subset and process eqtm2015 data
eqtm2015_sub <- eqtm2015[, c(2, 3, 4, 17, 20)]
names(eqtm2015_sub)[5] <- "Beta"

# Split Beta column
eqtm2015_sub <- eqtm2015_sub %>%
  separate(Beta, into = c("Beta1", "Beta2", "Beta3", "Beta4"), sep = ";") %>%
  mutate(
    Beta1 = gsub("\\s*\\(.*\\)", "", Beta1),
    Beta2 = gsub("\\s*\\(.*\\)", "", Beta2),
    Beta3 = gsub("\\s*\\(.*\\)", "", Beta3),
    Beta4 = gsub("\\s*\\(.*\\)", "", Beta4)
  )

# Retain Beta1 column and clean up
eqtm2015_sub$Beta <- eqtm2015_sub$Beta1
eqtm2015_sub <- eqtm2015_sub %>% select(-Beta1)
names(eqtm2015_sub)[1] <- "CpG"
eqtm2015_sub$P <- eqtm2015$PValue
eqtm2015_sub <- eqtm2015_sub[, c(1, 4, 8, 9)]
names(eqtm2015_sub) <- c("CpG", "Gene", "beta", "P")

# Load additional datasets
eqtm2018_gtp <- fread(paste0(Imp, "2018_BMC_GTP.cohort_blood_eQTM.txt"))
eqtm2018_gtp_sub <- eqtm2018_gtp[, c(2, 7, 15, 13)]
names(eqtm2018_gtp_sub) <- c("CpG", "Gene", "beta", "P")

eqtm2018_mesa <- fread(paste0(Imp, "2018_BMC_MESA.cohort_blood_eQTM.txt"))
eqtm2018_mesa_sub <- eqtm2018_mesa[, c(1, 7, 15, 13)]
names(eqtm2018_mesa_sub) <- c("CpG", "Gene", "beta", "P")

eqtm2022_kidney <- fread(paste0(Imp, "2022_NG_Kidney.eQTM.FDR0.05.txt"))
eqtm2022_kidney_sub <- eqtm2022_kidney[, c(1, 3, 5, 7)]
names(eqtm2022_kidney_sub) <- c("CpG", "Gene", "beta", "P")

# Load cis/trans EQTM 2023 data
cis_eqtm_2023 <- fread(paste0(Imp, "2023_SR_FHS_Top10000cis_eQTM.txt"), sep = "\t")
cis_eqtm_2023$class <- "cis"
trans_eqtm_2023 <- fread(paste0(Imp, "2023_SR_FHS_Top10000trans_eQTM.txt"), sep = "\t")
trans_eqtm_2023$class <- "trans"
eqtm_2023 <- rbind(cis_eqtm_2023, trans_eqtm_2023)
eqtm2023_sub=eqtm_2023[,c(1,5,11,14)]
# Clean up eqtm_2023 data

names(eqtm2023_sub) <- c("CpG", "Gene", "beta", "P")

# Rename columns for merging

eqtm2015_sub$method = "eqtm2015"
eqtm2018_gtp_sub$method = "eqtm2018_gtp"
eqtm2018_mesa_sub$method = "eqtm2018_mesa"
eqtm2022_kidney_sub$method = "eqtm2022_kidney"
eqtm2023_sub$method = "eqtm2023"

names(eqtm2015_sub)[2:4] <- c("Gene", "beta", "P")
names(eqtm2018_gtp_sub)[2:4] <- c("Gene", "beta", "P")
names(eqtm2018_mesa_sub)[2:4] <- c("Gene", "beta", "P")
names(eqtm2022_kidney_sub)[2:4] <- c("Gene", "beta", "P")
names(eqtm2023_sub)[2:4] <- c("Gene", "beta", "P")

eqtm <- rbind(eqtm2015_sub, eqtm2018_gtp_sub, eqtm2018_mesa_sub, eqtm2022_kidney_sub, eqtm2023_sub)


# Save the merged data to a file
write.table(eqtm, paste0(inp, "All_Published_eQTM_summary.txt"), row.names = FALSE, sep = "\t", quote = FALSE)


##加载TWAS 结果
base="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Database/"

twas=fread(paste0(base,"kidney_functon_TWAS.txt"))
 twas_sub=data.frame(twas_Tarit=twas$V1,twas_Tissue=twas$V2,Gene=twas$V4,twas_beta=twas$V11,twas_Pvalue=twas$V16)
 twas_sub=subset( twas_sub, twas_Pvalue <= 0.05/(nrow(gene))) ####这个阈值是代表一共注释到这个基因数量：nrow(gene)

merged_twas=twas_sub

write.table(merged_twas, paste0(inp, "All_eGFR_TWAS_summary.txt"), row.names = FALSE, sep = "\t", quote = FALSE)




###将eQTM数据和TWAS数据，以gene作为键，连接两个数据


######################加载具有trio的CpG位点
gene <- fread(paste0(out, "Final.ESRD_DMP_Gene.txt"), sep = "\t")
DMP=read.csv(paste0(index,"Full.Rep.DMP.Results"),sep=" ")
DMP=DMP[,c(1,24,25,17,19)]

merged_twas=fread( paste0(inp, "All_eGFR_TWAS_summary.txt"),sep="\t") #twas
eqtm=fread(paste0(inp, "All_Published_eQTM_summary.txt"),sep="\t")  #eqtm


# 提取source列中包含"eqtm"的行
filtered_gene <- gene %>%
  filter(grepl("eqtm", source)) 

sub=subset(filtered_gene, TWAS == "*")

df=subset(eqtm,CpG %in% sub$CpG  & Gene %in% sub$gene)
df=df[order(df$CpG),]
ID=which(df$method == "eqtm2022_kidney")
blood=df[-ID,]


##########合并DMP的EWAS值和meQTM值
data=merge(blood,DMP,by="CpG",all.x=T)

########然后合并上述与twas的结果

all=merge(data,merged_twas,by="Gene")
all$eqtm_beta=as.numeric(all$eqtm_beta)
all$Consistant_effect=all$eqtm_beta * all$twas_beta * all$Bacon_b_meta  



write.table(all,paste0(out,"Trio_DMP_Gene_eGFR.txt"),sep="\t",quote=F,row.names=F)








#######explore the association between ESRD shared DMP and kidney trait
Imp="/share/home/liujianjun/ESRD_Project/Important_files/"
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
out="/share/home/liujianjun/ESRD_Project/Table/Output/"


library(data.table)

P=fread(paste0(Imp,"Kidney_Trait_CpGs_Pvalue.txt"))  ##既往发表的EWAS的DMP的P值
names(P)[c(2,4)]=c("Pval","Trait")
beta=fread(paste0(Imp,"Kidney_Trait_CpGs.txt"))  ##既往发表的EWAS的DMP的beta值，注意beta表格和P表格的长度不一致，可能有部分没有粘贴上去
data=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")


list=intersect(data$CpG,P$CpG)  #能overlap的DMP

sub_beta=subset(beta,CpG %in% list)
sub_P=subset(P,CpG %in% list)
df=merge(sub_beta,sub_P,all=TRUE,by=c("CpG","Trait","PMID"))
df[is.na(df)] <- NA


sub=data[,c(1,2,5,8,11,14,24,25,27)] ##提取我们研究中的ESRD效应值

dt=merge(df,sub,by="CpG",all.x=TRUE)
all=unique(dt)

all$beta=as.numeric(all$beta)
all$Pval=as.numeric(all$Pval)

##长格式变为宽格式
all_wide <- dcast(
  all,
  CpG  + Bacon_b_DM + Bacon_b_SLEN + Bacon_b_IgAN + Bacon_b_PKD + 
  Bacon_b_HRD + CHR + MAPINFO + UCSC_RefGene_Name ~ Trait,
  value.var = c("beta", "Pval", "PMID"),
  fun.aggregate = function(x) paste(unique(x), collapse = ";")
)

all_wide[is.na(all_wide)] <- NA
all_wide[all_wide == ""] <- NA
all_wide[all_wide == -99] <- NA

write.table(all_wide,paste0(out,"Table.Overlapped_eGFR_DMP.txt"),quote=F,row.names=F,sep="\t")   ####注意有部分beta的原始输入文件是残缺的，所以结果导出来有很多的残缺


##合并Overlap DMP 位点与association的结果
out="/share/home/liujianjun/ESRD_Project/Table/Output/"
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
path1="/share/home/liujianjun/ESRD_Project/Table/Input/"

library(dplyr)
library(data.table)

eGFR=fread(paste0(out,"Table.Overlapped_eGFR_DMP.txt"),fill=TRUE,sep="\t")
eGFR=subset(eGFR, PMID_eGFR != "NA")
data=read.csv(paste0(index,"Full.Rep.DMP.Results"),sep=" ")
eGFR=unique(eGFR[,1])

output=subset(data,CpG %in% eGFR$CpG)

output1=output %>% 
      select ( CpG,CHR,MAPINFO,Strand,UCSC_RefGene_Name,
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


write.table(output1,paste0(out,"eGFR_overlap.DMP_ESRD_association.txt"),row.names=F,quote=F,sep="\t")

##合并Overlap 的DMP和 prioritized gene


#Prioritization gene
setwd(path1)
gene=fread("Final.ESRD_DMP_Gene.txt")
result=output1

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
write.table(result1,paste0(out,"eGFR_overlap.DMP_result_eQTMgene.txt"),row.names=F,quote=F)

###选取enhancer 位点
new_dt2=subset(new_dt1,Enhancer_linked_gene > 0 )


new_dt3=setDT(new_dt2)[, .(gene = paste(gene, collapse = ";")), by = CpG]
names(new_dt3)[2]="Enhancer_linked_gene"

##合并
result1=merge(result,new_dt3,by="CpG",all.x=T)
 result1=result1[,c("CpG","Enhancer_linked_gene")]
 write.table(result1,paste0(out,"eGFR_overlap.DMP_result_Enhancer_linked_gene.txt"),row.names=F,quote=F)


























#####Merge candidate gene and MR result

path1="/share/home/liujianjun/ESRD_Project/Table/Input/"
path2="/share/home/liujianjun/ESRD_Project/Table/Output/"
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"


library(data.table)
library(dplyr)

setwd(path1)
#MR results
ESRD=fread("ESRD_ivw_MR_Result.txt")
eGFR=fread("eGFR1_ivw_wald_Result.txt")
names(ESRD)[2]="Outcome"
names(eGFR)[2]="Outcome"
ESRD$Outcome="ESRD"
eGFR$Outcome="eGFR"
disease=rbind(ESRD,eGFR)
#result=subset(disease, pval < 0.05)
result$id.exposure=c()
result$outcome=c()
names(result)[2]="CpG"

#Prioritization gene
gene=fread("Final.ESRD_DMP_Gene.txt")

#Colocalized gene
meqtl_eqtl=fread("MeQTL_eQTL_all_tissue.Results.txt")
meqtl_eqtl=subset(meqtl_eqtl, tissue== "blood")
sub_meqtl_eqtl=subset(meqtl_eqtl, PP.H4.abf > 0.75)
sub_meqtl_eqtl <- subset(sub_meqtl_eqtl, select = -c(gene, nsnps, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf, block))
sub_meqtl_eqtl$SNP=c()
sub_meqtl_eqtl$tissue=c()
sub_meqtl_eqtl=unique(sub_meqtl_eqtl)
dim(sub_meqtl_eqtl)
sub_meqtl_eqtl$coloc=1



#gene colocalized with eGFR
eGFR=fread("Blood_Gene_eGFR_Coloc_Results.txt")
eGFR= subset(eGFR,PP.H4.abf > 0.75)
eqtl_eGFR2=eGFR[,c("SYMBOL")]
eqtl_eGFR2$eqtl_eGFR=1
names(eqtl_eGFR2)[1]="gene"

#eGFR TWAS gene
TWAS=fread("Full_Summary_eGFR_TWAS.txt")
TWAS$TWAS=1
names(TWAS)[3]="gene"



##合并MR 位点与prioritized gene

data1=merge(result,gene,by="CpG",all.x=T)

data1$eQTM_blood <- ifelse(rowSums(data1[, c("eqtm_2015", "eqtm_2018_GTP", "eqtm_2018_MESA", "eqtm_2023")]) >= 1, 1, 0)
data1$eQTM_kidney <- ifelse(rowSums(data1[, c("eqtm_2022_kidney")]) >= 1, 1, 0)
data1$Enhancer_linked_gene <- ifelse(rowSums(data1[, c("atac", "epimap", "hic","ABC")]) > 1, 1, 0)
data1$Nearest_gene <- ifelse(data1$epic_nearest > 1, 1, 0)


new_dt1=subset(data1, select = -c(atac, epic_nearest, epimap, hic, eqtm_2015, eqtm_2018_GTP, eqtm_2018_MESA, eqtm_2022_kidney, eqtm_2023, ABC, group, sum))




######进一步合并甲基化共定位基因
sub_meqtl_eqtl1=subset(sub_meqtl_eqtl, CpG %in% new_dt1$CpG)
 names(sub_meqtl_eqtl1)[2]="gene"
new_dt2=merge(new_dt1,sub_meqtl_eqtl1,by=c("CpG","gene"),all=T)
length(unique(new_dt2$gene))
thres=0.05/length(unique(new_dt2$gene))


####评估new_dt2的基因是否共定位到eGFR或者被报道为eGFR TWAS 基因
new_dt3=merge(new_dt2,eqtl_eGFR2,by="gene",all.x=T)
TWAS=subset(TWAS, twas_Pvalue <= thres)
new_dt4=merge(new_dt3,TWAS,by="gene",all.x=T)
new_dt4$Outcome=c()

unique(subset(new_dt4, select= c(gene,CpG,Nearest_gene,Enhancer_linked_gene,eQTM_blood ,eQTM_kidney ,coloc ,eqtl_eGFR,twas_Tarit))) #查看结果



#####合并ESRD EWAS结果
data=read.csv(paste0(index,"Full.Rep.DMP.Results"),sep=" ")

# Step 1: Remove columns containing "_meta"
data <- data[, !grepl("_meta", names(data))]

# Step 2: Combine columns starting with "Bacon_b_"
data$Beta <- apply(data[, grepl("Bacon_b_", names(data))], 1, function(x) paste(x, collapse = "; "))

# Step 3: Combine columns starting with "Bacon_P_"
data$Pvalue <- apply(data[, grepl("Bacon_P_", names(data))], 1, function(x) paste(x, collapse = "; "))

# Step 4: Remove original "Bacon_b_*" and "Bacon_P_*" columns
data <- data[, !grepl("Bacon_b_|Bacon_P_", names(data))]

# View the modified dataframe
data=subset(data,select=c(CpG,Beta,Pvalue ))   #######疾病顺序DM, SLEN, IgAN, PKD, HRD


##合并
all=merge(new_dt4,data,by="CpG",all.x=TRUE)
all[,c("CpG","b","Beta")]


write.csv(all,paste0(out,"Table.MR.result.txt"),row.names=F,quote=F)

####查看association 和 MR estimate方向是否一致
df=merge(result,data,by="CpG")
df$causal.estimate=df$b
df$Pvalue=c()
df$trait=df$Outcome
df$exposure=df$CpG















######Gene regulated by multiple DMPs based on gene prioritization results

out="/share/home/liujianjun/ESRD_Project/Table/Output/"
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
path1="/share/home/liujianjun/ESRD_Project/Table/Input/"

#Prioritization gene
setwd(path1)
gene=fread("Final.ESRD_DMP_Gene.txt")



data=data.frame(unlist(table(gene$gene)))
names(data)[1]="Gene"
target=subset(data,Freq >= 2)

sub_gene=subset(gene, gene %in% target$Gene)


write.csv(sub_gene,paste0(out,"Gene_regulated.by_Multiple.DMPs.csv"),row.names=F,quote=F)






















#####Table: meQTL and eGFR GWAS


#meqtl and eGFR
qtl="/share/home/liujianjun/ESRD_Project/DataBase/01.eQTLGen/"
obtain="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/"
con="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/Conditional.Results/"
col="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/meQTL_data/"
col_gwas="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/GWAS_data/"
out="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/Results/"
gtex_1="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/eQTL_data/GTEx/"



library(data.table)
library(dplyr)
library(coloc)
library(data.table)
library(locuscomparer) #####作图用

result=fread(paste0(col,"DMP_meQTL.txt"))  #meQTL
names(result)[7:8]=c("t_stat","p_value")
result_meqtl= result %>%  select(rsID, CHR,bp,beta,t_stat,p_value,MAF,NCHROBS,Group,block,gene)
result_meqtl$se=result_meqtl$beta/result_meqtl$t_stat
result_meqtl$N=result_meqtl$NCHROBS/2
result_meqtl$CpG=result_meqtl$gene
result_meqtl$gene=c()
result_meqtl$t_stat=c()
result_meqtl$NCHROBS=c()
names(result_meqtl)[1]="rs_id"
result_meqtl=unique(result_meqtl)


##去除meqtl数量小于100的block否则会报错
table=as.data.frame(table(result_meqtl$block))
ID=which(table$Freq <= 100)
 result_meqtl=subset(result_meqtl, ! block %in% table$Var1[ID])

lead=fread(paste0(con,"All.Lead.meQTL.SNP.txt"))#lead SNP
lead$block=paste0("block",1:nrow(lead))
lead=subset(lead,block %in% result_meqtl$block)

#给lead增加rs_id
 meqtl=result_meqtl[,c(1:3)]
 names(lead)[1]="CHR"
 lead_rsid=unique(merge(lead,meqtl,by=c("CHR","bp"),all.x=T))
 dim(lead_rsid)
 nrow(lead_rsid) == nrow(lead)


####GWAS: eGFR
gwas=fread(paste0(col_gwas,"Formatted.eGFR.GWAS.tsv"),fill=TRUE) #gwas

gwas_new=gwas %>% select ("rs_id","MAF","beta","se","p_value","N","A1","A2")
names(gwas_new)=c("rs_id","MAF","beta","se","p_value","N","effect_allele","other_allele")

gwas_new=unique(gwas_new)
if (length(which(gwas_new$beta == 0)) > 0) {gwas_new=gwas_new[-c(which(gwas_new$beta == 0)),]};
gwas_new=subset(gwas_new, rs_id %in% result_meqtl$rs_id)
gwas_new$N=406504


#sentienl meQTL SNP在meQTL中的统计量
lead_meqtl=subset(result_meqtl, rs_id %in% lead_rsid$rs_id)


#sentienl meQTL SNP在eGFR GWAS中的统计量
gwas_new_meqtl=subset(gwas_new, rs_id %in% lead_rsid$rs_id)



frame=c()
for ( i in unique(result_meqtl$block))  {
result_meqtl_sub= result_meqtl %>% filter(block == i)
all=merge( result_meqtl_sub,gwas_new, by="rs_id", all.x=T,suffixes=c("_meqtl","_gwas"))  ##连接方式一定要是all=TRUE 否则总是报错
all$block=i
all=unique(all)
   if ( length(which(all$MAF == 0)) > 0 ) {all=all[-which(all$MAF == 0),]}
input=all
dim(input)
 if (length(which(duplicated(input$rs_id))) > 0) {input=input[-(which(duplicated(input$rs_id))),]};
 input$N_gwas=406504;
if (length(which(is.na(input$p_value_gwas))) > 0) {input=input[-(which(is.na(input$p_value_gwas))),]};
if (length(which(is.na(input$p_value_beta))) > 0) {input=input[-(which(is.na(input$p_value_beta))),]};

 input$varbeta_meqtl=input$se_meqtl*input$se_meqtl;
  input$varbeta_gwas=input$se_gwas*input$se_gwas;


 # 获取 sentiel SNP
 sentienl_rs_id = lead_rsid$rs_id[which(lead_rsid$block == i)]

# coloc
result <- coloc.abf(
  dataset1 = list(
    pvalues = input$p_value_meqtl,
    snp = input$rs_id,
    beta = input$beta_meqtl,
    varbeta = input$varbeta_meqtl,
    type = "quant", 
    N = input$N_meqtl,
    MAF = input$MAF_meqtl
  ),
  dataset2 = list(
    pvalues = input$p_value_gwas,
    beta = input$beta_gwas,
    varbeta = input$varbeta_gwas,
    type = "quant",
    N = input$N_gwas,
    MAF = input$MAF_gwas,
    snp = input$rs_id
  )
)

  # 提取结果 summary，添加 gene, CpG, block, SNP 信息
      sum         = as.data.frame(t(result$summary))
      sum$Trait    = "eGFR"
      sum$CpG     = input$CpG[1]
      sum$block   = input$block[1]
      sum$SNP     = lead_rsid$SNP[which(lead_rsid$block == i)]
      sum$rs_id   = lead_rsid$rs_id[which(lead_rsid$block == i)]

 # 提取该 sentiel SNP 的统计量
      ID=which(all$rs_id == sum$rs_id)
            df =all[ID, c("rs_id", "CHR", "bp", "beta_meqtl", "p_value_meqtl", "MAF_meqtl", "block","beta_gwas", "se_gwas", "p_value_gwas")]  
            
             min_ID=which(input$p_value_gwas==min(input$p_value_gwas))  ##input就是过滤了NA的all
            df$p_value_gwas=min(input$p_value_gwas)
            df$beta_gwas=input$beta_gwas[min_ID]
            df$se_gwas=input$se_gwas[min_ID]

#将共定位结果和sentiel SNP 的统计量横向合并
      sum_all = cbind(sum, df)

##纵向合并所有coloc结果
      frame   = rbind(frame, sum_all)

}



frame_eGFR=subset(frame, (PP.H3.abf > 0.6 | PP.H4.abf > 0.6))   ###过滤条件





frame_eGFR$Group_backup=frame_eGFR$Group
frame_eGFR$Group=c()
frame_eGFR=unique(frame_eGFR)

n1=length(unique(frame$CpG))  #n1=28
n2=length(unique(frame_eGFR$CpG))  #n2=11
n1
n2

index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"

sentinels <- read.csv(paste0(index, "Full.Rep.DMP.Results"), sep = " ", header = TRUE)
sentinels=sentinels[,c(1,24:28)]

data_out=merge(frame_eGFR,sentinels,by="CpG")
data_out$dis=abs(data_out$MAPINFO - data_out$bp)
data_out$Group=ifelse( data_out$CHR.y == data_out$CHR.x, "cis_or_longcis", "trans")
data_out=subset(data_out, Group  == "cis_or_longcis")

table="/share/home/liujianjun/ESRD_Project/Table/Output/"
write.table(data_out,paste0(table,"Table.Coloc.CpG.eGFR.txt"),row.names=F,quote=F)

write.csv(data_out,paste0(table,"Table.Coloc.CpG.eGFR.csv"),row.names=F,quote=F)








































































#####Table: meQTL and BUN GWAS


#meqtl and BUN
qtl="/share/home/liujianjun/ESRD_Project/DataBase/01.eQTLGen/"
obtain="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/"
con="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/Conditional.Results/"
col="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/meQTL_data/"
col_gwas="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/GWAS_data/"
out="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/Results/"
gtex_1="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/eQTL_data/GTEx/"
path_1="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Causal_Effect/Published.GWAS.Datasets/"




library(data.table)
library(dplyr)
library(coloc)
library(data.table)
library(locuscomparer) #####作图用

result=fread(paste0(col,"DMP_meQTL.txt"))  #meQTL
names(result)[7:8]=c("t_stat","p_value")
result_meqtl= result %>%  select(rsID, CHR,bp,beta,t_stat,p_value,MAF,NCHROBS,Group,block,gene)
result_meqtl$se=result_meqtl$beta/result_meqtl$t_stat
result_meqtl$N=result_meqtl$NCHROBS/2
result_meqtl$CpG=result_meqtl$gene
result_meqtl$gene=c()
result_meqtl$t_stat=c()
result_meqtl$NCHROBS=c()
names(result_meqtl)[1]="rs_id"
result_meqtl=unique(result_meqtl)


##去除meqtl数量小于100的block否则会报错
table=as.data.frame(table(result_meqtl$block))
ID=which(table$Freq <= 100)
 result_meqtl=subset(result_meqtl, ! block %in% table$Var1[ID])

lead=fread(paste0(con,"All.Lead.meQTL.SNP.txt"))#lead SNP
lead$block=paste0("block",1:nrow(lead))
lead=subset(lead,block %in% result_meqtl$block)

#给lead增加rs_id
 meqtl=result_meqtl[,c(1:3)]
 names(lead)[1]="CHR"
 lead_rsid=unique(merge(lead,meqtl,by=c("CHR","bp"),all.x=T))
 dim(lead_rsid)
 nrow(lead_rsid) == nrow(lead)


####GWAS: BUN
gwas=fread(paste0(col_gwas,"BUN.mean.GWAS.tsv"),fill=TRUE) #gwas

gwas_new=gwas %>% select ("rsid","effect_allele_frequency","beta","standard_error","p_value","n","effect_allele","other_allele")

names(gwas_new)=c("rs_id","MAF","beta","se","p_value","N","effect_allele","other_allele")

gwas_new=unique(gwas_new)
if (length(which(gwas_new$beta == 0)) > 0) {gwas_new=gwas_new[-c(which(gwas_new$beta == 0)),]};
gwas_new=subset(gwas_new, rs_id %in% result_meqtl$rs_id)



##full meQTL 和full BUN GWAS共有的SNP

gwas_new=subset(gwas_new, rs_id %in% result_meqtl$rs_id)

#sentienl meQTL SNP在meQTL中的统计量
lead_meqtl=subset(result_meqtl, rs_id %in% lead_rsid$rs_id)


#sentienl meQTL SNP在BUN GWAS中的统计量
gwas_new_meqtl=subset(gwas_new, rs_id %in% lead_rsid$rs_id)



frame=c()
for ( i in unique(result_meqtl$block))  {
result_meqtl_sub= result_meqtl %>% filter(block == i)
all=merge( result_meqtl_sub,gwas_new, by="rs_id", all.x=T,suffixes=c("_meqtl","_gwas"))  ##连接方式一定要是all.x=TRUE 否则总是报错
all$block=i
all=unique(all)
   if ( length(which(all$MAF == 0)) > 0 ) {all=all[-which(all$MAF == 0),]}
input=all
dim(input)
 if (length(which(duplicated(input$rs_id))) > 0) {input=input[-(which(duplicated(input$rs_id))),]};
if (length(which(is.na(input$p_value_gwas))) > 0) {input=input[-(which(is.na(input$p_value_gwas))),]};
if (length(which(is.na(input$p_value_beta))) > 0) {input=input[-(which(is.na(input$p_value_beta))),]};

 input$varbeta_meqtl=input$se_meqtl*input$se_meqtl;
  input$varbeta_gwas=input$se_gwas*input$se_gwas;


 # 获取 sentiel SNP
 sentienl_rs_id = lead_rsid$rs_id[which(lead_rsid$block == i)]

# coloc
result <- coloc.abf(
  dataset1 = list(
    pvalues = input$p_value_meqtl,
    snp = input$rs_id,
    beta = input$beta_meqtl,
    varbeta = input$varbeta_meqtl,
    type = "quant", 
    N = input$N_meqtl,
    MAF = input$MAF_meqtl
  ),
  dataset2 = list(
    pvalues = input$p_value_gwas,
    beta = input$beta_gwas,
    varbeta = input$varbeta_gwas,
    type = "quant",
    N = input$N_gwas,
    MAF = input$MAF_gwas,
    snp = input$rs_id
  )
)

      # 提取结果 summary，添加 gene, CpG, block, SNP 信息
      sum         = as.data.frame(t(result$summary))
      sum$Trait    = "BUN"
      sum$CpG     = input$CpG[1]
      sum$block   = input$block[1]
      sum$SNP     = lead_rsid$SNP[which(lead_rsid$block == i)]
      sum$rs_id   = lead_rsid$rs_id[which(lead_rsid$block == i)]

 # 提取该 sentiel SNP 的统计量
      ID=which(all$rs_id == sum$rs_id)
            df =all[ID, c("rs_id", "CHR", "bp", "beta_meqtl", "p_value_meqtl", "MAF_meqtl", "block","effect_allele","other_allele")]  
            
            min_ID=which(input$p_value_gwas==min(input$p_value_gwas))  ##input就是过滤了NA的all
            df$p_value_gwas=min(input$p_value_gwas)
            df$beta_gwas=input$beta_gwas[min_ID]
            df$se_gwas=input$se_gwas[min_ID]
          


#将共定位结果和sentiel SNP 的统计量横向合并
      sum_all = cbind(sum, df)

##纵向合并所有coloc结果
      frame   = rbind(frame, sum_all)

}



frame_BUN=subset(frame, (PP.H3.abf > 0.6 | PP.H4.abf > 0.6) )   ###过滤条件


n1=length(unique(frame$CpG))  #n1=28
n2=length(unique(frame_BUN$CpG))  #n2=10

frame_BUN$Group=c()
frame_BUN=unique(frame_BUN)



index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"

sentinels <- read.csv(paste0(index, "Full.Rep.DMP.Results"), sep = " ", header = TRUE)
sentinels=sentinels[,c(1,24:28)]

data_out=merge(frame_BUN,sentinels,by="CpG")
data_out$dis=abs(data_out$MAPINFO - data_out$bp)

data_out$Group=ifelse( data_out$CHR.y == data_out$CHR.x, "cis_or_longcis", "trans")

table="/share/home/liujianjun/ESRD_Project/Table/Output/"


write.csv(data_out,paste0(table,"Table.Coloc.CpG.BUN.csv"),row.names=F,quote=F)

























####Table: colocalization: cis-meQTL and cis-eQTL  (参考文献https://pmc.ncbi.nlm.nih.gov/articles/PMC7617265/； 文章中关键词搜索：“Coloc analysis was subsequently performed for loci with a potentially causal relationship between DNA methylation levels and gene expression in cis ”)

#meqtl
qtl="/share/home/liujianjun/ESRD_Project/DataBase/01.eQTLGen/"
obtain="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/"
con="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/Conditional.Results/"
col="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/meQTL_data/"
col_gwas="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/GWAS_data/"
out="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/Results/"
gtex_1="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/eQTL_data/GTEx/"



library(data.table)
library(dplyr)
library(coloc)
library(data.table)
library(locuscomparer) #####作图用

result=fread(paste0(col,"DMP_meQTL.txt"))  #meQTL
names(result)[7:8]=c("t_stat","p_value")
result_meqtl= result %>%  select(rsID, CHR,bp,beta,t_stat,p_value,MAF,NCHROBS,Group,block,gene)
result_meqtl$se=result_meqtl$beta/result_meqtl$t_stat
result_meqtl$N=result_meqtl$NCHROBS/2
result_meqtl$CpG=result_meqtl$gene
result_meqtl$gene=c()
result_meqtl$t_stat=c()
result_meqtl$NCHROBS=c()
names(result_meqtl)[1]="rs_id"
result_meqtl=unique(result_meqtl)


##去除meqtl数量小于100的block否则会报错
table=as.data.frame(table(result_meqtl$block))
ID=which(table$Freq <= 100)
 result_meqtl=subset(result_meqtl, ! block %in% table$Var1[ID])

lead=fread(paste0(con,"All.Lead.meQTL.SNP.txt"))#lead SNP
lead$block=paste0("block",1:nrow(lead))
lead=subset(lead,block %in% result_meqtl$block)


#给lead增加rs_id
 meqtl=result_meqtl[,c(1:3)]
 names(lead)[1]="CHR"
 lead_rsid=unique(merge(lead,meqtl,by=c("CHR","bp"),all.x=T))
 dim(lead_rsid)
 nrow(lead_rsid) == nrow(lead)


 #eqtl


eqtl=fread(paste0(qtl,"Final.blood.eQTL.txt")) #eQTL
names(eqtl)
eqtl_new=eqtl %>% select ( SNP,Freq,Gene,b,SE,p,N)
names(eqtl_new)=c("rs_id","MAF","Gene","beta","se","p_value","N")
eqtl_new=unique(eqtl_new)
eqtl_new=subset(eqtl_new,rs_id %in% result_meqtl$rs_id)


#Coloc : 
#注意： 我参考了John Chamber那篇文章，文章的补充材料表格主要是在有共定位证据（PPH3 or PPH4  > 0.8) 的基础上分别提取了，each lead meQTL SNP of ESRD DMP在meQTL和eQTL中的summary,然后汇总到表格中。
#注意： 有些loci虽然显示了共定位的证据，但是lead meQTL SNP of ESRD DMP在eqtl中并不显著，因此我们增加了筛选条件，也就是必须满足该lead meQTL SNP在eQTL中必须满足P<5e-4（suggestive significant  threshold)
final = c()

for (i in unique(result_meqtl$block)) {
  
  result_meqtl_sub = result_meqtl %>% filter(block == i)
  
  all = merge(
    result_meqtl_sub, 
    eqtl_new, 
    by = "rs_id", 
    all.x = TRUE,
    suffixes = c("_meqtl", "_eqtl")
  )
  
  all$block = i
  all = unique(all)
  
  # 去除 MAF 为 0 的数据
  if (length(which(all$MAF_eqtl == 0)) > 0) {
    all = all[-which(all$MAF_eqtl == 0), ]
  }

  # 去除缺失 MAF 的行
  if (length(which(is.na(all$MAF_eqtl)))) {
    all = all[-which(is.na(all$MAF_eqtl)), ]
  }

  gene = unique(all$Gene)
  frame = c()

  for (j in gene) {

    input = subset(all, Gene == j)

    # 添加方差列
    input$varbeta_meqtl = input$se_meqtl * input$se_meqtl
    input$varbeta_eqtl = input$se_eqtl * input$se_eqtl

    # 去除重复 SNP
    if (length(which(duplicated(input$rs_id))) > 0) {
      input = input[-which(duplicated(input$rs_id)), ]
    }

    # 获取 sentiel SNP
    sentienl_rs_id = lead_rsid$rs_id[which(lead_rsid$block == i)]

    # 如果 sentienl SNP 存在于该基因的 eQTL 数据中
    if (any(input$rs_id %in% sentienl_rs_id)) {

      result <- coloc.abf(
        dataset1 = list(
          pvalues = input$p_value_meqtl,
          snp     = input$rs_id,
          beta    = input$beta_meqtl,
          varbeta = input$varbeta_meqtl,
          type    = "quant", 
          N       = input$N_meqtl,
          MAF     = input$MAF_meqtl
        ),
        dataset2 = list(
          pvalues = input$p_value_eqtl,
          snp     = input$rs_id,
          beta    = input$beta_eqtl,
          varbeta = input$varbeta_eqtl,
          type    = "quant",
          N       = input$N_eqtl,
          MAF     = input$MAF_eqtl
        )
      )

      # 提取结果 summary，添加 gene, CpG, block, SNP 信息
      sum         = as.data.frame(t(result$summary))
      sum$gene    = j
      sum$CpG     = input$CpG[1]
      sum$block   = input$block[1]
      sum$SNP     = lead_rsid$SNP[which(lead_rsid$block == i)]
      sum$rs_id   = lead_rsid$rs_id[which(lead_rsid$block == i)]

      # 提取该 sentiel SNP 的统计量
      ID = which(input$rs_id == sum$rs_id)
      df = input[ID, c("rs_id", "CHR", "bp", "beta_meqtl", "p_value_meqtl", "MAF_meqtl", "Group",
      "Gene", "beta_eqtl", "se_eqtl", "p_value_eqtl")]  

      sum_all = cbind(sum, df)
      frame   = rbind(frame, sum_all)
    }
  }

  final = rbind(final, frame)
}

final_gene=subset(final, (PP.H3.abf > 0.8 | PP.H4.abf > 0.8) & p_value_eqtl < 0.0005)   ###过滤条件


##将该结果和eGFR整合：目的：上述基因final_gene与甲基化共定位，下面的步骤是为了看看其中哪些基因和eGFR共定位

#该结果是eQTLGEN eQTL和eGFR GWAS共定位的结果 （full)
out="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/Results/"  
eGFR_col=fread(paste0(out,"Gene_eGFR_Coloc_Results.txt"))
eGFR_col=eGFR_col[,c("Gene","PP.H4.abf","trait")]
names(eGFR_col)[2]="PP.H4.eGFR"

#合并
meth_gene_eGFR=merge(final_gene,eGFR_col,by="Gene",all.x=T)

meth_gene_eGFR$gene_cloc_eGFR=ifelse(meth_gene_eGFR$PP.H4.eGFR > 0.7,"*","-")


sub=meth_gene_eGFR

#转化基因为symbol
 library(ChIPseeker)
library(org.Hs.eg.db)
library(GenomicFeatures)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library('clusterProfiler')


library(dplyr)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

all=all %>%  mutate(Gene = gsub("\\.\\d+$", "", Gene)) 

gene<- bitr(sub$gene, fromType = "ENSEMBL", toType=c("SYMBOL"),OrgDb = org.Hs.eg.db)
names(gene)[1]="gene"


all=merge(sub,gene,by="gene")

#write.table(all,paste0(out,"Methylation_Gene_Coloc_ResultsV1.txt"),row.names=F,quote=F)


index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"

sentinels <- read.csv(paste0(index, "Full.Rep.DMP.Results"), sep = " ", header = TRUE)
sentinels=sentinels[,c(1,24:28)]

data_out=merge(all,sentinels,by="CpG")



###加上基因座标
path="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/Results/"

###这个bulid 是在metascape网站上下载的
build=fread(paste0(path,"Coloc.eQTLGEN.Gene.Build"))
df=build
names(df)=c("SYMBOL","Gene_Location")

library(tidyr)

# Split Gene_Location into Chromosome, Start, and End
df_split <- df %>%
  separate(Gene_Location, 
           into = c("Chromosome", "Position"), 
           sep = ":", 
           remove = FALSE) %>%
  separate(Position, 
           into = c("Start", "End"), 
           sep = "-")

# View the result
print(df_split)

data_out1=merge(data_out,df_split,by="SYMBOL",all.x=T)

####判断eQTL是否为cis
data_out1$Start=as.numeric(data_out1$Start)
data_out1$eQTL.distance=abs(data_out1$bp - data_out1$Start)
data_out1$eQTL.Group= ifelse(data_out1$eQTL.distance > 50000 ,  "trans","cis")


table="/share/home/liujianjun/ESRD_Project/Table/Output/"
write.table(data_out1,paste0(table,"Table.Methylation_Gene_Coloc_Results.txt"),row.names=F,quote=F)
write.csv(data_out1,paste0(table,"Table.Methylation_Gene_Coloc_Results.csv"),row.names=F,quote=F)































##########Main Table 1

table="/share/home/liujianjun/ESRD_Project/Table/Output/"
out="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Result/"
mr_result="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Results/"


library(data.table)
library(dplyr)


coloc.gene=fread(paste0(table,"Table.Methylation_Gene_Coloc_Results.csv"),sep=",") #共定位基因

prio.gene=fread(paste0(out,"Final.ESRD_DMP_Gene.txt"),sep="\t")  #优选基因


eGFR=fread(paste0(table,"Table.Overlapped_eGFR_DMP.txt"),sep="\t")  #look up in eGFR EWAS的结果
eGFR0 <- eGFR %>% select(CpG,  beta_eGFR, Pval_eGFR, PMID_eGFR)

ESRD=fread(paste0(table,"eGFR_overlap.DMP_ESRD_association.txt"),sep="\t") #DMP在多个ESRD Group中的summary statistics

 
target.DMP=fread(paste0(table,"Target.list.txt"))    #与eGFR：既有lookup证据又有共定位或者MR证据的DMP位点




#MR result

BUN=fread(paste0(mr_result,"BUN_All_MR_Result.txt"))
eGFR1=fread(paste0(mr_result,"eGFR1_All_MR_Result.txt"))


sig.BUN=subset(BUN,pval < 0.05 & dir == "yes")
sig.eGFR=subset(eGFR1,pval < 0.05 & dir == "no")

all_mr=rbind(sig.BUN,sig.eGFR)

all_mr1 <- all_mr %>%
  select(CpG, outcome, nsnp, b, se, pval) %>%
  rename(mr_b = b, mr_se = se, mr_pval = pval)



##Coloc Result

coloc.eGFR=fread(paste0(table,"Table.Coloc.CpG.eGFR.txt"),fill=TRUE)
coloc.eGFR1=coloc.eGFR %>%  select (CpG,nsnps,PP.H3.abf,PP.H4.abf,rs_id,p_value_meqtl,Trait)


coloc.BUN=fread(paste0(table,"Table.Coloc.CpG.BUN.csv"),fill=TRUE)
coloc.BUN1=coloc.BUN %>%  select (CpG,nsnps,PP.H3.abf,PP.H4.abf,rs_id,p_value_meqtl,Trait)

coloc1=rbind(coloc.eGFR1,coloc.BUN1)

coloc1=subset(coloc1,PP.H3.abf > 0.6 | PP.H4.abf > 0.6)



###以下步骤是获取所有的基因：包括共定位和优选基因
head(coloc.gene)
head(prio.gene)

coloc.gene1=coloc.gene[,c(2,1)]
names(coloc.gene1)=c("CpG","Gene")

prio.gene1=prio.gene[,c(1,2)]
names(prio.gene1)=c("CpG","Gene")


Gene=rbind(coloc.gene1,prio.gene1)

sub=unique(subset(Gene, CpG %in% target.DMP$CpG))

setDT(sub)

# Aggregate by CpG, concatenate genes
result<- sub[, .(Gene = toString(unique(Gene))), by = CpG]





##合并所有结果

all <- result %>%
  left_join(ESRD, by = "CpG") %>%
  left_join(eGFR0, by = "CpG") %>%
  left_join(all_mr1, by = "CpG") 

 all_result= merge(all,coloc1,by="CpG",all.x=T)



write.table(all_result,paste0(table,"Main.Table1.Association.kidney.function.txt"),sep="/")















































#########################All significant meQTL SNPs: add rsID

#meqtl
qtl="/share/home/liujianjun/ESRD_Project/DataBase/01.eQTLGen/"
obtain="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/"
con="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/Conditional.Results/"
col="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/meQTL_data/"
col_gwas="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/GWAS_data/"
out="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/Results/"
gtex_1="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/eQTL_data/GTEx/"
table="/share/home/liujianjun/ESRD_Project/Table/Output/"



library(data.table)
library(dplyr)
library(coloc)
library(data.table)
library(locuscomparer) #####作图用

result=fread(paste0(col,"DMP_meQTL.txt"))  #meQTL
names(result)[7:8]=c("t_stat","p_value")
result_meqtl= result %>%  select(rsID, CHR,bp,beta,t_stat,p_value,MAF,NCHROBS,Group,block,gene)
library(data.table)


raw=fread("/share/home/liujianjun/ESRD_Project/tmp/meQTL.tmp.txt")  ###手动准备：这是补充材料表的原始版本，将其复制粘贴到集群上，需要手动准备

names(raw)[c(11,13)]=c("chr","bp")
result0=result[,c(1,2,16)]
result0=unique(result0)
final=merge(raw,result0,by=c("chr","bp"),all.x=T)


final$chr=paste0("chr",final$chr)
head(final)
setDT(final)
final0=final[, Group := fifelse(CpG.chr != chr, 
                         "trans", 
                         fifelse(abs(CpG_position - bp) > 1e6, 
                               "long cis", 
                                 "cis"))]

write.table(final0, paste0(table, "Table.ALL.significant.meQTL.txt"), quote = FALSE, sep = "\t", row.names = FALSE)















##############Lead SNP

qtl="/share/home/liujianjun/ESRD_Project/DataBase/01.eQTLGen/"
obtain="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/"
con="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/Conditional.Results/"
col="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/meQTL_data/"
col_gwas="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/GWAS_data/"
out="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/Results/"
gtex_1="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/eQTL_data/GTEx/"
table="/share/home/liujianjun/ESRD_Project/Table/Output/"



###手动准备：将补充材料表的“Table 13. lead meQTL SNP”的结果复制粘贴到/share/home/liujianjun/ESRD_Project/tmp/lead.meQTL.tmp.txt

raw=fread("/share/home/liujianjun/ESRD_Project/tmp/lead.meQTL.tmp.txt")


all=fread(paste0(table,"Table.ALL.significant.meQTL.txt"))

library(data.table)
library(dplyr)
library(coloc)
library(data.table)
library(locuscomparer) #####作图用

result=fread(paste0(col,"DMP_meQTL.txt"))  #meQTL
names(result)[7:8]=c("t_stat","p_value")
result_meqtl= result %>%  select(rsID, CHR,bp,beta,t_stat,p_value,MAF,NCHROBS,Group,block,gene)
library(data.table)




names(raw)[c(1,3)]=c("chr","bp")
result0=result[,c(1,2,16)]
result0=unique(result0)
final=merge(raw,result0,by=c("chr","bp"),all.x=T)

sub_lead=subset(all, rsID %in% final$rsID  & CpG %in% final$CpG)

write.table(sub_lead,paste0(table,"Table.ALL.Top.meQTL.txt"),quote = FALSE, sep = "\t", row.names = FALSE)































########################################Table: Genetic intruments


mr_ex="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Exposure_meQTL/"
mr_out="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Outcome_GWAS/"
mr_result="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Results/"
mr_code="/share/home/liujianjun/ESRD_Project/Part12.Causal_MR/Code/"
wd1="/share/home/liujianjun/GWAS.OF.ESRD/Important.file/"
col="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/meQTL_data/"
table="/share/home/liujianjun/ESRD_Project/Table/Output/"



###################################Harmonized data: eGFR GWAS summary were downloaded from https://www.nature.com/articles/s41467-024-53516-7#data-availability
library(ieugwasr)
library(TwoSampleMR)
library(data.table)
library(dplyr)

path_1="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Causal_Effect/Published.GWAS.Datasets/"

final=fread(paste0(mr_ex,"Pruned.cis-meQTL.exposure.txt")) #exposure



########meQTL SNP: BP CHR
location=fread(paste0(mr_ex,"cis-meQTL.rsID.Results.txt"))


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

dat_eGFR1=subset(Fdata, Fvalue >=10)
names(dat_eGFR1)[1]="rsid"

instrument_eGFR=merge(dat_eGFR1,location,by="rsid",all.x=T)


index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"

sentinels <- read.csv(paste0(index, "Full.Rep.DMP.Results"), sep = " ", header = TRUE)
sentinels=sentinels[,c(1,24:28)]

instrument_eGFR1=merge(instrument_eGFR,sentinels,by="CpG",all.x=T)


write.table(instrument_eGFR1,paste0(table,"Table.eGFR.Instruments.txt"),sep="\t",row.names=FALSE,quote=FALSE)







###########F statistics:BUN


#BUN
out_BUN<- read_outcome_data(filename=paste0(path_1,"BUN.mean.GWAS.tsv"),samplesize_col="n",snps=NULL,sep = "\t",snp_col ="rsid",beta_col = "beta", se_col = "standard_error",effect_allele_col = "effect_allele",other_allele_col = "other_allele", pval_col = "p_value",eaf_col="effect_allele_frequency")
out_BUN$outcome="BUN"

dat_BUN<- harmonise_data(exposure_dat = final,outcome_dat = out_BUN)
dat_BUN=subset(dat_BUN,mr_keep == "TRUE")




#F-statistics
dat_BUN$samplesize.exposure=923
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
names(dat_BUN)[1]="rsid"

instrument_BUN=merge(dat_BUN,location,by="rsid",all.x=T)



index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"

sentinels <- read.csv(paste0(index, "Full.Rep.DMP.Results"), sep = " ", header = TRUE)
sentinels=sentinels[,c(1,24:28)]

instrument_BUN1=merge(instrument_BUN,sentinels,by="CpG",all.x=T)

write.table(instrument_BUN1,paste0(table,"Table.BUN.Instruments.txt"),sep="\t",row.names=FALSE,quote=FALSE)






















############################################################Gene 与eGFR和BUN的关系 (cis-meQTL- cis-eQTL MR)



###############Step1. 选取有强烈证据等级的基因，进行MR和coloc分析

library(data.table)

tmp="/share/home/liujianjun/ESRD_Project/tmp/"
table="/share/home/liujianjun/ESRD_Project/Table/Output/"
input="/share/home/liujianjun/ESRD_Project/Table/Input/"


list=fread(paste0(tmp,"Gene.list.txt"))

#将SYMBOL 转化为 ENSG
 library(ChIPseeker)
library(org.Hs.eg.db)
library(GenomicFeatures)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library('clusterProfiler')


library(dplyr)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene


gene<- bitr(list$SYMBOL, fromType = "SYMBOL", toType=c("ENSEMBL"),OrgDb = org.Hs.eg.db)

all=merge(list,gene,by="SYMBOL")


write.table(all,paste0(input,"Target.Gene.list.txt"),row.names=F,quote=F,sep="\t")


######load eQTL data


qtl="/share/home/liujianjun/ESRD_Project/DataBase/01.eQTLGen/"
obtain="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/"
con="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/Conditional.Results/"
col="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/meQTL_data/"
col_gwas="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/GWAS_data/"
out="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/Results/"
gtex_1="/share/home/liujianjun/ESRD_Project/Part14.Colcalization/00.Input_files/eQTL_data/GTEx/"


eqtl=fread(paste0(qtl,"Final.blood.eQTL.txt")) #eQTL
names(eqtl)
eqtl_new=eqtl %>% select ( SNP,Freq,Gene,b,SE,p,N)
names(eqtl_new)=c("rs_id","MAF","Gene","beta","se","p_value","N")
eqtl_new=unique(eqtl_new)
eqtl_new=subset(eqtl_new,rs_id %in% result_meqtl$rs_id)




#####extract eQTL of the target gene


#!/usr/bin/perl
use strict;
use warnings;

# Define file paths
my $file1_path = '/share/home/liujianjun/ESRD_Project/DataBase/01.eQTLGen/Final.blood.eQTL.txt';
my $file2_path = '/share/home/liujianjun/ESRD_Project/Table/Input/Target.Gene.list.txt';
my $output_file = '/share/home/liujianjun/ESRD_Project/Table/Input/Target.Gene.eQTL.txt';  # Output file path

# Open both files for reading
open my $file1, '<', $file1_path or die "Cannot open $file1_path: $!";
open my $file2, '<', $file2_path or die "Cannot open $file2_path: $!";
open my $output_fh, '>', $output_file or die "Cannot create $output_file: $!";  # Open output file for writing

# Read the ENSEMBL IDs from the second file into a hash for fast lookup
my %ensembl_ids;
while (my $line = <$file2>) {
    chomp $line;
    next if $. == 1;  # skip header if present
    my ($symbol, $ensembl) = split /\s+/, $line;
    $ensembl_ids{$ensembl} = 1;
}
close $file2;

# Print header to output file
my $header = <$file1>;
print $output_fh $header;

# Process the first file and print matching rows to output file
while (my $line = <$file1>) {
    chomp $line;
    my @fields = split /\s+/, $line;
    my $gene = $fields[9];  # Assuming gene column is in the 10th position (index 9)
    if (exists $ensembl_ids{$gene}) {
        print $output_fh "$line\n";
    }
}

close $file1;
close $output_fh;

print "Matching rows have been written to $output_file.\n";










############Step2. Gene 与eGFR和BUN的关系 (cis-meQTL- cis-eQTL MR)




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
input="/share/home/liujianjun/ESRD_Project/Table/Input/"
path_1="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Causal_Effect/Published.GWAS.Datasets/"





library(data.table)
library(ieugwasr)
library(TwoSampleMR)


##load eQTL :exposure data

exp <- read_exposure_data(filename=paste0(input,"Target.Gene.eQTL.txt"),sep = "\t",snp_col ="SNP",beta_col = "b", se_col = "SE",effect_allele_col = "A1",other_allele_col = "A2", pval_col = "p",phenotype_col="Gene",eaf="Freq")
exp=subset(exp, pval.exposure < 5e-4)
exp=subset(exp,mr_keep.exposure=="TRUE")

a=dplyr::tibble(rsid=exp$SNP, pval=exp$pval.exposure, id=exp$exposure)


##########Pruning (HPC)

list=unique(a$id)
result=c()
for ( i in 1:length(list)) { data=subset(a,id==list[i]);data_clump=ld_clump(dat=data,clump_kb = 1000,clump_r2 = 0.2,plink_bin = "/share/home/liujianjun/software/plink",bfile =paste0(panel,"EUR"));result=rbind(result,data_clump)}

#可尝试将clump_r2 = 0.2 改为clump_r2 = 0.4


names(result)=c("SNP","pval","exposure")
LD=result[,c(1,3)]
final=merge(LD,exp,by=c("SNP","exposure"))
#save
write.table(final,paste0(mr_ex,"Pruned.cis-eQTL.exposure.txt"),row.names=FALSE,quote=FALSE)




#####Load GWAS data


#ESRD

out<- read_outcome_data(filename=paste0(mr_out,"ESRD.GWAS.rsID.Results.txt"),snps=NULL,sep = " ",snp_col ="rsid",beta_col = "b", se_col = "se",effect_allele_col = "A1",other_allele_col = "A2",eaf_col = "MAF", pval_col = "P")
out$outcome="ESRD"

final=fread(paste0(mr_ex,"Pruned.cis-eQTL.exposure.txt"))


dat <- harmonise_data(exposure_dat = final,outcome_dat = out)

#write.table(dat,paste0(mr_out,"cis-eQTL-ESRD-Harmonised.Data.txt"),row.names=FALSE,quote=FALSE)


######MR:eQTL and ESRD GWAS
hm=subset(dat ,mr_keep.outcome == "TRUE")


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


all_ESRD=rbind(ESRD_ivw,ESRD_wald)


###ENSEMBL to SYMBOL

#转化基因为symbol
 library(ChIPseeker)
library(org.Hs.eg.db)
library(GenomicFeatures)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library('clusterProfiler')


library(dplyr)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

# 假设 all_ESRD$exposure 包含 ENSEMBL ID
ensembl_ids <- all_ESRD$exposure

# 使用 biomaRt 进行映射
gene <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)


names(gene)[1]="exposure"


all=merge(all_ESRD,gene,by="exposure")

write.table(all,paste0(mr_result,"Gene_ESRD_MR_Result.txt"),quote=FALSE,row.names=FALSE,sep="\t")






#eGFR Gabriel 

out_eGFR1<- read_outcome_data(filename=paste0(path_1,"UCSF.Gabriel.eGFR.GWAS.tsv"),snps=NULL,sep = "\t",snp_col ="rsid",beta_col="beta", se_col = "standard_error",effect_allele_col = "effect_allele",other_allele_col = "other_allele",eaf_col="effect_allele_frequency", pval_col = "p_value")
out_eGFR1$outcome="eGFR1"


final=fread(paste0(mr_ex,"Pruned.cis-eQTL.exposure.txt"))

dat_eGFR1 <- harmonise_data(exposure_dat = final,outcome_dat = out_eGFR1)
 dat_eGFR1=subset(dat_eGFR1,mr_keep == "TRUE")


eGFR1_ivw <- mr(dat_eGFR1 , method_list = "mr_ivw")  #IVW


###ENSEMBL to SYMBOL

#转化基因为symbol
 library(ChIPseeker)
library(org.Hs.eg.db)
library(GenomicFeatures)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library('clusterProfiler')


library(dplyr)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

# 假设 all_ESRD$exposure 包含 ENSEMBL ID
ensembl_ids <- eGFR1_ivw$exposure

# 使用 biomaRt 进行映射
gene <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)



names(gene)[1]="exposure"


all=merge(eGFR1_ivw,gene,by="exposure")

write.table(all,paste0(mr_result,"Gene_eGFR_MR_Result.txt"),quote=FALSE,row.names=FALSE,sep="\t")










 
#BUN
out_BUN<- read_outcome_data(filename=paste0(path_1,"BUN.mean.GWAS.tsv"),samplesize_col="n",snps=NULL,sep = "\t",snp_col ="rsid",beta_col = "beta", se_col = "standard_error",effect_allele_col = "effect_allele",other_allele_col = "other_allele", pval_col = "p_value",eaf_col="effect_allele_frequency")
out_BUN$outcome="BUN"


final=fread(paste0(mr_ex,"Pruned.cis-eQTL.exposure.txt"))

dat_BUN<- harmonise_data(exposure_dat = final,outcome_dat = out_BUN)
dat_BUN=subset(dat_BUN,mr_keep == "TRUE")




BUN_ivw <- mr(dat_BUN , method_list = "mr_ivw")  #IVW


###ENSEMBL to SYMBOL

#转化基因为symbol
 library(ChIPseeker)
library(org.Hs.eg.db)
library(GenomicFeatures)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library('clusterProfiler')


library(dplyr)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

# 假设 all_ESRD$exposure 包含 ENSEMBL ID
ensembl_ids <- BUN_ivw$exposure

# 使用 biomaRt 进行映射
gene <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)


names(gene)[1]="exposure"


all=merge(BUN_ivw,gene,by="exposure")

write.table(all,paste0(mr_result,"Gene_BUN_MR_Result.txt"),quote=FALSE,row.names=FALSE,sep="\t")































#######MR: DNA methylation and gene expression



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
input="/share/home/liujianjun/ESRD_Project/Table/Input/"
path_1="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Causal_Effect/Published.GWAS.Datasets/"





library(data.table)
library(ieugwasr)
library(TwoSampleMR)



library(data.table)

#DNA methylation

final=fread(paste0(mr_ex,"Pruned.cis-meQTL.exposure.txt"))


#Gene expression
exp <- read_outcome_data(filename=paste0(input,"Target.Gene.eQTL.txt"),sep = "\t",snp_col ="SNP",beta_col = "b", se_col = "SE",effect_allele_col = "A1",other_allele_col = "A2", pval_col = "p",phenotype_col="Gene",eaf="Freq")



dat_eqtl <- harmonise_data(exposure_dat = final,outcome_dat = exp)
 dat_eqtl=subset(dat_eqtl,mr_keep == "TRUE")


eqtl_ivw <- mr(dat_eqtl , method_list = "mr_ivw")  #IVW




##Wald ratio

freq_table <- table(dat_eqtl$exposure) 
data<- data.frame(
  CpG = names(freq_table),
  Freq = as.numeric(freq_table)
)

sub=subset(data,Freq ==1)  #工具变量只有1个的CpG位点

list=unique(dat_eqtl$exposure)  ##Wald ratio
file=c()
for ( i in 1:length(list)) { data=subset(dat_eqtl,exposure ==list[i]);new=data[which(data$pval.exposure==min(data$pval.exposure)),];file=as.data.frame(rbind(file,new))}
file=subset(file, exposure %in% sub$CpG)
eqtl_wald<- mr(file , method_list = "mr_wald_ratio")


all_eqtl=rbind(eqtl_ivw,eqtl_wald)




###ENSEMBL to SYMBOL

#转化基因为symbol
 library(ChIPseeker)
library(org.Hs.eg.db)
library(GenomicFeatures)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library('clusterProfiler')


library(dplyr)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

# 假设 all_ESRD$exposure 包含 ENSEMBL ID
ensembl_ids <- all_eqtl$outcome

# 使用 biomaRt 进行映射
gene <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)


names(gene)[1]="outcome"


all=merge(all_eqtl,gene,by="outcome")

write.table(all,paste0(mr_result,"Methylation_Gene_MR_Result.txt"),quote=FALSE,row.names=FALSE,sep="\t")













##########################合并结果


#eGFR
meth_gene=fread(paste0(mr_result,"Methylation_Gene_MR_Result.txt"))




gene_eGFR=fread(paste0(mr_result,"Gene_eGFR_MR_Result.txt"))
name=names(gene_eGFR)
names(gene_eGFR)=paste0("eGFR_",name)
names(gene_eGFR)[10]="SYMBOL"
gene_eGFR$eGFR.fdr=p.adjust(gene_eGFR$eGFR_pval, method = "fdr")

index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")
names(sentinels)[1]="CpG"
assoc=sentinels[,c("CpG","Bacon_b_DM","UCSC_RefGene_Name","UCSC_RefGene_Group","CHR","MAPINFO")]
names(meth_gene)[4]="CpG"

result=merge(meth_gene,sentinels,by="CpG")
 result$gene.fdr=p.adjust(result$pval, method = "fdr")


my_result=merge(result,gene_eGFR,by="SYMBOL",all.x=T)


library(dplyr)

my_result1= my_result %>% select (SYMBOL,CpG,CHR,MAPINFO,nsnp,b,se,pval,Bacon_b_DM,eGFR_b,eGFR_se,eGFR_pval,gene.fdr,eGFR.fdr)


my_result1$dir_final=my_result1$b * my_result$eGFR_b * my_result1$Bacon_b_DM  ###必须＜0

my_result2=subset(my_result1, dir_final <0 & gene.fdr < 0.05 & eGFR.fdr < 0.05)

write.table(my_result2,paste0(table,"Trio_Meth_gene_eGFR"))





#BUN
table="/share/home/liujianjun/ESRD_Project/EWAS_Manuscript_Table/"


meth_gene=fread(paste0(mr_result,"Methylation_Gene_MR_Result.txt"))


gene_BUN=fread(paste0(mr_result,"Gene_BUN_MR_Result.txt"))
name=names(gene_BUN)
names(gene_BUN)=paste0("BUN_",name)
names(gene_BUN)[10]="SYMBOL"


index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")
names(sentinels)[1]="CpG"
assoc=sentinels[,c("CpG","Bacon_b_DM","UCSC_RefGene_Name","UCSC_RefGene_Group","CHR","MAPINFO")]
names(meth_gene)[4]="CpG"

result=merge(meth_gene,sentinels,by="CpG")


my_result=merge(result,gene_BUN,by="SYMBOL",all.x=T)


library(dplyr)

my_result1= my_result %>% select (SYMBOL,CpG,CHR,MAPINFO,nsnp,b,se,pval,Bacon_b_DM,BUN_b,BUN_se,BUN_pval)

my_result1$dir_gene_BUN=my_result1$b * my_result$BUN_b

my_result1$dir_final=my_result1$dir_gene_BUN * my_result1$Bacon_b_DM  ###必须＜0

my_result2=subset(my_result1, dir_final >0 & pval < 0.05 & BUN_pval < 0.05)

write.table(my_result2,paste0(table,"Trio_Meth_gene_BUN"))

























#ESRD

gene_ESRD=fread(paste0(mr_result,"Gene_ESRD_MR_Result.txt"))
name=names(gene_ESRD)
names(gene_ESRD)=paste0("ESRD_",name)
names(gene_ESRD)[10]="SYMBOL"


index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")
names(sentinels)[1]="CpG"
assoc=sentinels[,c("CpG","Bacon_b_DM","UCSC_RefGene_Name","UCSC_RefGene_Group","CHR","MAPINFO")]
names(meth_gene)[4]="CpG"

result=merge(meth_gene,sentinels,by="CpG")


my_result=merge(result,gene_ESRD,by="SYMBOL",all.x=T)


library(dplyr)

my_result1= my_result %>% select (SYMBOL,CpG,CHR,MAPINFO,nsnp,b,se,pval,Bacon_b_DM,ESRD_b,ESRD_se,ESRD_pval)


my_result1$dir_final=my_result1$b * my_result$ESRD_b * my_result1$Bacon_b_DM  ###必须>0

my_result2=subset(my_result1, dir_final >0 & pval < 0.05 & ESRD_pval < 0.05)

write.table(my_result2,paste0(table,"Trio_Meth_gene_ESRD"))
























################################extract the genes of each strategy for DMPs


out="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Result/"
table="/share/home/liujianjun/ESRD_Project/Table/Output/"

###需要被注释的位点
CpG <- c("cg07054804", "cg16980393", "cg26171235", "cg13854688", "cg19526450")
result=data.frame(CpG=CpG)

#Prioritization gene
setwd(out)
gene=fread("All.ESRD_shared_DMP_Gene.txt")


data1=gene
data1$eQTM_blood <- ifelse(rowSums(data1[, c("eqtm_2015", "eqtm_2018_GTP", "eqtm_2018_MESA", "eqtm_2023")]) >= 1, 1, 0)
data1$eQTM_kidney <- ifelse(rowSums(data1[, c("eqtm_2022_kidney")]) >= 1, 1, 0)
data1$Enhancer_linked_gene <- ifelse(rowSums(data1[, c("atac", "epimap", "hic","ABC")]) >= 1, 1, 0)
data1$coloc <- ifelse(rowSums(data1[, c("meQTL_eQTL")]) >= 1, 1, 0)

new_dt1=subset(data1, select = -c(atac, epic_nearest, epimap, hic, eqtm_2015, eqtm_2018_GTP, eqtm_2018_MESA, eqtm_2022_kidney, eqtm_2023, ABC, meQTL_eQTL,group, sum))

###选取eQTM 位点
new_dt2=subset(new_dt1,eQTM_blood > 0 )

new_dt3=setDT(new_dt2)[, .(gene = paste(gene, collapse = ";")), by = CpG]
names(new_dt3)[2]="eQTM_gene"

###合并
result1=merge(result,new_dt3,by="CpG",all.x=T)
 result1=result1[,c("CpG","eQTM_gene")]
write.table(result1,paste0(table,"Complication.DMP_result_eQTMgene.txt"),row.names=F,quote=F)

###选取enhancer 位点
new_dt2=subset(new_dt1,Enhancer_linked_gene > 0 )


new_dt3=setDT(new_dt2)[, .(gene = paste(gene, collapse = ";")), by = CpG]
names(new_dt3)[2]="Enhancer_linked_gene"

##合并
result1=merge(result,new_dt3,by="CpG",all.x=T)
 result1=result1[,c("CpG","Enhancer_linked_gene")]
 write.table(result1,paste0(table,"Complication.DMP_result_Enhancer_linked_gene.txt"),row.names=F,quote=F)


###选取coloc 位点
new_dt2=subset(new_dt1,coloc > 0 )


new_dt3=setDT(new_dt2)[, .(gene = paste(gene, collapse = ";")), by = CpG]
names(new_dt3)[2]="Colocalized_genes"

##合并
result1=merge(result,new_dt3,by="CpG",all.x=T)
 result1=result1[,c("CpG","Colocalized_genes")]
 write.table(result1,paste0(table,"Complication.DMP_result_Colocalized_genes.txt"),row.names=F,quote=F)
