
###############################################Gene Prioritization
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
base="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Database/"
out="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Result/"
Imp="/share/home/liujianjun/ESRD_Project/Important_files/"
fi="/share/home/liujianjun/ESRD_Project/Part13.DMP_to_Protein/Final_Results/"


library(data.table)

sentinels <- read.csv(paste0(index, "Full.Rep.DMP.Results"), sep = " ", header = TRUE)

anno=read.table(paste0(meth,"3_col_annotation_EPIC.txt"),header=T,sep=" ")

DMP=subset(anno,SNP %in% sentinels$CpG)
names(DMP)=c("CpG","CHR","BP")

##############Step1. ABC Map: ABC max  hg37
abc=fread(paste0(base,"ABC.Map.whole_blood.txt"))

library(dplyr)

# Assuming your data frames are named DMP and abc

# Step 1: Filter abc based on DMP coordinates
filtered_abc <- abc %>%
  inner_join(DMP, by = c("chr" = "CHR")) %>%
  filter(BP >= start, BP <= end)

# Step 2: Find rows with maximum ABC.Score for each CpG
max_abc <- filtered_abc %>%
  group_by(CpG) %>%
  filter(ABC.Score == max(ABC.Score))
## & ABC.Score >= 0.1

# Step 3: Select relevant columns for output
output <- max_abc %>%
  select(CpG, chr, BP, TargetGene, ABC.Score,CellType)  # Add 'cell type' column if available

# Print the resulting dataframe
print(output)

names(output)=c("CpG", "chr", "BP", "gene", "score", "group")
output$data="ABC"

write.table(output,paste0(out,"ABC.Result.txt"),row.names=F,quote=F)



################Step2. Epimap: score > 0.3 (blood t cell and hsc b cell)  hg37
blood=fread(paste0(base,"links_by_group.blood_t_cell.tsv"))
hsc=fread(paste0(base,"links_by_group.hsc_b_cell.tsv"))


#######将基因版本改为Gene Symbol
library(ChIPseeker)
library(org.Hs.eg.db)
library(GenomicFeatures)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library('clusterProfiler')


library(dplyr)

# 假设 DMP 和 blood 已加载到数据框中

# 步骤 1: 根据 DMP 的 CpG 定位到 blood 数据
 blood=blood[,c(1,3,4,5,6,7)]
result <- blood %>%
  inner_join(DMP, by =c("chr" = "CHR"), relationship = "many-to-many") %>%  # 匹配 chr
  filter(start <= BP, end >= BP) %>%  # BP 位于 start 和 end 范围内
  filter(score >= 0.3) %>%  # 筛选 score >= 0.5
  select(CpG, chr, BP, gene, score, group)  # 选择需要的列

  gene<- bitr(result$gene, fromType = "ENSEMBL", toType=c("SYMBOL"),OrgDb = org.Hs.eg.db)  ###转换Gene Symbol
  
  library(dplyr)

# 使用 left_join 将 result 和 gene 数据框连接
result_with_blood <- result %>%
  left_join(gene, by = c("gene" = "ENSEMBL")) %>%  # 连接两个数据框
  mutate(gene = SYMBOL) %>%  # 用 SYMBOL 列替换 gene 列
  select(-SYMBOL)  # 移除 SYMBOL 列，保留替换后的 gene 列

# 输出结果
print(result_with_blood)

###############hsc

 hsc=hsc[,c(1,3,4,5,6,7)]
result <- hsc %>%
  inner_join(DMP, by =c("chr" = "CHR"), relationship = "many-to-many") %>%  # 匹配 chr
  filter(start <= BP, end >= BP) %>%  # BP 位于 start 和 end 范围内
  filter(score >= 0.3) %>%  # 筛选 score >= 0.5
  select(CpG, chr, BP, gene, score, group)  # 选择需要的列

  gene<- bitr(result$gene, fromType = "ENSEMBL", toType=c("SYMBOL"),OrgDb = org.Hs.eg.db)  ###转换Gene Symbol
  
  library(dplyr)

# 使用 left_join 将 result 和 gene 数据框连接
result_with_hsc <- result %>%
  left_join(gene, by = c("gene" = "ENSEMBL")) %>%  # 连接两个数据框
  mutate(gene = SYMBOL) %>%  # 用 SYMBOL 列替换 gene 列
  select(-SYMBOL)  # 移除 SYMBOL 列，保留替换后的 gene 列

# 输出结果
print(result_with_hsc)


#########合并blood和hsc的结果

all=as.data.frame(rbind(result_with_blood,result_with_hsc))
names(all)=c("CpG", "chr", "BP", "gene", "score", "group")
all$data="Epimap"

all <- all %>%
  mutate(group = gsub(" & ", "_", group))
write.table(all,paste0(out,"Epimap.Result.txt"),row.names=F,quote=F)


#######################Hi-C annotation : clusterPostProb > 0.78  hg37

#########需用到两个输入文件：1）增强子和启动子名称的文件：PCHi-C_17BloodCell_ActivePromoterEnhancerLinks.tsv  2）增强子名称对应的基因的文件：PCHiC_peak_matrix_cutoff5.tsv
########################注意oeSt和oeEnd是增强子区域，baitSt和baitEnd是启动子区域！！！！！诱饵序列
#########Step1: 先找到每个DMP所对应的oeID和baitID

hic=fread(paste0(base,"PCHi-C_17BloodCell_ActivePromoterEnhancerLinks.tsv"))
matrix=fread(paste0(base,"PCHiC_peak_matrix_cutoff5.tsv"))

library(data.table)
library(dplyr)

names(hic)[9]="Celltype"
ID= hic %>%
       inner_join(DMP,by=c("oeChr" = "CHR"),relationship = "many-to-many") %>%
       filter(oeSt <= BP, oeEnd >= BP) %>%
       select(CpG,oeChr,BP,oeID,baitID,Celltype)

result= ID %>%
        inner_join(matrix,by=c("oeID","baitID")) %>%
        filter(clusterPostProb > 0.78) %>%
        select(CpG,oeChr.x,BP,baitName,clusterPostProb,Celltype)

names(result)=c("CpG", "chr", "BP", "gene", "score", "group")
result$data="hic"

write.table(result,paste0(out,"Hic.Result.txt"),row.names=F,quote=F)

##############################################snATAC annotation : corr score > 0.45

atac=fread(paste0(base,"Hematopoeisis-Cicero-E-P_TxDb_hg19.tsv"))

library(data.table)

library(dplyr)

result=atac %>%
        inner_join(DMP,by=c("seqnames"="CHR"),relationship = "many-to-many") %>%
        filter(BP>=start,BP<=end, cor >= 0.45) %>%
        select(CpG,seqnames,BP,linkedGene,cor,nearestGene)

names(result)=c("CpG", "chr", "BP", "gene", "score", "group")
result$group="hematopoiesis"
result$data="atac"

write.table(result,paste0(out,"ATAC.Result.txt"),row.names=F,quote=F)





################################Nearest gene annotation
near=fread(paste0(Imp,"Slim.MethylationEPIC_v-1-0_B4.txt"))

names(near)[1]="CpG"

library(dplyr)

 result= near %>%
         inner_join(DMP,by="CpG")  %>%
         select(CpG,CHR.x,BP,UCSC_RefGene_Name,Strand,UCSC_RefGene_Group)
        
names(result)=c("CpG", "chr", "BP", "gene", "score", "group")
result$score="missing"
result$data="epic_nearest"

#result <- result %>%
#  mutate(gene = ifelse(gene == "", NA, gene)) %>%
# mutate(group=ifelse(gene == "",NA,group))

write.table(result,paste0(out,"EPIC_Nearest.Result.txt"),row.names=F,quote=F)




#############################eQTM annotation

##########2015 CE eQTM
eqtm_2015=fread(paste0(base,"2015_09_02_cis_eQTMsFDR0.05-CpGLevel_blood.txt"))
names(eqtm_2015)[2]="CpG"


library(dplyr)

 result= eqtm_2015 %>%
         inner_join(DMP,by="CpG")  %>%
         select(CpG,CHR,BP,HGNCName,IncludedDatasetsCorrelationCoefficient,DatasetsNrSamples)
names(result)=c("CpG", "chr", "BP", "gene", "score", "group")
result$data="eqtm_2015"


write.table(result,paste0(out,"Blood_eQTM_2015.Result.txt"),row.names=F,quote=F)



##########2018 BMC GTP cohort eQTM
eqtm_2018=fread(paste0(base,"2018_BMC_GTP.cohort_blood_eQTM.txt"))
names(eqtm_2018)[2]="CpG"


library(dplyr)

 result= eqtm_2018 %>%
         inner_join(DMP,by="CpG")  %>%
         select(CpG,CHR,BP,annot.gene,beta,status)
names(result)=c("CpG", "chr", "BP", "gene", "score", "group")
result$data="eqtm_2018_GTP"


write.table(result,paste0(out,"Blood_eQTM_2018GTP.Result.txt"),row.names=F,quote=F)


##########2018 BMC MESA cohort eQTM
eqtm_2018=fread(paste0(base,"2018_BMC_MESA.cohort_blood_eQTM.txt"))
names(eqtm_2018)[1]="CpG"


library(dplyr)

 result= eqtm_2018 %>%
         inner_join(DMP,by="CpG")  %>%
         select(CpG,CHR,BP,annot.gene,beta,status)
names(result)=c("CpG", "chr", "BP", "gene", "score", "group")
result$data="eqtm_2018_MESA"
ID=which(nchar(result$group) > 8)
result$group[ID] <- "missing"


write.table(result,paste0(out,"Blood_eQTM_2018MESA.Result.txt"),row.names=F,quote=F)

##########2023 SR eQTM (science report)
cis_eqtm_2023=fread(paste0(base,"2023_SR_FHS_Top10000cis_eQTM.txt"),sep="\t")
cis_eqtm_2023$class="cis"
trans_eqtm_2023=fread(paste0(base,"2023_SR_FHS_Top10000trans_eQTM.txt"),sep="\t")
trans_eqtm_2023$class="trans"
eqtm_2023=rbind(cis_eqtm_2023,trans_eqtm_2023)

names(eqtm_2023)[1]="CpG"
names(eqtm_2023)[9]="Gene"
names(eqtm_2023)[11]="Estimate"
names(eqtm_2023)[14]="Pvalue"

library(dplyr)

 result= eqtm_2023 %>%
         inner_join(DMP,by="CpG")  %>%
         select(CpG,CHR,BP,Gene,Estimate,class)
names(result)=c("CpG", "chr", "BP", "gene", "score", "group")
result$data="eqtm_2023"


write.table(result,paste0(out,"Blood_eQTM_2023.Result.txt"),row.names=F,quote=F)




##########2022 kidney eQTM
eqtm_2022_kidney=fread(paste0(base,"2022_NG_Kidney.eQTM.FDR0.05.txt"),sep="\t")


library(dplyr)

 result= eqtm_2022_kidney %>%
         inner_join(DMP,by="CpG")  %>%
         select(CpG,CHR,BP,GeneSymbol,Beta,Pvalue)
names(result)=c("CpG", "chr", "BP", "gene", "score", "group")
result$data="eqtm_2022_kidney"
result$group=as.character(result$group)
result$group <- "missing"


write.table(result,paste0(out,"Kidney_eQTM_2022.Result.txt"),row.names=F,quote=F)






################Coloc meQTL-eQTL gene
##Input meQTL-eQTL files:
#CpG	gene	score	group	data
#cg23024467	ABHD12B	1.000 	eQTLGEN	meQTL_eQTL
#cg04816311	ADAP1	1.000 	eQTLGEN	meQTL_eQTL


coloc=fread(paste0(out,"meQTL.eQTL.Input.txt"))


library(dplyr)

 result= coloc %>%
         inner_join(DMP,by="CpG")  %>%
         select(CpG,CHR,BP,gene,score,group,data)
names(result)=c("CpG", "chr", "BP", "gene", "score", "group","data")



write.table(result,paste0(out,"meQTL_eQTL.Result.txt"),row.names=F,quote=F)






##############################Merge all the result
setwd("/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Result")

# 获取当前目录下所有Result.txt文件的路径
files <- list.files(pattern = "Result.txt$")

# 读取所有文件并使用rbind合并
data_list <- lapply(files, read.csv, sep=" ",header = TRUE,quote = "")  # 假设文件有列名
merged_data <- do.call(rbind, data_list)

# 查看合并后的数据
head(merged_data)


back_up=merged_data
back_up$data <- ifelse(grepl("^eqtm", back_up$data, ignore.case = TRUE), "eQTM", back_up$data)

inp="/share/home/liujianjun/ESRD_Project/Plot.dir/Input/"
write.table(back_up,paste0(inp,"Shared_ESRD_Gene_Input.txt"),row.names=F,quote=F,sep="\t")   ##该文件用来作基因注释的图

############### 按分号分隔并将基因拆成多行
library(tidyr)

merged_data <- merged_data %>%
  separate_rows(gene, sep = ";")  # 按分号分隔并将基因拆成多行

# 查看结果
head(merged_data)

###########重构表格
final=merged_data
source=unique(final$data)

library(dplyr)


# 处理数据
merged_data_processed <- merged_data %>%
  separate_rows(data, sep = ";") %>%
  group_by(CpG, gene) %>%
  summarise(
    atac = as.integer(any(data == "atac")),
    epic_nearest = as.integer(any(data == "epic_nearest")),
    epimap = as.integer(any(data == "Epimap")),
    hic = as.integer(any(data == "hic")),
    eqtm_2015 = as.integer(any(data == "eqtm_2015")),
    eqtm_2018_GTP = as.integer(any(data == "eqtm_2018_GTP")),
    eqtm_2018_MESA = as.integer(any(data == "eqtm_2018_MESA")),
    eqtm_2022_kidney = as.integer(any(data == "eqtm_2022_kidney")),
    eqtm_2023 = as.integer(any(data == "eqtm_2023")),
    ABC = as.integer(any(data == "ABC")),
    meQTL_eQTL = as.integer(any(data == "meQTL_eQTL")),  # 新增的列
    source = paste(
      ifelse(any(data == "atac"), "atac", ""),
      ifelse(any(data == "epic_nearest"), "epic_nearest", ""),
      ifelse(any(data == "Epimap"), "Epimap", ""),
      ifelse(any(data == "hic"), "hic", ""),
      ifelse(any(data == "eqtm_2015"), "eqtm_2015", ""),
      ifelse(any(data == "eqtm_2018_GTP"), "eqtm_2018_GTP", ""),
      ifelse(any(data == "eqtm_2018_MESA"), "eqtm_2018_MESA", ""),
      ifelse(any(data == "eqtm_2022_kidney"), "eqtm_2022_kidney", ""),
      ifelse(any(data == "eqtm_2023"), "eqtm_2023", ""),
      ifelse(any(data == "ABC"), "ABC", ""),
      ifelse(any(data == "meQTL_eQTL"), "meQTL_eQTL", ""),  # 在 source 中添加 meQTL_eQTL
      sep = ";"
    ),
    group = dplyr::first(group),
    .groups = "drop"
  )



# 查看结果
final=as.data.frame(merged_data_processed)

final <- final %>%
  mutate(sum = rowSums(select(., atac:meQTL_eQTL), na.rm = TRUE))

# 查看结果
head(final)
final=as.data.frame(final)
final$gene <- ifelse(final$gene == "", "missing", final$gene)
final$group <- ifelse(final$group == "", "missing", final$group)
sub=subset(final,! gene %in% "missing")
sub$source <- gsub(";+", ";", sub$source)
sub$source <- sub("^;", "", sub$source)
sub$source <- sub(";$", "", sub$source)

write.table(sub,paste0(out,"All.ESRD_shared_DMP_Gene.txt"),row.names=F,quote=F,sep="\t")






##################Relationship between annotated gene and eGFR (Published TWAS) and Protein

base="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Database/"
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
base="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Database/"
out="/share/home/liujianjun/ESRD_Project/Part7.Annotation_Genes/Result/"
Imp="/share/home/liujianjun/ESRD_Project/Important_files/"
fi="/share/home/liujianjun/ESRD_Project/Part13.DMP_to_Protein/Final_Results/"

library(data.table)
library(dplyr)

final=fread(paste0(out,"All.ESRD_shared_DMP_Gene.txt"))
twas=fread(paste0(base,"kidney_functon_TWAS.txt"))
list=fread(paste0(Imp,"Kidney_Trait_CpGs.txt"))  #既往报道的与肾功能或肾病相关的CpG


twas1 = twas %>%
  select(V1, V2, V4, V11, V15, V16, V17) %>%
  rename(
    trait = V1,
    tissue = V2,
    gene = V4,
    beta = V11,
    direction_same = V15,
    P = V16,
    FDR = V17
  )
head(twas1)

#twas_sig=subset(twas1,FDR == "TRUE")
twas_sig=subset(twas1,P < 0.05/nrow(DMP))

for ( i in 1:nrow(final)) { ifelse (final$gene[i] %in% twas_sig$gene, final$TWAS[i] <- "*" , final$TWAS[i] <- "-")}
head(final)
for ( i in 1:nrow(final)) { ifelse (final$CpG[i] %in% list$CpG, final$EWAS[i] <- "Known" , final$EWAS[i] <- "Novel")}
head(final)
for ( i in 1:nrow(final)) { ifelse (final$gene[i] %in% pro$Protein, final$PWAS[i] <- "*" , final$PWAS[i] <- "-")}
head(final)


slim=fread(paste0(Imp,"Slim.MethylationEPIC_v-1-0_B4.txt"))
slim=slim[,1:3]
names(slim)[1]="CpG"
all=merge(final,slim,by="CpG",all.x=TRUE)

write.table(all,paste0(out,"Final.ESRD_DMP_Gene.txt"),row.names=F,quote=F,sep="\t")




