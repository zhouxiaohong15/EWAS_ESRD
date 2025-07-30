############################ Functional analysis
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
func="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/"
inp="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/Input_files/"
epi="/share/home/liujianjun/BGI_data/BGItrans/PartH.Target_gene_annotation/epimap_link/Epimap.by.Pergroup/"
set="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/Important_files/Backgroud.1000CpGs.sets/"
code="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/01.Code/"
batch="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/01.Code/Batch.Code/"

###########################Step1. Calculate mean and sd for each CpG 
#usr/bin/Rscript

library(data.table)
library(matrixStats)

bval=as.matrix(fread(paste0(meth,"校正scan和batch后的bval")),rownames=1)


# 计算每行的均值
row_means <- rowMeans2(bval)

# 计算每行的标准差
row_sds <- rowSds(bval)

result=data.frame(CpG=rownames(bval),mean=row_means,sd=row_sds)

#########读取所有探针的注释文件
anno=read.table(paste0(meth,"3_col_annotation_EPIC.txt"),header=T,sep=" ")

names(anno)[1]="CpG"

Site_anno=anno[which(anno$CpG %in% result$CpG),]

final=merge(Site_anno,result,by="CpG")
rownames(final)=final$CpG

write.table(final,paste0(meth,"All.Anno.CpG.Mean.SD.txt"),row.names=T,quote=F)

#####################Step2 . Generate similar 1000 CpG sets to matching the DMP list

###############Ref:https://github.com/WRScottImperial/Human-adipocyte-5mC-obesity/blob/master/03_Genomic_Enrichment_Analyses.R 

##############################生成1000个背景的CpG set
##############加载DMP

sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ",sep=" ")
names(sentinels)[1]="CG"


############加载所有位点的平均值和标准差
anno=read.table(paste0(meth,"All.Anno.CpG.Mean.SD.txt"),header=T)
 names(anno)=c("CG","chr","pos","mean","sd")
rownames(anno)=anno$CG
msd=anno[,c(1,4:5)]
colnames(msd)=c("CG","meanMeth","sdMeth")
rownames(msd)=anno$CG
#msd$CG=final$CG


# generate a background table of all cgs, without i. sentinels and ii. cgs within 5kb of a sentinel

hits=subset(msd,CG %in% sentinels$CG)
#hits$CG=as.character(rownames(hits))

################背景中去除和DMP共定位在10kb范围内的CpG
cis.cgs=c()
for(c in 1:length(hits$CG)){
  
  print(c)
  cg = hits$CG[c]
  cg.coord = anno[anno$CG==cg,]
  cis.cg.coord = anno[anno$chr == cg.coord$chr & abs(anno$pos - cg.coord$pos)<10000,]
  cis.cgs = c(cis.cgs,cis.cg.coord$CG) }

cis.cgs=unique(cis.cgs)

backg=data.frame(msd[!(rownames(msd) %in% c(hits$CG,cis.cgs)),])

#############JJ说后续采用该文件作为背景甲基化位点
write.table(backg,paste0(meth,"All.Shared.DMP.background.CpGs"),quote=F,row.names=T)

backg=read.table(paste0(meth,"All.Shared.DMP.background.CpGs"),header=T,row.names=T)

# set permutation parameters for sliding mean and standard deviation (starting low ),生成10*12的120行矩阵
parameters = expand.grid(seq(0, 0.02, 0.002),seq(0,0.05,0.005),stringsAsFactors = FALSE); 
colnames(parameters)=c("SD_threshold","Mean_threshold")

# generate permutation set (nrows = number of sentinels, ncols =1000 )

matches = matrix(nrow=nrow(hits), ncol=1000)
colnames(matches)=as.character(seq(from =1 , to =1000, 1))
rownames(matches)=hits$CG
matchParameters = matrix(nrow=nrow(hits), ncol=3)
colnames(matchParameters)=c("mean_threshold","sd_threshold","n_matches")
rownames(matchParameters)=hits$CG



# calculate number of matches across the full range of mean sd thresholds
 
for(h in 1:nrow(hits)){
  cpg.match=NA
  len=NA
  hit=hits[h,]
  print(h)  
  for(para in 1:nrow(parameters)){
    n = length(sample(rownames(backg[abs(backg$meanMeth-hit$meanMeth)<parameters[para,2] & abs(backg$sdMeth-hit$sdMeth)<parameters[para,1],]))  )
    if(para==1) {sel = c(parameters[para,2],parameters[para,1],n)}
    if(para>1) {sel = rbind(sel,c(parameters[para,2],parameters[para,1],n))}
  }
  colnames(sel) = c("mean","sd","n_matches")
  
  # select mean sd thresholds;  n >1000 and use top (smallest) thresholds
   sel = sel[which(sel[,3]>1000),]
  meanthres = sel[1,1]
  sdthres = sel[1,2]
  nthresh = sel[1,3]
  
  # identify matches for each hit, filling by columns without replacement
  for(p in 1:1000){
    set.seed(p)
    print(p)
    if(exists('cgm')){rm(cgm)}
    cgms=sample(rownames(backg[abs(backg$meanMeth-hit$meanMeth)<meanthres & abs(backg$sdMeth-hit$sdMeth)<sdthres,])) 
    cgm=setdiff(cgms, cpg.match)[1]
    if(!is.na(cgm)){
      cpg.match=c(cpg.match,cgm)
    }else{
      cpg.match=c(cpg.match,(as.character(rownames(hit))))  # temporary solution
      print('no matching marker found')
    }
  }
  matches[h,]=na.omit(cpg.match)  
  matchParameters[h,] = c(meanthres,sdthres,nthresh)
  }

setwd(inp)
save(matches,file="matches.Rdata")
save(matchParameters,file="matchParameters.Rdata")

###########################Step3 Generate batch R script files for permutation + binomial test

####group="adipose	blood_t_cell	bone	brain	digestive	endocrine	endothelial	epithelial    esc	es	eye	heart	hsc_b_cell	ipsc	kidney	liver	lung	lymphoblastoid	mesench	musclemyosat	neurosph	pancreas	placenta	pns	reproductive	smmuscle	spleen	stromathymus	urinary";


 #!/usr/bin/env Rscript
library(data.table)
setwd(inp)

sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")
bed=fread(paste0(epi,"links_by_group.$i.tsv"))  #####$i替换为$group, group="adipose	blood_t_cell	bone	brain	digestive	endocrine	endothelial	epithelial    esc	es	eye	heart	hsc_b_cell	ipsc	kidney	liver	lung	lymphoblastoid	mesench	musclemyosat	neurosph	pancreas	placenta	pns	reproductive	smmuscle	spleen	stromathymus	urinary";

all= read.table(paste0(meth,"All.Anno.CpG.Mean.SD.txt"),header=T)
setwd(inp)
load("matches.Rdata")
matches=as.data.frame(matches)

#############################Background CpG sites
result=c()
file=c()
length=c()


########读取每一个background CpG set, 并注释位置信息
for (j in 1:ncol(matches)) { 
result=c(); 
print (paste0("Background.CpG.set",j));
anno=subset(all, CpG %in% matches[,j] )[,1:3];  
names(anno)=c("CpG","CHR","BP")
name=paste0("Backgroud.CpG.list",j)
write.table(anno,paste0(set,name),row.names=F,quote=F)
 }


############第一个for:每次读取5个 background CpG sets (i in 1:5)
####################第二个for:对每个background CpG set中的nrow(sentinels)个CpG进行 open chromatin的注释(h in 1:nrow(sentinels))，并得到与open chromatin overlap的数量;

index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
func="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/"
inp="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/Input_files/"
epi="/share/home/liujianjun/BGI_data/BGItrans/PartH.Target_gene_annotation/epimap_link/Epimap.by.Pergroup/"
set="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/Important_files/Backgroud.1000CpGs.sets/"
back="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/02.Result_files/Background_CpG_Results/"

library(data.table)
setwd(inp)

sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")
##group="adipose	blood_t_cell	bone	brain	digestive	endocrine	endothelial	epithelial    esc	es	eye	heart	hsc_b_cell	ipsc	kidney	liver	lung	lymphoblastoid	mesench	musclemyosat	neurosph	pancreas	placenta	pns	reproductive	smmuscle	spleen	stromathymus	urinary";
bed=fread(paste0(epi,"links_by_group.$i.tsv"))  #####要替换，s/\$i/${group}/g,group为组织的list,如上
all= read.table(paste0(meth,"All.Anno.CpG.Mean.SD.txt"),header=T)
setwd(inp)
load("matches.Rdata")
matches=as.data.frame(matches)

result=c()
file=c()
length=c()

numbers <- sapply(seq(1, 491, by = 10), function(x) paste0(x, ":", x + 9))  #####一共100个 (100*10=1000)
# 处理 numbers 的每一个范围，要替换,s/\$g/$g/g,g为1到100，g=seq(1:100)
first_range <- unlist(strsplit(numbers[$g], ":"))  # 例如将 "1:10" 拆分为 c("1", "10"),要替换,s/\$g/$g/g
# 将范围转换为数字并生成序列
seq_range <- as.numeric(first_range[1]):as.numeric(first_range[2])
 
result <- c()
lengths <- c()

for (g in seq_range) {   # seq_range 例如 1:10, 11:20, ..., 991:1000
  name <- paste0("Backgroud.CpG.list", g)
  anno <- as.data.frame(fread(paste0(set, name)))
  result <- data.frame()  # 初始化 result 数据框

  for (h in 1:nrow(anno)) {
    CpG.open <- subset(bed, chr %in% anno$CHR[h] & start < anno$BP[h] & end > anno$BP[h])
    if (nrow(CpG.open) > 1) {
      result <- rbind(result, CpG.open[1, ])  # 添加第一行数据
    } else {
      result <- rbind(result, CpG.open)       # 添加全部数据
    }
  }

  # 记录当前组的 CpG 数量
  nCpG <- nrow(result)
  lengths <- c(lengths, nCpG)
}

# 创建最终数据框
file <- as.data.frame(t(lengths))  # 转置为数据框
rownames(file) <- "Counts of CpGs within enhancer"
colnames(file) <- paste0("set", seq_range)  # 列名与 seq_range 对应

# 写入文件
modified_string <- gsub(":", "to", as.character(numbers[$g]))
name <- paste0("$i.Back_Num_Set_", modified_string)
write.table(file, paste0(back, name), quote = FALSE, row.names = TRUE)




#############################Target CpG sites

index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
func="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/"
inp="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/Input_files/"
epi="/share/home/liujianjun/BGI_data/BGItrans/PartH.Target_gene_annotation/epimap_link/Epimap.by.Pergroup/"
set="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/Important_files/Backgroud.1000CpGs.sets/"
back="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/02.Result_files/Background_CpG_Results/"
tar="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/02.Result_files/Target_CpG_Results/"

library(data.table)
setwd(inp)


###########################Step 4. Calculate the counts of DMP within active enhancers across multiple tissues


# 读取必要的文件
sentinels <- read.csv(paste0(index, "Full.Rep.DMP.Results"), sep = " ", header = TRUE)
epic <- read.table(paste0(meth, "3_col_annotation_EPIC.txt"), header = TRUE, sep = " ")
all <- read.table(paste0(meth, "All.Anno.CpG.Mean.SD.txt"), header = TRUE)
bed=fread(paste0(epi,"links_by_group.$i.tsv")) 
setwd(inp)

# 提取感兴趣的子集
sub <- subset(epic, CpG %in% sentinels$CpG)
anno <- sub[, 1:3]

# 定义 group 变量
group <- "adipose blood_t_cell bone brain digestive endocrine endothelial epithelial esc es_deriv eye heart hsc_b_cell ipsc kidney liver lung lymphoblastoid mesench muscle myosat neurosph pancreas placenta_eem pns reproductive smmuscle spleen stromal thymus urinary"

# 将字符串分割为一个向量
group_list <- unlist(strsplit(group, "\\s+"))

# 初始化结果
length <- c()
# 遍历每一行进行匹配
for (g in group_list) {
  # Read the data file for the current group
  bed <- fread(paste0(epi, "links_by_group.", g, ".tsv"))
  result <- data.frame()
  # Loop through each row in sentinels
  for (h in seq_len(nrow(anno))) {
    # Filter rows matching the conditions
    CpG.open <- bed[chr == anno$CHR[h] & 
                    start < anno$BP[h] & 
                    end > anno$BP[h]]
    
    # Add the first matching row to the result if there are matches
    if (nrow(CpG.open) > 0) {
      result <- rbind(result, CpG.open[1, ], fill = TRUE)
    }
  }
 length=c(length,nrow(result))
}

# 统计结果
num <- data.frame(
  Sample.ID = group_list,
  Open.num = length,
  Close.num = nrow(sentinels) - length
)

# 写入结果
write.table(num, paste0(tar, "Open.DMP.All_Tissue.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)




###########Calculate: Binomial test

index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
func="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/"
inp="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/Input_files/"
epi="/share/home/liujianjun/BGI_data/BGItrans/PartH.Target_gene_annotation/epimap_link/Epimap.by.Pergroup/"
set="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/Important_files/Backgroud.1000CpGs.sets/"
back="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/02.Result_files/Background_CpG_Results/"
tar="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/02.Result_files/Target_CpG_Results/"


library(data.table)

sentinels <- read.csv(paste0(index, "Full.Rep.DMP.Results"), sep = " ", header = TRUE)
 
num=read.table(paste0(enr,"Open.DMP.All_Tissue.txt"),header=T)
names(num)[1]="Tissue"
num$Probability=num$Open.num/nrow(sentinels)
file=read.table(paste0(enr,"Background.output.txt"),header=T)

file$Probability = file$Open.Counts/(nrow(sentinels)*500)  ####500个CpG sets

file=as.data.frame(file)

# 初始化存储结果的变量
result <- data.frame(Tissue = character(), 
                     Enrichment = numeric(), 
                     P_value = numeric(), 
                     stringsAsFactors = FALSE)

# 遍历每行数据
for (i in 1:nrow(num)) {
    # 计算二项检验
    binomial <- binom.test(num$Open.num[i], n = nrow(sentinels), p = file$Probability[i])
    
    # 计算富集值
    enrichment <- (num$Open.num[i] / nrow(sentinels)) / file$Probability[i]
    
    # 将结果存储为一个新数据框
    new_row <- data.frame(Tissue = file$Tissue[i],
                          Enrichment = enrichment,
                          P_value = binomial$p.value,
                          stringsAsFactors = FALSE)
    
    # 合并新数据行到结果数据框
    result <- rbind(result, new_row)
}

# 查看结果
print(result)


write.table(result,paste0(enr,"Enrichment.Active.enhancer.result.txt"),quote=F) 











########################Multiple traits enrichment analysis (input file formats: CpG,Beta,SE,Pvalue) 
func="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/eGFR_EWAS_EPIC/"  ##要改
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
set="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/Important_files/Backgroud.1000CpGs.sets/"
inp="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/Input_files/"
back="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/02.Result_files/Background_CpG_Results/"
tar="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/02.Result_files/Target_CpG_Results/"
enr="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/02.Result_files/Enrichment_Results/"

library(data.table)
setwd(inp)

sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")
##trait="eGFR1 eGFR2 Albuminuria CRP basic_diabetes basic_heart_disease HbA1c PAH basic_alzheimers basic_breast_cancer basic_chronic_pain basic_CKD basic_colorectal_cancer basic_COPD  basic_lung_cancer basic_osteoarthritis basic_parkinsons basic_prostate_cancer basic_rheumatoid_arthritis basic_stroke	CKD Gluse CRP AD ... "; 记得补充完整
bed=fread(paste0(func,"Data.EWAS.$i.txt"))  #####要替换，s/\$i/${trait}/g,group为trait的list,如上
#all= read.table(paste0(meth,"All.Anno.CpG.Mean.SD.txt"),header=T)
setwd(inp)
load("matches.Rdata")
matches=as.data.frame(matches)

result=c()
file=c()
length=c()


numbers <- sapply(seq(1, 491, by = 10), function(x) paste0(x, ":", x + 9))  #####一共100个 (100*10=1000)
# 处理 numbers 的每一个范围，要替换,s/\$g/$g/g,g为1到100，g=seq(1:100)
first_range <- unlist(strsplit(numbers[$g], ":"))  # 例如将 "1:10" 拆分为 c("1", "10"),要替换,s/\$g/$g/g
# 将范围转换为数字并生成序列
seq_range <- as.numeric(first_range[1]):as.numeric(first_range[2])
 
result <- c()
lengths <- c()
lengths.analized=c()


for (g in seq_range) {   # seq_range 例如 1:10, 11:20, ..., 991:1000
  name <- paste0("Backgroud.CpG.list", g)
  anno <- as.data.frame(fread(paste0(set, name)))
  result <- data.frame()  # 初始化 result 数据框
  CpG.analized <- subset(anno,CpG %in% bed$CpG)

  for (h in 1:nrow(anno)) {
    CpG.open <- subset(bed, CpG %in% anno$CpG[h] & Pvalue < 0.05/nrow(sentinels))
    
    if (nrow(CpG.open) > 1) {
      result <- rbind(result, CpG.open[1, ])  # 添加第一行数据
    } else {
      result <- rbind(result, CpG.open)       # 添加全部数据
    }
  }

  # 记录当前组的 CpG 数量
  nCpG <- nrow(result)
  lengths <- c(lengths, nCpG)
  nCpG.analized <- nrow(CpG.analized)
  lengths.analized <- c(lengths.analized,nCpG.analized)
}

# 创建最终数据框
file <- as.data.frame(t(lengths))  # 转置为数据框
file=rbind(file,as.data.frame(t(lengths.analized)))
rownames(file) <- c("Counts of overlaped CpGs", "Counts of analyzed CpGs")
colnames(file) <- paste0("set", seq_range)  # 列名与 seq_range 对应

# 写入文件
modified_string <- gsub(":", "to", as.character(numbers[$g]))
name <- paste0("$i.Back_Num_Set_", modified_string)
write.table(file, paste0(back, name), quote = FALSE, row.names = TRUE)




###########################Step 4. Calculate the counts of DMP associated with with other trait
func="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/eGFR_EWAS_EPIC/"  ##要改
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
set="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/Important_files/Backgroud.1000CpGs.sets/"
inp="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/Input_files/"
back="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/02.Result_files/Background_CpG_Results/"
tar="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/02.Result_files/Target_CpG_Results/"
enr="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/02.Result_files/Enrichment_Results/"


# 读取必要的文件
sentinels <- read.csv(paste0(index, "Full.Rep.DMP.Results"), sep = " ", header = TRUE)
epic <- read.table(paste0(meth, "3_col_annotation_EPIC.txt"), header = TRUE, sep = " ")
#all <- read.table(paste0(meth, "All.Anno.CpG.Mean.SD.txt"), header = TRUE)
setwd(inp)

anno <- sentinels

# 定义 trait 变量
trait <- "eGFR1 eGFR2 Albuminuria CRP basic_diabetes basic_heart_disease HbA1c PAH basic_alzheimers basic_breast_cancer basic_chronic_pain basic_CKD basic_colorectal_cancer basic_COPD  basic_lung_cancer basic_osteoarthritis basic_parkinsons basic_prostate_cancer basic_rheumatoid_arthritis basic_stroke"
#######要改


# 将字符串分割为一个向量
group_list <- unlist(strsplit(trait, "\\s+"))

# 初始化结果
length <- c()
lengths.analized=c()
compare <- c()


# 遍历每一行进行匹配
for (g in group_list) {
  # Read the data file for the current group
  bed <- fread(paste0(func, "Data.EWAS.", g, ".txt"))
  nCpG.analized <- nrow(subset(anno,CpG %in% bed$CpG))
  result <- data.frame()
  # Loop through each row in sentinels
  for (h in seq_len(nrow(anno))) {
    # Filter rows matching the conditions
    CpG.open <- bed[ CpG %in% anno$CpG[h] & Pvalue < 0.05/nrow(sentinels)]
    CpG.open$ID=g
    # Add the first matching row to the result if there are matches
    if (nrow(CpG.open) > 0) {
      result <- rbind(result, CpG.open[1, ])
    }
  }
 length=c(length,nrow(result))
 lengths.analized <- c(lengths.analized,nCpG.analized)

effect=merge(result, anno[,c("CpG","Bacon_b_meta","Bacon_P_meta")],by="CpG")
name=paste0("Compare.Beta.",g)
write.table(effect, paste0(enr, name), quote = FALSE, row.names = FALSE, col.names = TRUE)

    }


# 统计结果
num <- data.frame(
  Sample.ID = group_list,
  Open.num = length,
  Close.num = lengths.analized - length
)

# 写入结果
write.table(num, paste0(tar, "Overlap.DMP.All_Trait.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)




###########Calculate: Binomial test
func="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/eGFR_EWAS_EPIC/"  ##要改
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
set="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/Important_files/Backgroud.1000CpGs.sets/"
inp="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/Input_files/"
back="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/02.Result_files/Background_CpG_Results/"
tar="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/02.Result_files/Target_CpG_Results/"
enr="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/02.Result_files/Enrichment_Results/"


library(data.table)

sentinels <- read.csv(paste0(index, "Full.Rep.DMP.Results"), sep = " ", header = TRUE)
 
num=read.table(paste0(tar,"Overlap.DMP.All_Trait.txt"),header=T)
names(num)[1]="Trait"
num$all=num$Open.num + num$Close.num
num$Probability=num$Open.num/num$all


file=read.table(paste0(back,"Background.output.txt"),header=T)
file$all=file$Overlaped.CpGs + file$Analyzed.CpGs
file$Probability = file$Overlaped.CpGs/file$all  ####500个CpG sets

file=as.data.frame(file)
names(file)[1]="Trait"

# 初始化存储结果的变量
result <- data.frame(Trait = character(), 
                     RR = numeric(), 
                     P_value = numeric(), 
                     stringsAsFactors = FALSE)

# 遍历每行数据
for (i in 1:nrow(num)) {
    # 计算二项检验
    binomial <- binom.test(num$Open.num[i], num$Close.num[i], p = file$Probability[i])
    
    # 计算富集值
    enrichment <- (num$Open.num[i] /  num$all[i]) / file$Probability[i]
    
    # 将结果存储为一个新数据框
    new_row <- data.frame(Trait = file$Trait[i],
                          RR = enrichment,
                          P_value = binomial$p.value,
                          stringsAsFactors = FALSE)
    
    # 合并新数据行到结果数据框
    result <- rbind(result, new_row)
}

# 查看结果
print(result)
result=result[order(result$P_value),]

result=result[order(result$Pvalue),]


write.table(result,paste0(enr,"Enrichment.Active.enhancer.result.txt"),quote=F) 


######################################整理所有的下载的EWAS数据，将名称全部改掉
func="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/eGFR_EWAS_EPIC/"  ##要改



for file in ${func}prevalent_*_gs.csv; do
    base_name=$(basename "$file")  # 获取文件名（去掉路径）
    new_name=$(echo "$base_name" | sed -e 's/prevalent_//g' -e 's/_gs.csv/.txt/g')
    mv "$file" "${func}Data.EWAS.$new_name"
done


  # 清空 Trait.List.txt
echo "" > ${func}Trait.List.txt

# 遍历文件并提取名称
for file in ${func}Data.EWAS.*.txt; do
  # 提取文件名（不含路径）
  base_name=$(basename "$file")
  # 移除前缀和后缀
  single=$(echo "${base_name}" | sed -e 's/Data.EWAS.//g' -e 's/.txt//g')
  # 追加到 Trait.List.txt，不换行
  echo -n " ${single}" >> ${func}Trait.List.txt
done

#######把列表里的 Albuminuria basic_alzheimers basic_breast_cancer basic_chronic_pain basic_CKD  放到最前面

HbA1c： 来自 2020 PNAS sutuzak
Albuminuria: 来自 2020 PNAS sutuzak
fasting glucose fasting insulin : 来自zendo
CRP: 来自zendo  文件crp_full_model_hillaryetal.csv
PAH： 来自Pubmed

####注意格式都必须是： CpG Beta Pvalue Trait


#################################################Protein enrichment analyis

pro="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/compressed-protein-ewas/"

cd ${pro}
for file in MWAS_*.csv; do
  # 使用正则提取文件名中需要的部分
  new_name=$(echo $file | sed -E 's/^MWAS_[0-9]+-[0-9]+_([A-Za-z0-9]+)\.csv$/Data.EWAS.\1.txt/')
  # 重命名文件
  mv "$file" "$new_name"
done

####第一行改为CpG,Beta,SE,Pvalue
cd ${pro}
for file in Data.EWAS.*.txt; do
  # 替换文件的第一行，并保存到临时文件
  sed '1s/.*/CpG,Beta,SE,Pvalue/' "$file" > temp.csv
  # 用临时文件替换原文件
  mv temp.csv "$file"
done


######输出所有基因名字到列表中
cd ${pro}
for file in Data.EWAS.*.txt; do
  # 提取基因名部分并添加到文件 List 中
  gene_name=$(echo "$file" | sed -E 's/^Data\.EWAS\.(.*)\.txt/\1/')
  echo -n "$gene_name " >> Protein.List.txt
done















############Explore the association between DMP and Proteins
#!/bin/bash

# 定义路径变量
func="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/eGFR_EWAS_EPIC/"  ##要改
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
set="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/Important_files/Backgroud.1000CpGs.sets/"
inp="/share/home/liujianjun/ESRD_Project/Part8.Function_Analysis/01.chromatin_accessibility/Input_files/"
back="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/02.Result_files/Background_CpG_Results/"
tar="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/02.Result_files/Target_CpG_Results/"
enr="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/02.Result_files/Enrichment_Results/"
pro="/share/home/liujianjun/ESRD_Project/Part11.Multiple_Trait_EWAS/compressed-protein-ewas/"
protein="/share/home/liujianjun/ESRD_Project/Part13.DMP_to_Protein/Code/"
pro_result="/share/home/liujianjun/ESRD_Project/Part13.DMP_to_Protein/Results/"
fi="/share/home/liujianjun/ESRD_Project/Part13.DMP_to_Protein/Final_Results/"

for i in $(cat "${pro}Protein.List.txt"); do
    echo "
file1=\"${pro}Data.EWAS.${i}.txt\"
file2=\"${index}Full.Rep.DMP.Results\"

# 计算 file2 的行数
rows=\$(awk 'NR>1' \"\$file2\" | wc -l)
threshold=\$(echo \"scale=10; 0.05 / \$rows\" | bc -l)

# 读取 file2 中的 CpG 位点
awk 'NR>1 {print \$1}' \"\$file2\" > ${index}cpg_list.txt
sed  '1iCpG'  ${index}cpg_list.txt > ${index}tmp.txt && mv ${index}tmp.txt ${index}cpg_list.txt

# 筛选 file1 中 Pvalue < threshold 且 CpG 位点在 file2 中，并添加 Protein 列
awk -F',' -v thresh=\"\$threshold\" -v protein=\"$i\" '
    NR==1 {print \$0, \"Protein\"; next}  # 打印列名并添加 Protein 列
    FILENAME==\"${index}cpg_list.txt\" {cpgs[\$1]=1; next} 
    (\$1 in cpgs) && (\$4+0 < 0.0010204081) {print \$0, protein}' \"${index}cpg_list.txt\" \"\$file1\" > ${pro_result}${i}_ESRD_DMP_signal.txt  
    sed -i \"s/ /,/g\" ${pro_result}${i}_ESRD_DMP_signal.txt  
" > ${protein}${i}_ESRD_DMP_signal.sh   # 修正：使用大括号保证变量替换

    # 提交任务（可选，去掉注释以执行）
     bsub -e ${protein}error.$i.log -o ${protein}output.$i.log < ${protein}${i}_ESRD_DMP_signal.sh
done


###注意上述代码记得改阈值0.0010204081 很重要!!!!!!!!!



############将上述结果全部整理，筛选大于等于3的行，证明注释成功
cd ${pro_result}${fi}Annotation_Results.txt
echo "" > 
for file in *; do
  if [ -f "$file" ] && [ $(wc -l < "$file") -ge 3 ]; then
    echo "$file" >> ${fi}Annotation_list.txt
  fi
done

##去除CpG开头的行
echo "" > Merge.Protein.Results.txt
for filename in $(cat "${fi}Annotation_list.txt"); do
sed -i "/^CpG/d" ${pro_result}${filename}
cat ${pro_result}${filename} >> ${fi}Merge.Protein.Results.txt
done;

sed -i "/Protein/d" ${fi}Merge.Protein.Results.txt

sed -i "1iCpG,Beta,SE,Pvalue,Protein" ${fi}Merge.Protein.Results.txt
