#!/usr/bin/env Rscript

library(ff)
library(RnBeads.hg19)
library(methylumi)
library(foreach)
library(doParallel)
options(fftempdir="/share/home/liujianjun/tmp")
getOption("fftempdir")
library(RnBeads)
library(RnBeads.hg19)
rnb.options(disk.dump.big.matrices=TRUE, disk.dump.bigff=TRUE)
rnb.options(logging.disk=TRUE)
rnb.options(enforce.memory.management=TRUE)
rnb.options(identifiers.column="Sample_ID",import.idat.platform='probesEPIC',import.sex.prediction = TRUE)#指定列为Sample_ID



import.dir="/share/home/liujianjun/ESRD_Project/Part0.GDPH_CKD_ESRD_EWAS/00.Import/"

idat.dir="/share/home/liujianjun/EPIC_data_import_xujinjin/SY_ESRD_Rnbeads/00.raw.data/"

sample.annotation="/share/home/liujianjun/ESRD_Project/Part0.GDPH_CKD_ESRD_EWAS/Important_files/GDPH.CKD_ESRD.pheno.csv"


Important="/share/home/liujianjun/ESRD_Project/Part0.GDPH_CKD_ESRD_EWAS/Important_files/"

report.dir = file.path(import.dir, "import_reports")

QC.dir = "/share/home/liujianjun/ESRD_Project/Part0.GDPH_CKD_ESRD_EWAS/01.QC/"

index="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Subtype_Meta_analysis/00.EWAS_results/"




###################################Step1.import:将idat文件读取为rnb.set

#######指定数据读取
data.source <- c(idat.dir, sample.annotation)


######开始导入数据并保存
rnb.set.disk <- rnb.execute.import(data.source=data.source, data.type="idat.dir")
rnb.set=rnb.set.disk
save.rnb.set(rnb.set,paste0(import.dir,"Import_rnb.set"),archive=F)



######################Step2. Normalization and QC

library(ff)
library(RnBeads.hg19)
library(methylumi)
library(foreach)
library(doParallel)
#options(fftempdir="/share/home/liujianjun/tmp")
#getOption("fftempdir")
library(RnBeads)
library(RnBeads.hg19)
library(grid)
rnb.options(disk.dump.big.matrices=TRUE, disk.dump.bigff=TRUE)
rnb.options(logging.disk=TRUE)
rnb.options(enforce.memory.management=TRUE)
rnb.options(identifiers.column="Sample_ID",import.idat.platform='probesEPIC',import.sex.prediction = TRUE)#指定列为Sample_ID

##############Step1.探针质控和数据化标准化路径创建，加载样本质控后的rnb.set

report.dir = file.path(QC.dir, "Norm_reports")


rnb.set=load.rnb.set(paste0(import.dir,"Import_rnb.set"))


#############Step2.设置预处理的参数包括a.过滤不合格位点（greedycut.pvalue，SNP，context,cross.reactive等），b.,两类探针偏差校正(BMIQ标准化处理),c. 根据greedycut来推测call rate，从而评估并过滤不合格样本

rnb.options(filtering.sex.chromosomes.removal=TRUE,identifiers.column="Sample_ID",normalization.method="bmiq",filtering.cross.reactive=TRUE,filtering.snp = "3",filtering.greedycut.pvalue.threshold = 0.01,filtering.missing.value.quantile = 0.1)
rnb.set.unfiltered <- rnb.set
result <- rnb.run.preprocessing(rnb.set.unfiltered, dir.reports=report.dir)

rnb.set <- result$rnb.set
rnb.set_new=rnb.execute.imputation(rnb.set,method ="knn",update.ff = TRUE) #####对missing value进行填补

save.rnb.set(rnb.set_new,paste0(QC.dir,"02.rnb.set_Normalized"),archive=F)



mval=mval(rnb.set,row.names=TRUE)
mval[1:20,1:5]

bval=meth(rnb.set_new,row.names=TRUE)
mval[1:20,1:5]

write.table(mval,paste0(Important,"Rnbset01.Normalized_Mval.txt"))
write.table(mval,paste0(Important,"GDPH.CKD_ESRD.Normalized_Mval.txt"))
write.table(bval,paste0(Important,"Rnbset01.Normalized_Bval.txt"))



##################################Cell Type 细胞类型估计,主要是scanner和control probe PC来获得 residual mval，预测细胞组成。
#library(ff)
library(RnBeads.hg19)
library(methylumi)
library(foreach)
library(doParallel)
#options(fftempdir="/share/home/liujianjun/tmp")
#getOption("fftempdir")
library(RnBeads)
library(RnBeads.hg19)
library(grid)
rnb.options(disk.dump.big.matrices=TRUE, disk.dump.bigff=TRUE)
rnb.options(logging.disk=TRUE)
rnb.options(enforce.memory.management=TRUE)
rnb.options(identifiers.column="Sample_ID",import.idat.platform='probesEPIC',import.sex.prediction = TRUE)#指定列为Sample_ID



library(EpiDISH)
#data(library(EpiDISH))
library(data.table)
data(centDHSbloodCKDC.m)
parallel.isEnabled()




################################### ##################################Control probes PC的提取

rnb.set=load.rnb.set(paste0(QC.dir,"02.rnb.set_Normalized"))

#############4.1 提取质控探针信号
green=qc(rnb.set)$Cy3
red=qc(rnb.set)$Cy5

setwd(Important)
write.table(green,"gre_CtrlProbe.txt",sep=" ")
write.table(red,"red_CtrlProbe.txt",sep=" ")

green=read.table("gre_CtrlProbe.txt",sep=" ",header=T)
red=read.table("red_CtrlProbe.txt",sep=" ",header=T)

##############4.2 读取质控探针注释

anno=read.table(paste0(Important,"850kcontrolprobe_anntation.txt"),sep=" ",header=T)

#############4.3 读取表型文件
setwd(Important)
pheno=read.table("GDPH.CKD_ESRD.pheno.csv",header=T,sep=",")



###################4.4 去掉negtive探针位点,提取norm probes.
red_norm=red[-c(49:53,54:464),]
gre_norm=green[-c(49:53,54:464),]
anno=anno[-c(49:53,54:464),]

###################提取并合并negtive探针位点
red_neg=red[c(54:464),]
gre_neg=green[c(54:464),]
negative=rbind(red_neg,gre_neg)
ID=pheno$Sample_ID
colnames(negative)=ID

#########################合并219个质控探针和信号
newred=cbind(anno,red_norm)
newgre=cbind(anno,gre_norm)
my_red=subset(newred,Evaluate.Red=="+")
my_green=subset(newgre,Evaluate.Green=="+")

all=rbind(my_green,my_red)
all0=all[,-c(1:10)]

#########################对合并后的219个探针进行PCA分析
intensity=t(all0)
intensity_new=scale(intensity,center=F,scale=T)
pca=prcomp(intensity_new)
sum=summary(pca)
head(sum$importance)
pca_new=pca$x[,1:30]

##########################对合并后的Negtive探针进行PCA分析
neg_intensity=t(negative)
neg_intensity_new=scale(neg_intensity,center=F,scale=T)
neg_pca=prcomp(neg_intensity_new)
neg_sum=summary(neg_pca)
head(neg_sum$importance)
neg_pca_new=neg_pca$x[,1:30]

setwd(Important)
write.table(pca_new, "GDPH.CKD_ESRD.CtrlProbe_PC.txt",quote=F)
write.table(neg_pca_new, "GDPH.CKD_ESRD.NegCtrlProbe_PC.txt",quote=F)


##################################细胞类型估计,主要是scanner和control probe PC来获得 residual mval，预测细胞组成。
#library(ff)
library(RnBeads.hg19)
library(methylumi)
library(foreach)
library(doParallel)
#options(fftempdir="/share/home/liujianjun/tmp")
#getOption("fftempdir")
library(RnBeads)
library(RnBeads.hg19)
library(grid)
rnb.options(disk.dump.big.matrices=TRUE, disk.dump.bigff=TRUE)
rnb.options(logging.disk=TRUE)
rnb.options(enforce.memory.management=TRUE)
rnb.options(identifiers.column="Sample_ID",import.idat.platform='probesEPIC',import.sex.prediction = TRUE)#指定列为Sample_ID



library(EpiDISH)
#data(library(EpiDISH))
library(data.table)
data(centDHSbloodCKDC.m)
parallel.isEnabled()


############################读取标准化后的mval，并校正批次效应
setwd(Important)

pheno=read.table("GDPH.CKD_ESRD.pheno.csv",sep=",",header=T)
mval=as.matrix(fread(file="GDPH.CKD_ESRD.Normalized_Mval.txt"),rownames=1)

#########读取batch effect

scan=pheno$scanner


###############################读入219个control探针的PCA

setwd(Important)

ctrl=read.table("GDPH.CKD_ESRD.CtrlProbe_PC.txt",header=T,sep=" ")
ctrl=ctrl[,1:30]
PC=ctrl


########################################合并batch effect
cov=cbind(PC,scan)


#######校正batch effect
library(limma)
res_mval=removeBatchEffect(mval,covariates=cov)
res_bvals=rnb.mval2beta(res_mval)
bvals=res_bvals

setwd(Important)
write.table(res_bvals,"Cell.Estiamete.Residual.Bvals.txt")

###################利用校正后的bval进行细胞组成估计

bvals=res_bvals

out.l <- epidish(beta.m =bvals, ref.m = centDHSbloodCKDC.m, method = "RPC")

cell_type=out.l$estF

#write.table(cell_type,file=paste0(Important,"Rnbset01.Cell_Type.txt"))
write.table(cell_type,file=paste0(Important,"GDPH.CKD_ESRD.Cell.Type.txt")






#####################################################################Correct the methylation level

library(data.table)

mval=as.matrix(fread(file=paste0(Important,"GDPH.CKD_ESRD.Normalized_Mval.txt")),rownames=1)


####################################读入分组
covariance=read.table(paste0(Important,"GDPH.CKD_ESRD.pheno.csv"),sep=",",header=T)
Sample_Group=covariance$Sample_Group

###############################读入219个control探针的PCA
ctrl=read.table(paste0(Important,"GDPH.CKD_ESRD.CtrlProbe_PC.txt"),header=T,sep=" ")
PC=ctrl

##############################读入批次效应Date
Date=as.factor(as.character(covariance$Date))

###############################读入批次效应Sentrix ID
Sentrix_ID=as.factor(as.character(covariance$Sentrix_ID))

####################################读入Age和Gender
sex_age=covariance[,c(2,3)]####要改

##################################读入细胞类型
cell=read.table(paste0(Important,"GDPH.CKD_ESRD.Cell.Type.txt"),header=T,sep=" ",row.names=1)
cell=cell[,1:6]


##############################读入吸烟
Smoking=covariance$smokingScore

############################读入扫描仪
scan=covariance$scanner



##########合并
pheno=cbind(sex_age,PC)
pheno=cbind(pheno,cell)
pheno=cbind(pheno,Smoking)
pheno=cbind(pheno,scan)



######################################################Step2 使用linear model 对标准化后的mval进行校正上述斜变量,根据batch effect检测的热图，可以看出Date对标准化后的mval的影响最大
library(limma)
res_mval=removeBatchEffect(mval,covariates = pheno)
res_bvals=rnb.mval2beta(res_mval)

########################################对residual mval进行PCA分析
res_mval_t=t(res_mval)
res_mvals_pca=prcomp(res_mval_t,center=TRUE,scale=TRUE)
res_mvals_sum_pc=summary(res_mvals_pca)
res_importance=res_mvals_sum_pc$importance[,1:100]
res_mvals_pca=res_mvals_pca$x
res_mvals_pca=res_mvals_pca[,1:100]

#################################################保存最终可以用于后续差异分析的残差mval

res_PC=res_mvals_pca[,1:5]
pheno=cbind(pheno,res_PC)
new_mval=removeBatchEffect(mval[1:10000,],covariates = pheno)
res_mval5=new_mval
res_bvals5=rnb.mval2beta(res_mval5)

pheno=pheno[, 1:(ncol(pheno) - 5)]
res_PC=res_mvals_pca[,1:8]
pheno=cbind(pheno,res_PC)
new_mval=removeBatchEffect(mval[1:10000,],covariates = pheno)
res_mval8=new_mval
res_bvals8=rnb.mval2beta(res_mval8)


pheno=pheno[, 1:(ncol(pheno) - 8)]
res_PC=res_mvals_pca[,1:10]
pheno=cbind(pheno,res_PC)
new_mval=removeBatchEffect(mval[1:10000,],covariates = pheno)
res_mval10=new_mval
res_bvals10=rnb.mval2beta(res_mval10)


#####################################################提取和其它亚型共同的CpG位点 (要改）

library(data.table)
#data=fread(paste0(index,"All_merged_file.txt")) ####另外三个亚型PKD IgAN  HRD的结果

#res_mval=subset(res_mval,rownames(res_mval) %in% data$CpG)
#res_bvals=subset(res_bvals,rownames(res_mval) %in% data$CpG)


####################################保存残差数据
setwd(Important)

write.table(res_mvals_pca,file=paste0(Important,"GDPH.CKD_ESRD.residual_PCA.txt"),sep ="\t",col.names =TRUE,row.names =TRUE)
write.table(res_importance,file=paste0(Important,"GDPH.CKD_ESRD.residual_PCA_importance.txt"),sep ="\t",col.names =TRUE,row.names =TRUE)

write.table(res_mval10, file=paste0(Important,"test.residual10.txt"),row.names =TRUE,quote=F)
write.table(res_mval8, file=paste0(Important,"test.residual8.txt"),row.names =TRUE,quote=F)
write.table(res_mval5, file=paste0(Important,"test.residual5.txt"),row.names =TRUE,quote=F)





#############################################################检验不同残差t水平下的差异分析的inflation: 最后选择了10个PC，对应的inflation是1.3
library(data.table)
res_mval5=fread(paste0(Important,"test.residual5.txt"))
res_mval10=fread(paste0(Important,"test.residual10.txt"))
res_mval8=fread(paste0(Important,"test.residual8.txt"))

##############读入协变量
covariance=read.table(paste0(Important,"GDPH.CKD_ESRD.pheno.csv"),sep=",",header=T)




###########################细胞类型
cell=read.table(paste0(Important,"GDPH.CKD_ESRD.Cell.Type.txt"),sep=" ",header=T)
cell=cell[,1:6]



############################性别和年龄
Age=covariance$Age
Gender=covariance$Gender


############读入分组
Sample_group=covariance$Sample_Group


################读入吸烟
Smoking=covariance$Smoking

########################################################logistic回归分析


######################################测试residual PC 为5个
result1=c()

header=unlist(as.vector(res_mval5[,1]))
data1=as.data.frame(t(res_mval5))
data1=data1[-1,]
colnames(data1)=header
library(dplyr)
data1 <- data1 %>% mutate_all(as.numeric)

for ( i in 1:ncol(data1)) { fit = glm(Sample_group~data1[,i] +Age+Gender+Smoking+cell[,1]+cell[,2]+cell[,3]+cell[,4]+cell[,5]+cell[,6],family=binomial); result1=as.data.frame(rbind(result1,c(colnames(data1)[i],coef(summary(fit))[2,c(1,2,4)])));names(result1)=c("CpG","Estimate" ,"SE" ,"Pvalue")}

########科学计数法
result1$Pvalue=as.numeric(result1$Pvalue)
result1$SE=as.numeric(result1$SE)
result1$Estimate=as.numeric(result1$Estimate)
result1$Trait="CKD"

 p_value=result1$Pvalue
 inflation1=median(qchisq(p_value, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)

custom_format <- function(x) {
  if (abs(x) < 1e-3) {
    return(format(x, scientific = TRUE, digits = 3))
  } else if (abs(x) >= 0.001 && abs(x) < 0.005) {
    return(format(round(x, 3), nsmall = 3))
  } else {
    return(format(round(x, 2), nsmall = 2))
  }
}

result1$Pvalue <- as.numeric(sapply(result1$Pvalue, custom_format))
result1$Estimate <- as.numeric(sapply(as.numeric(result1$Estimate), custom_format))
result1$SE <- as.numeric(sapply(as.numeric(result1$SE), custom_format))
 result1=result1[order(result1$Pvalue),]
result1[1:20,]


paste0 ("the first 5 PCs of residual to correct methylation, the inflation is: ", inflation1)


######################################测试residual PC 为8个
result1=c()

header=unlist(as.vector(res_mval8[,1]))
data1=as.data.frame(t(res_mval8))
data1=data1[-1,]
colnames(data1)=header
library(dplyr)
data1 <- data1 %>% mutate_all(as.numeric)

for ( i in 1:ncol(data1)) { fit = glm(Sample_group~data1[,i] +Age+Gender+Smoking+cell[,1]+cell[,2]+cell[,3]+cell[,4]+cell[,5]+cell[,6],family=binomial); result1=as.data.frame(rbind(result1,c(colnames(data1)[i],coef(summary(fit))[2,c(1,2,4)])));names(result1)=c("CpG","Estimate" ,"SE" ,"Pvalue")}

########科学计数法
result1$Pvalue=as.numeric(result1$Pvalue)
result1$SE=as.numeric(result1$SE)
result1$Estimate=as.numeric(result1$Estimate)
result1$Trait="CKD"

 p_value=result1$Pvalue
 inflation1=median(qchisq(p_value, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)

custom_format <- function(x) {
  if (abs(x) < 1e-3) {
    return(format(x, scientific = TRUE, digits = 3))
  } else if (abs(x) >= 0.001 && abs(x) < 0.005) {
    return(format(round(x, 3), nsmall = 3))
  } else {
    return(format(round(x, 2), nsmall = 2))
  }
}

result1$Pvalue <- as.numeric(sapply(result1$Pvalue, custom_format))
result1$Estimate <- as.numeric(sapply(as.numeric(result1$Estimate), custom_format))
result1$SE <- as.numeric(sapply(as.numeric(result1$SE), custom_format))
 result1=result1[order(result1$Pvalue),]
result1[1:20,]


paste0 ("the first 8 PCs of residual to correct methylation, the inflation is: ", inflation1)


######################################测试residual PC 为10个
result1=c()

header=unlist(as.vector(res_mval10[,1]))
data1=as.data.frame(t(res_mval10))
data1=data1[-1,]
colnames(data1)=header
library(dplyr)
data1 <- data1 %>% mutate_all(as.numeric)


for ( i in 1:ncol(data1)) { fit = glm(Sample_group~data1[,i] +Age+Gender+Smoking+cell[,1]+cell[,2]+cell[,3]+cell[,4]+cell[,5]+cell[,6],family=binomial); result1=as.data.frame(rbind(result1,c(colnames(data1)[i],coef(summary(fit))[2,c(1,2,4)])));names(result1)=c("CpG","Estimate" ,"SE" ,"Pvalue")}

########科学计数法
result1$Pvalue=as.numeric(result1$Pvalue)
result1$SE=as.numeric(result1$SE)
result1$Estimate=as.numeric(result1$Estimate)
result1$Trait="CKD"

 p_value=result1$Pvalue
 inflation1=median(qchisq(p_value, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)

custom_format <- function(x) {
  if (abs(x) < 1e-3) {
    return(format(x, scientific = TRUE, digits = 3))
  } else if (abs(x) >= 0.001 && abs(x) < 0.005) {
    return(format(round(x, 3), nsmall = 3))
  } else {
    return(format(round(x, 2), nsmall = 2))
  }
}

result1$Pvalue <- as.numeric(sapply(result1$Pvalue, custom_format))
result1$Estimate <- as.numeric(sapply(as.numeric(result1$Estimate), custom_format))
result1$SE <- as.numeric(sapply(as.numeric(result1$SE), custom_format))
 result1=result1[order(result1$Pvalue),]
result1[1:20,]


paste0 ("the first 10 PCs of residual to correct methylation, the inflation is: ", inflation1)



