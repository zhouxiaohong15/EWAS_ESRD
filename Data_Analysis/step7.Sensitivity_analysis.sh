
#########################Sensitivity analysis: Rscript

imp="/share/home/liujianjun/ESRD_Project/Important_files/"
path1="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
out="/share/home/liujianjun/ESRD_Project/Part15.Sensitivity_analysis/Output/"



######load the clinic characteristics of 998 samples
library(data.table)
pheno=fread(paste0(path1,"GDPH_998Samples.pheno.txt"))
cov=fread(paste0(imp,"GDPH.1166.Samples.Pheno.txt"))

 cov0=cov[,c(2,26,27,28,29,30,31,37,38,42,43,44,32,33,35,36)]
names(cov0)=c("Sample_ID","Albumin","TC","Triglyceride","LDLC","HDLC","Uric_Acid","HCRP","CRP","BMI","HTN","HyperL","Ca","P","Hemoglobin","PTH")
 all=merge(pheno,cov0,by="Sample_ID",all.x=TRUE)
 all <- all[match(pheno$Sample_ID, all$Sample_ID), ]

 identical(all$Sample_ID,pheno$Sample_ID)

 #write.table(all,paste0(imp,"GDPH_998Samples.pheno.txt"),quote=F,row.names=F,sep="\t")  ##相比path1目录同名的文件，多了临床基线资料



 ###load the WCC of 998 samples

 cell=read.table(paste0(path1,"GDPH_998Samples.Cell.Type.txt"),row.names=1)
 identical(rownames(cell),all$Sample_ID)


###load residuals of shared DMP in 998 samples
res=fread(paste0(index,"Shared.DMP.Residuals.txt"))

res=as.data.frame(res)
rownames(res)=res$V1
res$V1=c()

residual=as.data.frame(t(res))
residual$Sample_ID=rownames(residual)

identical(residual$Sample_ID,all$Sample_ID)


####################compare the effect size between main modle and sensitivity model



#outlier
factors="risk"
summary(as.numeric(all$risk))
var <- as.numeric(all$risk)   ##换指标分析的话 改这三行内容就行了


# 找到包含 NA 的行号
na_rows <- which(is.na(var))
# 过滤掉 NA 值
var_filter <- as.numeric(var[-na_rows])

mean <- mean(var_filter)
sd <- sd(var_filter)
outlier_threshold <- 3 * sd
outliers <- var_filter[abs(var_filter - mean) > outlier_threshold]

ID_outlier=which(var == outliers) ##找到outlier的行号
ID_NA=which(is.na(var))  ##找到NA值的行号

ID_remove= unique(c(ID_outlier,ID_NA))
ID_remove= ID_remove[order(ID_remove)] ##需要在residual和pheno文件，以及细胞类型中去除的行

cell_filter=cell[-c(ID_remove),]
all_filter=all[-c(ID_remove),]
residual_filter=residual[-c(ID_remove),]

data1=as.data.frame(table(all_filter$Sample_Group))   ##计算ESRD 和control的样本数
names(data1)[1]="Sample_Group"
data1$factors=factors

####Main model
result=c()
for ( i in 1:(ncol(residual)-1)) { 
  fit = glm(all_filter$Sample_Group~ residual_filter[,i] + all_filter$Gender+all_filter$Age+all_filter$Gender+all_filter$Smoking,family=binomial);
 result=as.data.frame(rbind(result,c(colnames(residual_filter)[i],coef(summary(fit))[2,c(1,2,4)])));
 names(result)=c("CpG","Estimate" ,"SE" ,"Pvalue")
 }


head(result)
result$Pvalue
str(result)
result$Pvalue=as.numeric(result$Pvalue)
result$SE=as.numeric(result$SE)
result$Estimate=as.numeric(result$Estimate)
result$Z_score <- result$Estimate / result$SE
result$Trait=factors

result=result[order(result$Pvalue),]
head(result)
main=result


####Sensitivity model

result=c()
for ( i in 1:(ncol(residual)-1)) { 
  fit = glm(all_filter$Sample_Group~ residual_filter[,i] + as.numeric(all_filter$risk) + all_filter$Gender+all_filter$Age+all_filter$Gender+all_filter$Smoking,family=binomial);
 result=as.data.frame(rbind(result,c(colnames(residual_filter)[i],coef(summary(fit))[2,c(1,2,4)])));
 names(result)=c("CpG","Estimate" ,"SE" ,"Pvalue")
 }


head(result)
result$Pvalue
str(result)
result$Pvalue=as.numeric(result$Pvalue)
result$SE=as.numeric(result$SE)
result$Estimate=as.numeric(result$Estimate)
result$Z_score <- result$Estimate / result$SE
result$Trait=factors

result=result[order(result$Pvalue),]
head(result)

sensitivity=result
sensitivity=sensitivity[match(main$CpG,sensitivity$CpG),]

#####Compare
final=data.frame(CpG=names(residual)[1:52],Zscore_sen=sensitivity$Z_score,Zscore_main=main$Z_score)

final$N=sum(data1$Freq)
final$case=data1$Freq[2]
final$Factor=factors
final$dir=final$Zscore_sen*final$Zscore_main

name=paste0("Sensitivity.Results.",factors,".txt")
write.table(final,paste0(out,name),row.names=F,quote=F)




###############Batch Shell 
##将上述代码复制粘贴到${code}Sisitivity.Analysis.Raw.Code.R


imp="/share/home/liujianjun/ESRD_Project/Important_files/"
path1="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
out="/share/home/liujianjun/ESRD_Project/Part15.Sensitivity_analysis/Output/"
code="/share/home/liujianjun/ESRD_Project/Part15.Sensitivity_analysis/Code/"
batch="/share/home/liujianjun/ESRD_Project/Part15.Sensitivity_analysis/Code/Batch.code/"

# 定义变量数组
variables=("Albumin" "TC" "Triglyceride" "LDLC" "HDLC" "Uric_Acid")

# 循环遍历变量数组
for i in "${variables[@]}"
do
    # 替换代码文件中的 "risk" 为当前变量值 $i，并生成新的 R 脚本文件
    sed "s/risk/$i/g" "${code}Sisitivity.Analysis.Raw.Code.R" > "${batch}Sisitivity.Analysis.$i.R"
    
    # 生成 shell 脚本文件并写入 Rscript 命令
    echo "/share/apps/R/4.4.1/bin/Rscript ${batch}Sisitivity.Analysis.$i.R" > "${batch}Sisitivity.Analysis.$i.sh"
    
    cd "${batch}"
    bsub < "${batch}Sisitivity.Analysis.$i.sh"
done




######Construct the result
out="/share/home/liujianjun/ESRD_Project/Part15.Sensitivity_analysis/Output/"

# 定义变量数组
variables=("Albumin" "TC" "Triglyceride" "LDLC" "HDLC" "Uric_Acid")

# 循环遍历变量数组
for i in "${variables[@]}"
do
    # 使用 sed 删除以 "CPG" 开头的行，并保存更改到同一文件中
    sed -i "/^CpG/d" "${out}Sensitivity.Results.$i.txt"
done

# 清空或创建空白文件 All.Sensitivity.Result.txt
echo "" > "${out}All.Sensitivity.Result.txt"

# 将所有以 Sensitivity.Results 开头的文件内容追加到 All.Sensitivity.Result.txt 中
cat "${out}Sensitivity.Results."*".txt" >> "${out}All.Sensitivity.Result.txt"

# 使用 sed 在文件开头插入标题行

sed -i "1iCpG Zscore_sen Zscore_main N case Factor dir" "${out}All.Sensitivity.Result.txt"
sed -i "2d" ${out}All.Sensitivity.Result.txt





