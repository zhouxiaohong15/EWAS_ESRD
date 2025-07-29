#############################################Case CKD-ESRD, Control CKD
###############读入residual mval
row$i=read.table("${path1}residual_mvals_block$i.txt",header=F,sep=" ")
row.names=row$i[,1]
data$i=t(row$i[,-1])
colnames(data$i)=row.names
file=row$i


##############读入协变量
covariance=read.table("${cov}GDPH.CKD_ESRD.pheno.txt",sep=" ",header=T)



###########################细胞类型 #######要改
cell=read.table("${cov}GDPH.CKD_ESRD.Cell.Type.txt",sep=" ",header=T)

cell=cell[,1:6]



############################性别和年龄
Age=covariance$Age
Gender=covariance$Gender


############读入分组
Sample_group=covariance$Sample_Group


################读入吸烟
Smoking=covariance$Smoking

########################################################logistic回归分析
result=c()

for ( i in 1:nrow(file)) { fit = glm(Sample_group~data$i[,i] +Age+Gender+Smoking+cell[,1]+cell[,2]+cell[,3]+cell[,4]+cell[,5]+cell[,6],family=binomial); result=as.data.frame(rbind(result,c(colnames(data$i)[i],coef(summary(fit))[2,c(1,2,4)])));names(result)=c("CpG","Estimate" ,"SE" ,"Pvalue")}

########科学计数法
result$Pvalue=as.numeric(result$Pvalue)
result$SE=as.numeric(result$SE)
result$Estimate=as.numeric(result$Estimate)
result$Trait="CKD"

#custom_format <- function(x) {
#  if (abs(x) < 1e-3) {
#    return(format(x, scientific = TRUE, digits = 3))
#  } else if (abs(x) >= 0.001 && abs(x) < 0.005) {
#    return(format(round(x, 3), nsmall = 3))
#  } else {
#    return(format(round(x, 2), nsmall = 2))
#  }
#}

#result$Pvalue <- as.numeric(sapply(result$Pvalue, custom_format))
#result$Estimate <- as.numeric(sapply(as.numeric(result$Estimate), custom_format))
#result$SE <- as.numeric(sapply(as.numeric(result$SE), custom_format))
#result=result[order(result$Pvalue),]
#head(result)


write.table(result,file="${CKD2}residual_result_block$i",sep = " ", row.names = F,quote=F)





