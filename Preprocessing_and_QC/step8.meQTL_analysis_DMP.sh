###############################################meQTL analysis



############File 1: Covariants file (Top 10 genomic PCs)

#! usr/bin/bash

hg38="/share/home/liujianjun/Reference_bigdatacompute/"
chain="/share/home/liujianjun/BGI_data/BGItrans/ESRD_ASA/PhD.paper/02.ALL.QCed.hg38.data/Important.files/"
Code="/share/home/liujianjun/BGI_data/BGItrans/ESRD_ASA/PhD.paper/01.Aligment.hg37tohg38/Batch.Code/"
picard="/share/home/liujianjun/software/picard-3.3.0/"
input="/share/home/liujianjun/BGI_data/BGItrans/ESRD_ASA/PhD.paper/01.Aligment.hg37tohg38/Raw.Aligment.Data.PLINK/"
output="/share/home/liujianjun/BGI_data/BGItrans/ESRD_ASA/PhD.paper/01.Aligment.hg37tohg38/Aligment.Data.hg38.PLINK/"
ASA_dir="/share/home/liujianjun/BGI_data/BGItrans/ESRD_ASA/PhD.paper/Impotant.files/"
QTLme="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.Methylation/"
geno="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.hg38.PLINK/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
ChinaMap="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Imputation.ChinaMAP/"
wd7="/share/home/liujianjun/GWAS.OF.ESRD/Raw.Data/3.PLINK.QCed.after.UpdatedID/"
split="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.hg38.PLINK/Split.Chromosome/"
var="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.hg38.PLINK/Genotype/"


for i in {1..22}; do 
    plink --vcf ${output}All.Samples.chr$i.QCed.hg38.vcf  --make-bed  --double-id --update-ids ${output}IID.Update.txt  --out ${output}All.UpdatedID.chr$i.QCed.hg38
done;  

echo "" > ${output}Merge.list.txt;
for i in {2..22}; do
  echo "${output}All.UpdatedID.chr$i.QCed.hg38.bed ${output}All.UpdatedID.chr$i.QCed.hg38.bim ${output}All.UpdatedID.chr$i.QCed.hg38.fam" >> ${output}Merge.list.txt
done

plink --bfile ${output}All.UpdatedID.chr1.QCed.hg38 --make-bed --merge-list ${output}Merge.list.txt   --out ${output}Merge.All.hg38

plink --bfile ${output}Merge.All.hg38 --hwe 0.000001 --mind 0.1 --maf 0.01 --geno 0.05 --make-bed --out ${output}Merge.QCed.All.hg38



##############################Extract the Samples mearsured both DNA methylation and genotype:

#!/usr/bin/Rscript

meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
ASA_dir="/share/home/liujianjun/BGI_data/BGItrans/ESRD_ASA/PhD.paper/Impotant.files/"
pheno=fread(paste0(meth,"GDPH_998Samples.pheno.csv"))
head(pheno)
pheno=as.data.frame(pheno)
file=data.frame(FID=pheno$Sample_ID,IID=pheno$Sample_ID)
head(file)
write.table(file,paste0(ASA_dir,"GDPH.Samples.list.txt"),sep="\t",quote=F,row.names=F)





###############################Do the PCA alaysis of common samples  across ASA and EPIC array (991 Samples)
#!/usr/bin/bash

plink --bfile ${wd7}Merged.QCed.ALL  --keep ${ASA_dir}GDPH.Samples.list.txt --make-bed --out ${geno}Common.Sample.hg38

plink --bfile ${geno}Common.Sample.hg38 --recode A --out ${geno}Common.Sample.hg38

awk '{$2=""; $3=""; $4=""; $5=""; $6=""; print $0}' ${geno}Common.Sample.hg38.raw > ${geno}Common.Sample.hg38.txt

sed -i "s/      / /g" ${geno}Common.Sample.hg38.txt


##########################split into 22 chrs
for chr in {1..22} ; do
  plink --bfile ${geno}Common.Sample.hg38 --chr $chr --make-bed --out ${split}Common.Sample.hg38.chr${chr}
done




####转置每条染色体基因型文件，使得行名为SNP，列名为样本名

# 生成plink命令和数据处理的通用部分

    
for i in {1..22}; do
    # 使用plink生成转换后的文件
    plink --bfile ${split}Common.Sample.hg38.chr$i --recode A-transpose --out ${split}Common.Sample.hg38.chr$i
done

for i in {1..22}; do # 删除不需要的列
    awk '{$1=""; $3=""; $4=""; $5=""; $6=""; print $0}' ${split}Common.Sample.hg38.chr${i}.traw > ${split}Common.Sample.hg38.chr${i}.txt
done;


for i in {1..22};do
sed  "s/     / /g" ${split}Common.Sample.hg38.chr${i}.txt > ${var}chr${i}.Genotype.txt
sed -i "s/SNP //g" ${var}chr${i}.Genotype.txt
sed -i '1s/_[a-zA-Z0-9]\+ / /g' ${var}chr${i}.Genotype.txt
sed -i '1s/_[a-zA-Z0-9]\+$//g' ${var}chr${i}.Genotype.txt
done

####因为没有去掉首字母空格而报错了
for i in {1..22};do
sed -i "s/^ //g" ${var}chr${i}.Genotype.txt
done


##################生成每条染色体的Location文件
# 遍历1到22染色体
for i in {1..22}; do
    # 使用awk生成染色体的位置信息
    awk 'BEGIN {print "snpid\tchr\tpos"} 
         {printf "%s\tchr%s\t%s\n", $2, $1, $4}' ${split}Common.Sample.hg38.chr${i}.bim > ${var}Genotype.hg38.chr${i}.Location
done


####################确保基因型文件和位置文件的SNP一致
echo "" ${var}Checking.SNPlist
for i in {1..22};do
# 提取文件A的第一列并去掉首个字符
A_first_col=$(cut -f1 ${var}Genotype.hg38.chr${i}.Location| sed '1d')

# 提取文件B的第一行
B_first_row=$(cut -f1  ${var}chr${i}.Genotype.txt )

# 比较两个结果
echo "" > ${var}Checking.SNPlist
if [[ "$A_first_col" == "$B_first_row" ]]; then
    echo "chr$i  same" >> ${var}Checking.SNPlist
else
    echo "chr$i  different" >> ${var}Checking.SNPlist
fi
done;


############提取最终版基因型文件的Sample List, 用于保证methylation residual的sample id 顺序与之保持一致
i=22
awk '{print $1}'  ${var}chr${i}.Genotype.txt > ${var}Correct.Order.Sample.List




#####do the PCA analysis

#!/usr/bin/bash
plink --bfile ${geno}Common.Sample.hg38 --indep-pairwise 1500 150 0.2  --maf 0.01 --make-bed  --out ${geno}Common.Sample.hg38.Pre_Prune 
plink --bfile ${geno}Common.Sample.hg38.Pre_Prune   --extract ${geno}Common.Sample.hg38.Pre_Prune.prune.in  --make-bed --out ${geno}Common.Sample.hg38.Pruned
plink --bfile  ${geno}Common.Sample.hg38.Pruned --pca 10 --out  ${geno}Common.Sample.PCA

awk -F ' ' '{print $1,$2}' ${geno}Common.Sample.hg38.fam  >  ${geno}Common.Sample.list.txt
sed -i '1iFID IID' ${geno}Common.Sample.list.txt



########Methylation Residual file

#!/usr/bin/Rscript

index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
Important="/share/home/liujianjun/ESRD_Project/Important_files/"
Comp="/share/home/liujianjun/ESRD_Project/Part5Complications/Complication_Results/"
QTLme="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/Common.Sample.Methylation/"
geno="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.hg38.PLINK/"
ASA_dir="/share/home/liujianjun/BGI_data/BGItrans/ESRD_ASA/PhD.paper/Impotant.files/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
ChinaMap="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Imputation.ChinaMAP/"
QTLme2="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.Methylation/"
var="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.hg38.PLINK/Genotype/"


library(data.table)
data=read.csv(paste0(index,"Full.Rep.DMP.Results"),sep=" ")
head(data)
CpG=data$CpG
write.table(CpG,paste0(meth,"shared.DMP.list"),quote=F,row.names=F)

#############File1: Methylation residuals: extract residuals of shared DMPs
#!/usr/bin/Rscript
library(data.table)
library(limma)


res=as.data.frame(res)
PC=fread(paste0(geno,"Common.Sample.PCA.eigenvec"))  #### 10 Genomic PCs

head(PC)
PC$V1=c()
head(PC)
names(PC)[2:11]=paste0("PC",1:10)
names(PC)[1]="Sample_ID"
write.table(PC,paste0(geno,"Common.Sample.Top10.PCs.txt"),row.names=F,quote=F)

res=fread(paste0(meth,"GDPH_998Samples.residual_mvals.txt"))  ##methylation
sample= as.data.frame(fread(paste0(geno,"Common.Sample.list.txt")))  ###Common sample list between EPIC and ASA aray

res=as.data.frame(res)
rownames(res)=res$V1
res$V1=c()
ID=which(colnames(res) %in% sample$IID)
identical(colnames(res)[ID],sample$IID)
res_sub=res[,ID]

res_sorted <- res_sub[, match(sample$IID, colnames(res_sub))]  ###使得methylation 和genotype对应的Sample_ID 的顺序保持一致
identical(colnames(res_sorted),sample$IID)

covariance=read.table(paste0(meth,"GDPH_998Samples.pheno.csv"),sep=",",header=T)
# 确保 Sample_ID 和 IID 的数据类型一致（假设它们应该是字符型）
covariance$Sample_ID = as.character(covariance$Sample_ID)
sample$IID = as.character(sample$IID)
ID=which(covariance$Sample_ID %in% sample$IID)
cov=covariance[ID,]
cov=cov[order(match(cov$Sample_ID, sample$IID)),]
identical(cov$Sample_ID,sample$IID)  ###要确认是否为TRUE

library(limma)
PCA=PC[,c(2:11)]
Sample_Group=cov$Sample_Group
counf=cbind(PCA,Sample_Group)
res_mval_com=removeBatchEffect(res_sorted,covariates = counf)
res_mval_com=as.data.frame(res_mval_com)

rownames(res_mval_com)=rownames(res)

# 确保和检查residual 的sample ID 和 最终版本的genotype 文件的sample ID顺序保持一致：很重要
list=fread(paste0(var,"Correct.Order.Sample.List"))
names(list)="Sample_ID"
identical(colnames(res_mval_com),list$Sample_ID) #####注意辨别是否为TRUE

write.table(res_mval_com,file=paste0(QTLme,"Common.Samples.residual"),row.names =TRUE,quote=F)  #####residuals correcting batch effect, 10 PCs and ESRD of common samples
write.table(res_mval_com,file=paste0(QTLme2,"Common.Samples.residual"),row.names =TRUE,quote=F)  #####residuals correcting batch effect, 10 PCs and ESRD of common samples



#####################Extract residual of ESRD shared DMPs
res_mval_com=as.data.frame(fread(paste0(QTLme,"Common.Samples.residual")))
rownames(res_mval_com)=res_mval_com$V1
res_mval_com$V1=c()
DMP=read.table(paste0(meth,"shared.DMP.list"),header=T)
CpG=DMP$x
sub=subset(res_mval_com,rownames(res_mval_com) %in% CpG )

list=fread(paste0(var,"Correct.Order.Sample.List"))
names(list)="Sample_ID"
identical(colnames(sub),list$Sample_ID)

write.table(sub,paste0(QTLme,"shared.DMP.Residuals"),row.names=T,quote=F)
write.table(sub,paste0(QTLme2,"shared.DMP.Residuals"),row.names=T,quote=F)



###########File 2 Shared DMP location
index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
Important="/share/home/liujianjun/ESRD_Project/Important_files/"
Comp="/share/home/liujianjun/ESRD_Project/Part5.Complications/Complication_Results/"
meth="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/"
geno="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.hg38.PLINK/"
ASA_dir="/share/home/liujianjun/BGI_data/BGItrans/ESRD_ASA/PhD.paper/Impotant.files/"
QTLme="/share/home/liujianjun/ESRD_Project/Part01.GDPH_998Samples/Important_files/Common.Sample.Methylation/"
QTLme2="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.Methylation/"
ChinaMap="/share/home/liujianjun/BGI_data/BGItrans/ESRD_ASA/PhD.paper/01.Aligment.hg37tohg38/Aligment.Data.hg38.PLINK/"
Imp="/share/home/liujianjun/ESRD_Project/Important_files/"
geno_out="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Imputation.ChinaMAP/"
var="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Common.Sample.hg38.PLINK/Genotype/"



library(data.table)
DMP=read.csv(paste0(index,"Full.Rep.DMP.Results"),sep=" ")


near=fread(paste0(Imp,"Slim.MethylationEPIC_v-1-0_B4.txt"))
names(near)[1]="CpG"

sub=subset(near,CpG %in% DMP$CpG)
DMP=sub[,1:3]

names(DMP)[1]="geneid"
DMP$chr=paste0("chr",DMP$CHR)
DMP$left=DMP$MAPINFO
DMP$right=DMP$MAPINFO
DMP$MAPINFO=c()
DMP$CHR=c()

write.table(DMP,paste0(QTLme2,"shared.DMP.hg19.Coloc"),sep="\t",row.names=F,quote=F)


######################hg19 to hg38
http://genome.ucsc.edu/cgi-bin/hgLiftOver







#######################################File 2: Genotype Files
#!usr/bin/bsh


for i in {1..22}; do
plink --bfile ${ChinaMap}All.UpdatedID.chr$i.QCed.hg38 --recode A --out ${geno_out}Chr$i.genotype
done;


#######################################File 4: Genotype  Location
for i in {1..22}; do
awk 'BEGIN {print "snpid\tchr\tpos"} 
     {printf "%s\tchr%s\t%s\n", $2, $1, $4}' ${ChinaMap}All.UpdatedID.chr$i.QCed.hg38.bim  > ${geno_out}Chr$i.Location
done;




###################################MatrixEQTL analysis (Reference:https://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/R.html#cis)
#!usr/bin/Rscript
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



library("MatrixEQTL")
library(data.table)

# Genotype file name
SNP_file_name = paste0(var, "chr$i.Genotype.txt");
snps_location_file_name = paste0(var, "Genotype.hg38.chr$i.Location");
# Gene expression file name
expression_file_name = paste0(QTLme2, "shared.DMP.Residuals");
gene_location_file_name = paste0(QTLme2, "shared.DMP.hg38.Coloc");


# Covariates file name
# Set to character() for no covariates



# Output file name
output_file_name_cis = paste0(obtain,"cismeQTL.chr$i.txt");
output_file_name_tra = paste0(obtain,"transmeQTL.chr$i.txt");

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 10e-5;
pvOutputThreshold_tra = 10e-5;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");
# Distance for local gene-SNP pairs

cisDist = 1e6;


## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = " "; # the blank character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 2000; # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);


## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = " "; # the blank character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1; # one row of column labels
gene$fileSkipColumns = 1; # one column of row labels
gene$fileSliceSize = 2000; # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);


####Position
snpspos = read.table(snps_location_file_name, header = TRUE, sep="\t" , stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, sep="\t", stringsAsFactors = FALSE);

###Checking the samples of snp and methylation file
snps
gene

########Run
useModel = modelLINEAR

me = Matrix_eQTL_main(
snps = snps,
gene = gene,
output_file_name = output_file_name_tra,
pvOutputThreshold = pvOutputThreshold_tra,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos,
genepos = genepos,
cisDist = cisDist,
pvalue.hist = TRUE,
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);


#unlink(output_file_name_tra);
#unlink(output_file_name_cis);


## Results:
#cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
#cat('Detected local eQTLs:', '\n');
#show(me$cis$eqtls)
#cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)
## Make the histogram of local and distant p-values
#plot(me)
#dev.off()


#########Batch R code
code="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Code/"
batch="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Code/Batch_Code/"
for i in {1..22};do
sed "s/\$i/$i/g" ${code}step1.meQTL.R > ${batch}meQTL.chr$i.R
echo "/share/apps/R/4.4.1/bin/Rscript ${batch}meQTL.chr$i.R" > ${batch}meQTL.chr$i.sh
bsub -e ${batch}error.chr$i.log  -o ${batch}output.chr$i.log  < ${batch}meQTL.chr$i.sh
done;




########Construct the result of whole cis and trans meQTL results

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





# cis
echo "" > ${obtain}cismeQTL.All.Results.txt
cat ${obtain}cismeQTL.chr*.txt  >> ${obtain}cismeQTL.All.Results.txt
sed -i '/^$/d' ${obtain}cismeQTL.All.Results.txt
head -1 ${obtain}cismeQTL.All.Results.txt > ${obtain}header.txt
sed -i '/^SNP/d' ${obtain}cismeQTL.All.Results.txt
cat ${obtain}header.txt ${obtain}cismeQTL.All.Results.txt > ${obtain}temp && mv ${obtain}temp ${obtain}cismeQTL.All.Results.txt


# trans
echo "" > ${obtain}transmeQTL.All.Results.txt
cat ${obtain}transmeQTL.chr*.txt  >> ${obtain}transmeQTL.All.Results.txt
sed -i '/^$/d' ${obtain}transmeQTL.All.Results.txt
head -1 ${obtain}transmeQTL.All.Results.txt > ${obtain}header.txt
cat ${obtain}header.txt
sed -i '/^SNP/d' ${obtain}transmeQTL.All.Results.txt
cat ${obtain}header.txt ${obtain}transmeQTL.All.Results.txt > ${obtain}temp && mv ${obtain}temp ${obtain}transmeQTL.All.Results.txt








#########################Cojo conditional analysis

## Step1: input file format:
#SNP A1 A2 freq b se p N
#chrX_2781514_C_A C A 0.57 0.02 0.043 0.7
#chrX_2781927_A_G G A 0.64 0.03 0.046 0.43

#!usr/bin/Rscript


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

library(data.table)
freq=fread(paste0(geno,"Common.Sample.hg38.frq"))  #frequency of the SNPs

cis=fread(paste0(obtain,"cismeQTL.All.Results.txt"))
trans=fread(paste0(obtain,"transmeQTL.All.Results.txt"))
cis$Group="cis"
trans$Group="trans"


qtl=rbind(cis,trans)
sub=subset(freq, SNP %in% qtl$SNP)
all=merge(qtl,sub,by="SNP")

write.table(all,paste0(obtain,"Full.meQTL.Results.txt"),row.names=F,quote=F)

#Extract signifcant meQTL pairs, and add the information of position of meQTL-SNP and meQTL-CpG
library(data.table)
full=fread(paste0(obtain,"Full.meQTL.Results.txt"))

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

 ## Annotate Nearest Genes to meQTL-SNPs
 library(ChIPseeker)
library(org.Hs.eg.db)
library(GenomicFeatures)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library('clusterProfiler')


library(dplyr)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

# Generate input file for chipseeker: note: remove the header of bed files: # SNP CHR BP BP (the seperate must be tab instead of blank)
bed=sig[,c(2,8,13)] # SNP CHR BP BP
bed1=unique(bed)
bed1$CHR=paste0("chr",bed1$SNP.chr)
final=bed1[,c(4,3,3,1)]
write.table(final,paste0(obtain,"Significant.meQTL.SNP.bed"),row.names=F,quote=F,col.names=F,sep="\t")

peak <- readPeakFile(paste0(obtain,"Significant.meQTL.SNP.bed"))
peakAnno <- annotatePeak(peak, tssRegion = c(-3000, 3000), TxDb = txdb)

 peakAnno 
str(peakAnno)
as.GRanges(peakAnno)
final=as.data.frame(as.GRanges(peakAnno))

#######################将基因ID转换为基因symbol
library(org.Hs.eg.db)
library('clusterProfiler')

gene<- bitr(final$geneId, fromType = "ENTREZID", toType=c("SYMBOL"),OrgDb = org.Hs.eg.db)
head(gene)
 names(gene)[1]="geneId"

result=merge(final,gene,by="geneId") # 基因注释结果
 names(result)[7]="SNP"
head(sig)  ### [1] "gene"    "SNP"     "beta"    "t-stat"  "P"       "FDR"     "Group"  "SNP.chr" "A1"      "A2"      "MAF"     "NCHROBS" "SNP_pos" "CpG.chr"  "CpG_pos"
 frame=merge(sig,result,by="SNP")

 write.table(frame,paste0(obtain,"meQTL_SNPs_table2.txt"),row.names=F,sep="\t")



##Construct the meQTL SNP table
library(data.table)

qtl <- fread(paste0(obtain, "meQTL_SNPs_table2.txt"))

### cis
cis=subset(qtl,Group == "cis")
cis$Distance.SNP.CpG = abs(cis$SNP_pos -cis$CpG_pos)

# 遍历每个独特的基因
result <- data.frame()
for (i in unique(cis$gene)) {
  sub <- subset(cis, gene == i)  # 删除重复的 `sub <- subset(cis, gene == i)`
  
  nSNPs <- nrow(sub)
  tenKb <- sum(sub$Distance.SNP.CpG < 10000)
  hundredKb <- sum(sub$Distance.SNP.CpG < 100000)
  thousandKb <- sum(sub$Distance.SNP.CpG < 1000000)
  same_chr <- nSNPs  # 直接等于 nSNPs
  

  
  # 创建数据框
  data <- data.frame(
    CpG = i,
    chr = unique(sub$CpG.chr),
    BP = unique(sub$CpG_pos),
    nSNPs = nSNPs,
    tenKb = tenKb,
    hundredKb = hundredKb,
    thousandKb = thousandKb,
    same_chr = same_chr  # 修复拼写错误
  )
  
  # 追加数据到 result
  result <- rbind(result, data)
}

# 显示前几行
head(result)





# trans
trans=subset(qtl, Group == "trans")
n=as.data.frame(unlist(table(trans$gene)))
head(n)
names(n)[1]="CpG"
names(n)[2]="n_tansSNP"


######merge cis and trans

index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"
sentinels=read.csv(paste0(index,"Full.Rep.DMP.Results"),header=T,sep=" ")   # #add the nearest gene for CpGs
head(sentinels)
names(sentinels)
near=sentinels[,c(1,24)]


file_trans=n
file_cis=result 

summary=merge(file_cis,file_trans,by="CpG",all= TRUE)
out=merge(summary,near,by="CpG",all.x=TRUE)

slim=fread(paste0(Imp,"Slim.MethylationEPIC_v-1-0_B4.txt"))
slim=slim[,c(1,5)]
names(slim)[1]="CpG"

output=merge(out,slim,by="CpG",all.x=TRUE)


write.table(output,paste0(obtain,"meQTL_SNP_table1.txt"),row.names=F,quote=F)


# Format file all in to input files of cojo:
#SNP A1 A2 freq b se p N
#chrX_2781514_C_A C A 0.57 0.02 0.043 0.7
#chrX_2781927_A_G G A 0.64 0.03 0.046 0.43

 slim=all[,c(1,9,10,11,3,4,5)]
 names(slim)[7]="p"
 names(slim)[6]="se"
 names(slim)[4]="freq"
 names(slim)[5]="b"
 
  slim$se <- sqrt(((slim$b)^2)/qchisq(slim$p, 1, lower.tail = F))
  slim$N=all$NCHROBS/2
  slim$CpG=all$gene
  slim=subset(slim,CpG %in% sig$gene)

  write.table(slim,paste0(obtain,"meQTL.COJOInput.txt"),row.names=F,quote=F)

# Generate meQTLs by CpG to batch running cojo analysis
for (i in unique(slim$CpG)) {
  frame <- subset(slim, CpG == i)  # Use 'i' directly, not '$i'
  frame$CpG <- NULL  # Remove the CpG column
  name <- i
  write.table(frame, paste0(class, name), row.names = FALSE, quote = FALSE)
}





####Step2 : COJO analysis : conditional analysis to find independent SNPs (batch)
#!usr/bin/bash

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

indep="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/meQTL.Results/Conditional.Results/" 
con="/share/home/liujianjun/software/gcta-1.94.3-linux-kernel-3-x86_64/" ##GCTA软件的路径


code="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Code/"
batch="/share/home/liujianjun/ESRD_Project/Part10.meQTL/Important_files/Code/Batch_Code/"

for file in $(find ${class} -type f -name 'cg*'); do
  name=$(basename "$file")
 echo "${con}gcta64 --bfile ${geno}Common.Sample.hg38 --maf 0.01 --cojo-file ${file}  --cojo-p 10e-8  --cojo-wind 10000  --cojo-slct --out ${indep}Independent.${name}.txt" > ${batch}cojo.${name}.sh
 bsub -e ${batch}error.${name}.log   -o ${batch}output.${name}.log  < ${batch}cojo.${name}.sh
done;


for file in ${indep}Independent.cg*.txt.jma.cojo; do
    # 使用 sed 去掉 "Independent" 和 ".txt.jma.cojo"
    filename=$(basename "$file" .txt.jma.cojo)
    filename=$(echo "$filename" | sed 's/^Independent\.//')  # 去掉 "Independent" 前缀

    # 使用 awk 在文件的最后一列添加文件名，并将列名设为 CpG
    awk -v filename="$filename" 'BEGIN {OFS="\t"} NR == 1 {print $0, "CpG"} NR > 1 {print $0, filename}' "$file" > ${indep}temp_file && mv ${indep}temp_file "$file"
done

echo "" > ${indep}All.Lead.meQTL.SNP.txt
echo "" > ${indep}temp_file
cat ${indep}Independent.cg*.txt.jma.cojo >>  ${indep}All.Lead.meQTL.SNP.txt
sed -i '/^$/d'  ${indep}All.Lead.meQTL.SNP.txt
head -1 ${indep}All.Lead.meQTL.SNP.txt > ${indep}header.txt
sed -i '/^Chr/d' ${indep}All.Lead.meQTL.SNP.txt
sed -i '/^$/d'  ${indep}All.Lead.meQTL.SNP.txt
cat ${indep}header.txt  ${indep}All.Lead.meQTL.SNP.txt >> ${indep}temp_file && mv ${indep}temp_file ${indep}All.Lead.meQTL.SNP.txt 
sed -i '/^$/d'  ${indep}All.Lead.meQTL.SNP.txt







