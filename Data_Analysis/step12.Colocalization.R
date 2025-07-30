
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

gwas_new=gwas %>% select ("rs_id","MAF","beta","se","p_value","N")
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

  # 提取结果并添加基因和CpG信息
  sum = as.data.frame(t(result$summary))
  sum$CpG=input$CpG[1]
  sum$block=input$block[1]
  sum$trait="eGFR"
  sum$rs_id= sentienl_rs_id


  ##该共定位 区域内meqtl和eGFR的最小P值
  sum$min_p_value_meqtl=min(input$p_value_meqtl)
  sum$min_p_value_gwas=min(input$p_value_gwas)

  #提取sentienl meQTL SNP在meQTL中的统计量

  sentienl_in_meqtl=subset(lead_meqtl, rs_id %in% sentienl_rs_id)

  #提取sentienl meQTL SNP在eGFR GWAS中的统计量 (可忽略)
  #sentienl_in_gwas=subset(gwas_new_meqtl, rs_id %in% sentienl_rs_id)
  #names(sentienl_in_gwas)=paste0(names(gwas_new_meqtl),"_gwas")

  #合并上述统计量
  all=cbind(sum,sentienl_in_meqtl)
  frame = rbind(frame, all)
}



frame_eGFR=subset(frame, (PP.H3.abf > 0.7 | PP.H4.abf > 0.7) & min_p_value_meqtl < 0.0005)   ###过滤条件


n1=length(unique(frame$CpG))  #n1=28
n2=length(unique(frame_eGFR$CpG))  #n2=10


index="/share/home/liujianjun/ESRD_Project/Part4.Bacon_Subtype_Meta_analysis/"

sentinels <- read.csv(paste0(index, "Full.Rep.DMP.Results"), sep = " ", header = TRUE)
sentinels=sentinels[,c(1,24:28)]

data_out=merge(frame_eGFR,sentinels,by="CpG")
data_out$dis=abs(data_out$MAPINFO - data_out$bp)

table="/share/home/liujianjun/ESRD_Project/Table/Output/"
write.table(frame,paste0(table,"Coloc.CpG.eGFR.txt"),row.names=F,quote=F)











####GWAS: UACR
gwas=fread(paste0(col_gwas,"UACR.GWAS.tsv"),fill=TRUE) #gwas

gwas_new=gwas %>% select ("RSID","Freq1","Effect","StdErr","P-value","n_total_sum")
names(gwas_new)=c("rs_id","MAF","beta","se","p_value","N")
gwas_new=unique(gwas_new)
if (length(which(gwas_new$beta == 0)) > 0) {gwas_new=gwas_new[-c(which(gwas_new$beta == 0)),]};
gwas_new=subset(gwas_new, rs_id %in% result_meqtl$rs_id)

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

if (length(which(is.na(input$p_value_gwas))) > 0) {input=input[-(which(is.na(input$p_value_gwas))),]};
if (length(which(is.na(input$p_value_beta))) > 0) {input=input[-(which(is.na(input$p_value_beta))),]};

 input$varbeta_meqtl=input$se_meqtl*input$se_meqtl;
  input$varbeta_gwas=input$se_gwas*input$se_gwas;

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

   # 提取结果并添加基因和CpG信息
  sum = as.data.frame(t(result$summary))
  sum$CpG=input$CpG[1]
  sum$block=input$block[1]
  sum$trait="UACR"
  min_ID=which(input$p_value_meqtl == min(input$p_value_meqtl))
  df=input[min_ID,c(1,2,3,4,5,6,7,13,14,15)]  ###提取lead meQTL的summary statistics和 trait GWAS的summary statistics
  all=cbind(sum,df)
  all_sub=subset(all,p_value_gwas %in% min(all$p_value_gwas))
  frame = rbind(frame, all_sub)
}


frame_UACR=subset(frame, (PP.H3.abf > 0.7 | PP.H4.abf > 0.7) & p_value_gwas < 0.05)   ###过滤条件





###GWAS: BUN
gwas=fread(paste0(col_gwas,"BUN.mean.GWAS.tsv"),fill=TRUE) #gwas

gwas_new=gwas %>% select ("rsid","effect_allele_frequency","beta","standard_error","p_value","n")
names(gwas_new)=c("rs_id","MAF","beta","se","p_value","N")
gwas_new=unique(gwas_new)
if (length(which(gwas_new$beta == 0)) > 0) {gwas_new=gwas_new[-c(which(gwas_new$beta == 0)),]};
gwas_new=subset(gwas_new, rs_id %in% result_meqtl$rs_id)

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

if (length(which(is.na(input$p_value_gwas))) > 0) {input=input[-(which(is.na(input$p_value_gwas))),]};
if (length(which(is.na(input$p_value_beta))) > 0) {input=input[-(which(is.na(input$p_value_beta))),]};

 input$varbeta_meqtl=input$se_meqtl*input$se_meqtl;
  input$varbeta_gwas=input$se_gwas*input$se_gwas;

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

    # 提取结果并添加基因和CpG信息
  sum = as.data.frame(t(result$summary))
  sum$CpG=input$CpG[1]
  sum$block=input$block[1]
  sum$trait="BUN"
  min_ID=which(input$p_value_meqtl == min(input$p_value_meqtl))
  df=input[min_ID,c(1,2,3,4,5,6,7,13,14,15)]  ###提取lead meQTL的summary statistics和 trait GWAS的summary statistics
  all=cbind(sum,df)
  all_sub=subset(all,p_value_gwas %in% min(all$p_value_gwas))
  frame = rbind(frame, all_sub)
}

frame_BUN=subset(frame, (PP.H3.abf > 0.7 | PP.H4.abf > 0.7) & p_value_gwas < 0.05)   ###过滤条件








###GWAS: ESRD
gwas=fread(paste0(col_gwas,"Formatted.ESRD.GWAS.tsv"),fill=TRUE) #gwas


gwas_new=gwas %>% select ("rs_id","MAF","beta","se","p_value","N")
gwas_new=unique(gwas_new)
if (length(which(gwas_new$beta == 0)) > 0) {gwas_new=gwas_new[-c(which(gwas_new$beta == 0)),]};
gwas_new=subset(gwas_new, rs_id %in% result_meqtl$rs_id)

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

if (length(which(is.na(input$p_value_gwas))) > 0) {input=input[-(which(is.na(input$p_value_gwas))),]};
if (length(which(is.na(input$p_value_beta))) > 0) {input=input[-(which(is.na(input$p_value_beta))),]};

 input$varbeta_meqtl=input$se_meqtl*input$se_meqtl;
  input$varbeta_gwas=input$se_gwas*input$se_gwas;

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

  # 提取结果并添加基因和CpG信息
  sum = as.data.frame(t(result$summary))
  frame = rbind(frame, sum)
}





#GWAS: CAD
gwas=fread(paste0(col_gwas,"Formatted.CAD.GWAS.txt"),fill=TRUE) #gwas


gwas_new=gwas %>% select ("rs_id","MAF","beta","se","p_value","N")
gwas_new=unique(gwas_new)
if (length(which(gwas_new$beta == 0)) > 0) {gwas_new=gwas_new[-c(which(gwas_new$beta == 0)),]};
gwas_new=subset(gwas_new, rs_id %in% result_meqtl$rs_id)

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

if (length(which(is.na(input$p_value_gwas))) > 0) {input=input[-(which(is.na(input$p_value_gwas))),]};
if (length(which(is.na(input$p_value_beta))) > 0) {input=input[-(which(is.na(input$p_value_beta))),]};

 input$varbeta_meqtl=input$se_meqtl*input$se_meqtl;
  input$varbeta_gwas=input$se_gwas*input$se_gwas;

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

  # 提取结果并添加基因和CpG信息
  sum = as.data.frame(t(result$summary))
  frame = rbind(frame, sum)
}















#GWAS: Anemia
gwas=fread(paste0(col_gwas,"Formatted.Hyperparathyroidism.GWAS.txt"),fill=TRUE) #gwas


gwas_new=gwas %>% select ("rs_id","MAF","beta","se","p_value","N")
gwas_new=unique(gwas_new)
if (length(which(gwas_new$beta == 0)) > 0) {gwas_new=gwas_new[-c(which(gwas_new$beta == 0)),]};
gwas_new=subset(gwas_new, rs_id %in% result_meqtl$rs_id)

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

if (length(which(is.na(input$p_value_gwas))) > 0) {input=input[-(which(is.na(input$p_value_gwas))),]};
if (length(which(is.na(input$p_value_beta))) > 0) {input=input[-(which(is.na(input$p_value_beta))),]};

 input$varbeta_meqtl=input$se_meqtl*input$se_meqtl;
  input$varbeta_gwas=input$se_gwas*input$se_gwas;

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

  # 提取结果并添加基因和CpG信息
  sum = as.data.frame(t(result$summary))
  frame = rbind(frame, sum)
}











####GWAS: UACR
gwas=fread(paste0(col_gwas,"UACR.GWAS.tsv"),fill=TRUE) #gwas

gwas_new=gwas %>% select ("RSID","Freq1","Effect","StdErr","P-value","n_total_sum")
names(gwas_new)=c("rs_id","MAF","beta","se","p_value","N")
gwas_new=unique(gwas_new)
if (length(which(gwas_new$beta == 0)) > 0) {gwas_new=gwas_new[-c(which(gwas_new$beta == 0)),]};
gwas_new=subset(gwas_new, rs_id %in% result_meqtl$rs_id)

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

if (length(which(is.na(input$p_value_gwas))) > 0) {input=input[-(which(is.na(input$p_value_gwas))),]};
if (length(which(is.na(input$p_value_beta))) > 0) {input=input[-(which(is.na(input$p_value_beta))),]};

 input$varbeta_meqtl=input$se_meqtl*input$se_meqtl;
  input$varbeta_gwas=input$se_gwas*input$se_gwas;

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

   # 提取结果并添加基因和CpG信息
  sum = as.data.frame(t(result$summary))
  sum$CpG=input$CpG[1]
  sum$block=input$block[1]
  sum$trait="UACR"
  min_ID=which(input$p_value_meqtl == min(input$p_value_meqtl))
  df=input[min_ID,c(1,2,3,4,5,6,7,13,14,15)]  ###提取lead meQTL的summary statistics和 trait GWAS的summary statistics
  all=cbind(sum,df)
  all_sub=subset(all,p_value_gwas %in% min(all$p_value_gwas))
  frame = rbind(frame, all_sub)
}


frame_UACR=subset(frame, (PP.H3.abf > 0.7 | PP.H4.abf > 0.7) & p_value_gwas < 0.05)   ###过滤条件

