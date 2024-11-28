###################LD_Clumping the meQTL data (exposure): using Rstudio


library(data.table)
library(TwoSampleMR)
library(ieugwasr)

exp <- read_exposure_data(filename="C:\\Users\\84478\\Desktop\\Rstudio_file\\ESRD.1390.DMP.meQTL.rsID.txt",sep = " ",snp_col ="rsID",beta_col = "Beta", se_col = "se",effect_allele_col = "EA",other_allele_col = "OA",eaf_col = "Freq", pval_col = "Pval",phenotype_col="DNAm")

a=dplyr::tibble(rsid=exp$SNP, pval=exp$pval.exposure, id=exp$exposure)

b=ld_clump(dat=a,clump_kb = 1000,clump_r2 = 0.2,clump_p=0.05,pop = "EAS")

LD=b
names(LD)=c("SNP","pval","exposure")
LD=LD[,c(1,3)]
final=merge(LD,exp,by=c("SNP","exposure"))
exp=final

sub=subset(final,pval.exposure < 8e-5)
write.table(final,"meQTL.exposure.data.txt",row.names=F,quote=F)

