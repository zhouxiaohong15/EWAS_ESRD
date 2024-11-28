#!/usr/bin/env Rscript

# Load necessary libraries
library(ff)
library(RnBeads.hg19)
library(methylumi)
library(foreach)
library(doParallel)
library(RnBeads)
library(grid)


##############################################################Step1. Normalization
# Set temporary directory for large data files
options(fftempdir="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tmp")

# RnBeads configuration options
rnb.options(disk.dump.big.matrices=TRUE, disk.dump.bigff=TRUE)
rnb.options(logging.disk=TRUE)
rnb.options(enforce.memory.management=TRUE)
rnb.options(identifiers.column="Sample_ID", import.idat.platform='probesEPIC', import.sex.prediction=TRUE)

# Define file paths
analysis.dir = "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartD.ESRD.GDPH_ESRD_CKD/01.Raw_to_QC/01.import/"
idat.dir = "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/Part.Raw_Data/"
sample.annotation = "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartD.GDPH_ESRD_CKD/01.Raw_to_QC/Important_files/GDPH.DKD_ESRD.pheno.csv"
report.dir = file.path(analysis.dir, "import_reports")

# Step 1: Import data (IDAT files) into RnBeads object
data.source <- c(idat.dir, sample.annotation)
rnb.set.disk <- rnb.execute.import(data.source = data.source, data.type = "idat.dir")
rnb.set = rnb.set.disk

# Save the imported data
save.rnb.set(rnb.set, paste0(analysis.dir, "Import_rnb.set"), archive = FALSE)

# Step 2: Preprocessing - Load the imported data and set QC parameters
report.dir = file.path(analysis.dir, "Norm_reports")
rnb.set = load.rnb.set(paste0(analysis.dir, "Import_rnb.set"))

# Set filtering and normalization parameters
rnb.options(
  filtering.sex.chromosomes.removal = TRUE,
  identifiers.column = "Sample_ID",
  normalization.method = "bmiq",
  filtering.cross.reactive = TRUE,
  filtering.snp = "3",
  filtering.greedycut.pvalue.threshold = 0.01,
  filtering.missing.value.quantile = 0.5
)

# Preprocess data and run normalization
rnb.set.unfiltered <- rnb.set
result <- rnb.run.preprocessing(rnb.set.unfiltered, dir.reports = report.dir)
rnb.set <- result$rnb.set

# Step 3: Impute missing values using KNN
rnb.set_new = rnb.execute.imputation(rnb.set, method = "knn", update.ff = TRUE)

# Save the normalized data
save.rnb.set(rnb.set_new, paste0(analysis.dir, "02.rnb.set_Normalized"), archive = FALSE)

# Extract and save M-values and B-values
mval = mval(rnb.set, row.names = TRUE)
bval = meth(rnb.set_new, row.names = TRUE)

# Display first few values for M-values and B-values
mval[1:20, 1:5]
bval[1:20, 1:5]

# Write the normalized data to files
write.table(mval, paste0(analysis.dir, "Rnbset01.Normalized_Mval.txt"))
write.table(bval, paste0(analysis.dir, "Rnbset01.Normalized_Bval.txt"))






#########################################################Step2.identifying confounding factors affecting DNA methylation levels
#############confounders: including 1)the first 30 control probes PCs (batch effect), 2) Age, 3) Gender, 4)estimate cell propotion, 5) scanner and 5)smoking


#####################################1) Obtain the first 30 control probes PCs
#!/usr/bin/env Rscript

# Load necessary libraries
library(ff)
library(RnBeads.hg19)
library(methylumi)
library(foreach)
library(doParallel)
library(EpiDISH)
library(data.table)
library(grid)

# Set temporary directory for large data files
getOption("fftempdir")

# Configure RnBeads options
rnb.options(disk.dump.big.matrices=TRUE, disk.dump.bigff=TRUE)
rnb.options(logging.disk=TRUE)
rnb.options(enforce.memory.management=TRUE)
rnb.options(identifiers.column="Sample_ID", import.idat.platform='probesEPIC', import.sex.prediction=TRUE)

# Load phenotype and data files
pheno = read.table("GDPH.DKD_ESRD.pheno.txt", header=TRUE)
anno = read.table("850kcontrolprobe_annotation.txt", sep=" ", header=TRUE)

# Step 1: Extract signal for control probes (Cy3 and Cy5 channels)
rnb.set = load.rnb.set("path_to_normalized_data")

# Extract quality control probe signals
green = qc(rnb.set)$Cy3
red = qc(rnb.set)$Cy5

# Save control probe signals to files
write.table(green, "Rnbset01.gre_CtrlProbe.txt", sep=" ")
write.table(red, "Rnbset01.red_CtrlProbe.txt", sep=" ")

# Read back the control probe signals
green = read.table("Rnbset01.gre_CtrlProbe.txt", sep=" ", header=TRUE)
red = read.table("Rnbset01.red_CtrlProbe.txt", sep=" ", header=TRUE)

# Step 2: Filter out negative probe sites
red_norm = red[-c(49:53, 54:464),]
gre_norm = green[-c(49:53, 54:464),]
anno = anno[-c(49:53, 54:464),]

# Extract negative probe sites
red_neg = red[c(54:464),]
gre_neg = green[c(54:464),]
negative = rbind(red_neg, gre_neg)
colnames(negative) = pheno$Sample_ID

# Combine the 219 control probes and their signals
newred = cbind(anno, red_norm)
newgre = cbind(anno, gre_norm)

# Filter probes for red and green channels based on evaluation criteria
my_red = subset(newred, Evaluate.Red == "+")
my_green = subset(newgre, Evaluate.Green == "+")

# Combine the valid probes
all = rbind(my_green, my_red)
all0 = all[, -c(1:10)]  # Exclude first 10 columns

# Step 3: Perform PCA on the merged control probes (excluding negative probes)
intensity = t(all0)
intensity_new = scale(intensity, center=FALSE, scale=TRUE)
pca = prcomp(intensity_new)

# Summarize PCA results and extract the first 30 components
summary(pca)
pca_new = pca$x[, 1:30]

# Step 4: Perform PCA on negative probes (for comparison)
neg_intensity = t(negative)
neg_intensity_new = scale(neg_intensity, center=FALSE, scale=TRUE)
neg_pca = prcomp(neg_intensity_new)

# Summarize PCA results for negative probes and extract the first 30 components
summary(neg_pca)
neg_pca_new = neg_pca$x[, 1:30]

# Save PCA results for control and negative probes
write.table(pca_new, "GDPH.DKD_ESRD.CtrlProbe_PC.txt")
write.table(neg_pca_new, "GDPH.DKD_ESRD.NegCtrlProbe_PC.txt")



#########2)Obtain estimation of Cell type : using DNA methylation of specific CpGs adjusted for scanner and 30 control probe PCs
data(centDHSbloodDMC.m)

library(EpiDISH)
library(data.table)
library(limma)
library(EpiDISH)

# Load phenotype data and normalized M-values
pheno = read.table("GDPH.DKD_ESRD.pheno.txt", header = TRUE)
mval = as.matrix(fread("GDPH.DKD_ESRD.Normalized_Mval.txt"), rownames = 1)

# Extract batch effect information
scan = pheno$scanner

# Load PCA of control probes
ctrl = read.table("GDPH.DKD_ESRD.CtrlProbe_PC.txt", header = TRUE, sep = " ")
ctrl = ctrl[, 1:30]  # Use first 30 PCs

# Combine batch effect covariates
cov = cbind(ctrl, scan)

# Remove batch effect using limma
res_mval = removeBatchEffect(mval, covariates = cov)
res_bvals = rnb.mval2beta(res_mval)

# Save the corrected beta values
write.table(res_bvals, "Cell.Estimate.Residual.Bvals.txt")

# Perform cell type estimation using corrected beta values
out.l = epidish(beta.m = res_bvals, ref.m = centDHSbloodDMC.m, method = "RPC")
cell_type = out.l$estF

# Save the cell type composition estimates
write.table(cell_type, "Rnbset01.Cell_Type.txt")
write.table(cell_type, "GDPH.DKD_ESRD.Cell.Type.txt")

#########################################3) Obtain the smoking score
library(data.table)
library(limma)

# Set working directory and load phenotype data and normalized M-values
setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/01.Raw_to_QC/Important_files")
pheno = read.table("GDPH.A1.pheno.txt", header = TRUE)
mval = as.matrix(fread("GDPH.A1.Normalized_Mval.txt"), rownames = 1)

# Extract batch effect information
scan = pheno$scanner

# Load list of 1000 smoking-related CpGs
smoking = read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA.GDPH_SY_EWAS/1.Raw_QC_Normalized/01.Raw_to_QC/Important_files/Smoking_CpGs_EpiSmokEr.txt", header = FALSE)
Smoke.probe = smoking[, 1]

# Extract M-values for the 1000 smoking-related CpGs
CpGs1000 = subset(mval, rownames(mval) %in% Smoke.probe)

# Load PCA of control probes
setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/01.Raw_to_QC/Important_files")
ctrl = read.table("GDPH.A1.CtrlProbe_PC.txt", header = TRUE, sep = " ")
ctrl = ctrl[, 1:30]  # Use first 30 PCs
PC = ctrl

# Combine PCA and batch effect information
pheno = cbind(PC, scan)

# Load Sentrix Position information
Sentrix_Position = as.factor(as.character(covariance$Sentrix_Position))

# Remove batch effects from the 1000 smoking CpGs
res_mval = removeBatchEffect(CpGs1000, covariates = pheno)
res_bvals = rnb.mval2beta(res_mval)

# Save residual beta values
write.table(res_bvals, "Smoking.1000CpGs.Residual.Bvals.txt")

# Apply EpiSmokEr to predict smoking status and scores
suppressPackageStartupMessages({
  library(EpiSmokEr)
  library(IlluminaHumanMethylation450kmanifest)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(minfi)
  library(htmlTable)
  library(rmarkdown)
})

# Read phenotype and residual beta values in RStudio
pheno = read.table("C:\\Users\\84478\\Desktop\\Rstudio_file\\GDPH.A1.pheno.txt", header = TRUE)
pheno$sex = pheno$Gender
rownames(pheno) = pheno$Sample_ID

bval = read.table("C:\\Users\\84478\\Desktop\\Rstudio_file\\Smoking.1000CpGs.Residual.Bvals.txt", row.names = 1, header = TRUE)

# Predict Smoking Score (SSc)
result = epismoker(dataset = bval, samplesheet = pheno, method = "SSc")
pheno$smokingScore = result$smokingScore

# Predict Smoking Status (SSt)
result = epismoker(dataset = bval, samplesheet = pheno, method = "SSt")
pheno$PredictedSmokingStatus = result$PredictedSmokingStatus

# Save updated phenotype file with smoking scores
write.table(pheno, "C:\\Users\\84478\\Desktop\\GDPH.A1.pheno.txt", row.names = FALSE)

# Transfer phenotype file with Sample_ID and Smoking columns to cluster
smoke = read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA.GDPH_SY_EWAS/1.Raw_QC_Normalized/01.Raw_to_QC/Important_files/Smoking_status.txt", sep = "\t", header = TRUE)
covariance = read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA.GDPH_SY_EWAS/1.Raw_QC_Normalized/01.Raw_to_QC/Important_files/GDPH.pheno.txt", sep = " ", header = TRUE)

# Check if Sample_IDs match between phenotype files
identical(covariance$Sample_ID, smoke$Sample_ID)

# Add Smoking status to the covariance data
covariance$Smoking = smoke[, 2]

# Save updated phenotype data with Smoking status
write.table(covariance, "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA.GDPH_SY_EWAS/1.Raw_QC_Normalized/01.Raw_to_QC/Important_files/GDPH.pheno.txt", row.name = FALSE)




###############################Step3: obtain the residuals: The DNA methylation data was adjusted for the above confounders using a linear model.
library(data.table)
library(limma)

# Read normalized methylation values
mval = as.matrix(fread(file = paste0(Important, "GDPH.DKD_ESRD.Normalized_Mval.txt")), rownames = 1)

# Read phenotype data
covariance = read.table(paste0(Important, "GDPH.DKD_ESRD.pheno.txt"), sep = " ", header = T)
Sample_Group = covariance$Sample_Group

# Read control probe PCA data
ctrl = read.table(paste0(Important, "GDPH.DKD_ESRD.CtrlProbe_PC.txt"), header = T, sep = " ")
PC = ctrl

# Read batch effect variables
Date = as.factor(as.character(covariance$Date))
Sentrix_ID = as.factor(as.character(covariance$Sentrix_ID))

# Read age, gender, and other covariates
sex_age = covariance[, c(2, 3)]

# Read cell type data
cell = read.table(paste0(Important, "GDPH.DKD_ESRD.Cell.Type.txt"), header = T, sep = " ", row.names = 1)
cell = cell[, 1:6]

# Read smoking score
Smoking = covariance$smokingScore

# Read scanner data
scan = covariance$scanner

# Combine all covariates
pheno = cbind(sex_age, PC, cell, Smoking, scan)

# Adjust methylation values for confounders using linear model
res_mval = removeBatchEffect(mval, covariates = pheno)
res_bvals = rnb.mval2beta(res_mval)

# Perform PCA on residual methylation values
res_mval_t = t(res_mval)
res_mvals_pca = prcomp(res_mval_t, center = TRUE, scale = TRUE)
res_mvals_sum_pc = summary(res_mvals_pca)
res_importance = res_mvals_sum_pc$importance[, 1:100]
res_mvals_pca = res_mvals_pca$x[, 1:100]

# Save the residual methylation PCA data
res_PC = res_mvals_pca[, 1:5]
pheno = cbind(pheno, res_PC)
new_mval = removeBatchEffect(mval, covariates = pheno)
res_mval = new_mval
res_bvals = rnb.mval2beta(res_mval)

# Save the residual methylation results
setwd(Important)
write.table(res_mvals_pca, file = paste0(Important, "GDPH.DKD_ESRD.residual_PCA.txt"), sep = "\t", col.names = TRUE, row.names = TRUE)
write.table(res_importance, file = paste0(Important, "GDPH.DKD_ESRD.residual_PCA_importance.txt"), sep = "\t", col.names = TRUE, row.names = TRUE)

write.table(res_mval, file = paste0(Important, "residual.mval.txt"), row.names = TRUE, quote = FALSE)
