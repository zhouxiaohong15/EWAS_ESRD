########################################Bacon Analysis: Adjusted the bias of EWAS based on summary statistics

# Define file path
index = "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Subtype_Meta_analysis/00.EWAS_results/"

# Load data
library(data.table)
data = fread(paste0(index, "All_merged_file.txt"))

# Add an ID column to data
data$ID = 1:nrow(data)

# Extract effect sizes (es) and standard errors (se) for SG and SY cohorts
es = cbind(data[, 2], data[, 6], data[, 10])  # Effect sizes from columns 2, 6, 10
se = cbind(data[, 3], data[, 7], data[, 11])  # Standard errors from columns 3, 7, 11

# Perform BACON correction
library(bacon)
bc = bacon(NULL, es, se)

# Extract BACON-adjusted statistics (optional)
Bacon.se = se(bc)
Bacon.estimate = es(bc)
Bacon.pval = pval(bc)

# Prepare separate datasets for each subgroup
data_SLEN = data.frame(CpG = data$CpG, Bacon_b_SLEN = Bacon.estimate[, 1], Bacon_se_SLEN = Bacon.se[, 1], Bacon_P_SLEN = Bacon.pval[, 1])
data_IgAN = data.frame(CpG = data$CpG, Bacon_b_IgAN = Bacon.estimate[, 1], Bacon_se_IgAN = Bacon.se[, 1], Bacon_P_IgAN = Bacon.pval[, 1])
data_PKD = data.frame(CpG = data$CpG, Bacon_b_PKD = Bacon.estimate[, 2], Bacon_se_PKD = Bacon.se[, 2], Bacon_P_PKD = Bacon.pval[, 2])
data_HRD = data.frame(CpG = data$CpG, Bacon_b_HRD = Bacon.estimate[, 3], Bacon_se_HRD = Bacon.se[, 3], Bacon_P_HRD = Bacon.pval[, 3])

# Load dplyr for merging datasets
library(dplyr)

# Merge datasets for IgAN, PKD, and HRD into one dataframe
merged_data = bind_cols(data_IgAN, select(data_PKD, -CpG), select(data_HRD, -CpG))

# Check number of rows in the merged data
row = nrow(merged_data)
row

# Save inflation factors (bias correction output is part of the BACON result)
# Perform fixed-effect meta-analysis and inspect the results (e.g., bias and inflation factors)


############################################ BACON Meta-analysis and Result Export

# Perform meta-analysis using BACON results
bcm <- meta(bc)

# Save results to files
save(bcm, file = paste0(index, "Bacon.Subtype.Meta.Results.Rdata"))
save(bc, file = paste0(index, "Bacon_Subtype.Rdata"))

# Add meta-analysis results to merged data
merged_data$Bacon_b_meta = es(bcm)[, 4]
merged_data$Bacon_se_meta = se(bcm)[, 4]
merged_data$Bacon_P_meta = pval(bcm)[, 4]

# Apply scientific notation and rounding to certain columns
library(dplyr)
df <- merged_data
df[, 2:13] <- lapply(df[, 2:13], function(col) {
  ifelse(abs(col) < 1e-3, as.numeric(format(col, scientific = TRUE, digits = 3)),
         ifelse(abs(col) < 0.005, round(col, 3), round(col, 2)))
})

# Sort the results by meta-analysis p-value
all <- df[order(df$Bacon_P_meta), ]

# Determine directionality based on the sign of beta values for each subtype
all$dir <- ifelse(
  (all$Bacon_b_HRD < 0 & all$Bacon_b_PKD < 0 & all$Bacon_b_IgAN < 0) |
  (all$Bacon_b_HRD > 0 & all$Bacon_b_PKD > 0 & all$Bacon_b_IgAN > 0),
  1, 0
)

# Save sorted results to file
write.table(all, paste0(index, "All_Bacon_Result.txt"), row.names = FALSE, quote = FALSE)

# Calculate cutoff values for each subgroup (IgAN, PKD, HRD) based on p-value thresholds
cut.IgAN = 0.05 / 4 / nrow(subset(all, Bacon_P_IgAN < 1e-5 / 4))
cut.PKD = 0.05 / 4 / nrow(subset(all, Bacon_P_PKD < 1e-5 / 4))
cut.HRD = 0.05 / 4 / nrow(subset(all, Bacon_P_HRD < 1e-5 / 4))

# Identify CpGs that are significant for each subgroup
IgAN.CpG = subset(all, Bacon_P_IgAN < 1e-5 / 4 & Bacon_P_PKD < cut.PKD & Bacon_P_HRD < cut.HRD)$CpG
PKD.CpG = subset(all, Bacon_P_PKD < 1e-5 / 4 & Bacon_P_IgAN < cut.IgAN & Bacon_P_HRD < cut.HRD)$CpG
HRD.CpG = subset(all, Bacon_P_HRD < 1e-5 / 4 & Bacon_P_PKD < cut.PKD & Bacon_P_IgAN < cut.IgAN)$CpG

# Combine unique CpGs across all subgroups
can.CpG = unique(c(IgAN.CpG, PKD.CpG, HRD.CpG))

# Identify DMPs with meta p-value < threshold and directionality = 1
DMP = subset(all, CpG %in% can.CpG & Bacon_P_meta < 1.6e-08 & dir == 1)

# Save the final list of DMPs (shared across subtypes)
write.table(DMP, paste0(index, "ESRD.Shared.DMPs.txt"), row.names = FALSE, quote = FALSE)

# Count the number of significant CpGs for each subgroup and meta-analysis
dat = data.frame(IgAN = length(IgAN.CpG), PKD = length(PKD.CpG), HRD = length(HRD.CpG), meta = nrow(DMP))
rownames(dat) = "count"

# Prepare an empty dataframe for meta-analysis results (heterogeneity and effect size)
meta_results <- data.frame(
  CpG = character(0),
  I2 = numeric(0),
  Q = numeric(0),
  Phet = numeric(0),
  stringsAsFactors = FALSE
)

# Loop to perform meta-analysis and calculate heterogeneity (I²) for each CpG
for(i in 1:nrow(DMP)) {
  # Extract beta and standard error for the current CpG from each subtype
  yi <- c(DMP$Bacon_b_IgAN[i], DMP$Bacon_b_PKD[i], DMP$Bacon_b_HRD[i])
  sei <- c(DMP$Bacon_se_IgAN[i], DMP$Bacon_se_PKD[i], DMP$Bacon_se_HRD[i])
  
  # Perform fixed-effect meta-analysis
  meta_result <- rma(yi, sei = sei, method = "FE")
  
  # Store the meta-analysis results (I², Q statistics, and heterogeneity p-value)
  meta_results <- rbind(meta_results, data.frame(
    CpG = DMP$CpG[i],
    I2 = meta_result$I2,
    Q = meta_result$QE,
    Phet = meta_result$QEp
  ))
}

# Save meta-analysis results
write.table(meta_results, paste0(index, "Meta_Analysis_Results.txt"), row.names = FALSE, quote = FALSE)


#############################################Remove Non-Variant CpG Sites (90%-10% quantile range < 0.005)

# Define file paths
index = "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartC.Bacon_Subtype_Meta_analysis/00.EWAS_results/"
SY_A1 = "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartA1.GDPH_ESRD_CKD/01.Raw_to_QC/Important_files/"
SG = "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartB.AHPL_SG_EWAS/1.Raw_QC_Normalized/01.Raw_to_QC/Important_files/"

# Read significant CpG sites from meta-analysis results
meta = read.table(paste0(index, "ESRD.Shared.DMP.txt"), header = TRUE, sep = " ")
CpG = meta$CpG  # Extract CpG sites

# Read annotation file
anno = read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartH.Target_gene_annotation/epimap_link/4_col_annotation.txt", header = TRUE, sep = " ")

# Annotate significant CpG sites (chromosome and coordinates)
sigSite_anno = anno[anno$CpG %in% CpG, ]

# Read preprocessed b-values for GDPH cohort (SY)
library(data.table)
SY = as.matrix(fread(paste0(SY_A1, "GDPH.A1.Normalized_Bval.txt")), rownames = 1)
SYsig.bval = SY[rownames(SY) %in% CpG, ]
SYsig.bval.t = t(SYsig.bval)  # Transpose matrix
write.table(SYsig.bval, paste0(index, "GDPH.Meta.DMP.Bvals"))

# Read preprocessed b-values for AHPL cohort (SG)
SG = as.matrix(fread(paste0(SG, "AHPL.Normalized_Bval.txt")), rownames = 1)
SG = SG[, -224]  # Remove unwanted column (e.g., if extra column is present)
SGsig.bval = SG[rownames(SG) %in% CpG, ]
SGsig.bval.t = t(SGsig.bval)  # Transpose matrix
write.table(SGsig.bval, paste0(index, "AHPL.Meta.DMP.Bvals"))

# Alternative approach (can be used instead of the previous block)
# SYsig.bval = as.matrix(fread("path_to_file/GDPH.Meta.DMP.Bvals"), rownames = 1)
# SGsig.bval = as.matrix(fread("path_to_file/AHPL.Meta.DMP.Bvals"), rownames = 1)

# Merge both matrices (GDPH and AHPL) into one
all.sig.bval = cbind(SYsig.bval, SGsig.bval)
all.sig.bval = t(all.sig.bval)  # Transpose for analysis

# Analyze the range (90%-10% quantile difference) for each CpG
ncol = ncol(all.sig.bval)
result = data.frame()

for (i in 1:ncol) {
  ave = mean(all.sig.bval[, i])
  sdv = sd(all.sig.bval[, i])
  range = quantile(all.sig.bval[, i], c(0.1, 0.9))
  delta = range[2] - range[1]
  result = rbind(result, c(colnames(all.sig.bval)[i], ave, sdv, range, delta))
}

colnames(result) = c("CpG", "Mean", "SD", "10%", "90%", "Range")
write.table(result, paste0(index, "Whether.VarCpG.Result.txt"))

# Filter CpGs with range <= 0.05 (non-variable CpGs)
failed.site = subset(result, Range <= 0.05)
failed.CpG = failed.site$CpG
length(failed.CpG)

# Remove non-variable CpGs from the annotated list
retained.CpG = sigSite_anno[!sigSite_anno$CpG %in% failed.CpG, ]
dim(retained.CpG)

# Save the retained variable CpG sites
write.table(retained.CpG, paste0(index, "Variable.Bacon.DMP.txt"))
save(all.sig.bval, file = paste0(index, "All.SYSG.Sig.Bval.txt"))

# Annotate retained DMPs with chromosome and coordinate information
all = read.table(paste0(index, "Meta.DMP.FEandRE.Result.txt"), header = TRUE)  # Update file path if needed
names(all)[1] = "CpG"  # Ensure column name consistency

DMP = read.table(paste0(index, "Variable.Bacon.DMP.txt"), header = TRUE)
overlap = intersect(all$CpG, DMP$CpG)

# Merge annotation info for overlapping CpGs
anno = subset(DMP, CpG %in% overlap)
anno = anno[, 1:4]  # Keep only relevant columns
result = merge(anno, all, by = "CpG")
result = result[order(result$pval.adj.meta), ]  # Sort by adjusted p-value

# Save the annotated DMP results
write.table(result, paste0(index, "Bacon.Annotate.DMP.txt"), row.names = FALSE)
