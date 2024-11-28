############################################## Stroke MR Analysis

# Load required libraries
library(data.table)
library(TwoSampleMR)
library(ieugwasr)

# Outcome data format (GWAS results):
# chromosome | base_pair_location | effect_allele | other_allele | beta | standard_error | effect_allele_frequency | p_value | rsid | is_strand_flip

# Example data for outcome (Stroke)
# 4 3107715 G A 0.00840730 0.0065276 0.38120 0.19780 rs10015979 no
# 4 74272581 A G -0.00048932 0.0074825 0.22990 0.94790 rs10019396 no
# 4 74332086 T C -0.00248530 0.0127220 0.07013 0.84510 rs10023906 no

#  Example for outcome metadata:
# rs10015979 1286331 13 no 22 0.8494 --+--+-++-+++-+++--+-+ 68408 lo 10
# rs10019396 1216078 12 no 20 0.6196 -+-+??+-+++-----+-+--- 66089 lo 10
# rs10023906 1216078 12 no 20 0.2274 -+--??+----+++--++-+-+ 66089 lo 10

# Load exposure data (meQTL)
final = fread("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartD2.ESRD_MR_causal_analysis/Important.files/meQTL.exposure.data.txt")

# Read outcome data for Stroke
out1 <- read_outcome_data(
  filename = "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartD2.ESRD_MR_causal_analysis/ESRD.Complications/Stroke.outcome.data.txt",
  sep = "\t", snps = NULL, snp_col = "rsid", beta_col = "beta", se_col = "standard_error",
  effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value"
)
out1$outcome = "Stroke"  # Assign outcome label "Stroke"

# Harmonize exposure and outcome data
dat <- harmonise_data(exposure_dat = final, outcome_dat = out1)

# Save harmonized data
write.table(dat, "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartD2.ESRD_MR_causal_analysis/Important.files/Stroke.Harmonised.Data.txt", row.names = FALSE, quote = FALSE)

############################### MR Analysis: Wald Ratio

# Load harmonized data for Stroke
hm1 = fread("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartD2.ESRD_MR_causal_analysis/Important.files/Stroke.Harmonised.Data.txt")

# Filter data for p-value <= 0.05 for significant exposures
hm1 = subset(hm1, pval.exposure <= 0.05)

# Load biologically relevant CpGs for KEGG pathway analysis
biologic.CpGs = read.table("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartD2.ESRD_MR_causal_analysis/Important.files/Twopathway.CpGs.txt", header = FALSE)

# Filter harmonized data for exposures involved in KEGG pathways
hm1.kegg = subset(hm1, exposure %in% biologic.CpGs$V1)

# List of unique exposures for MR analysis
list = unique(hm1.kegg$exposure)

# Initialize container for results
file = c()

# Loop through each unique exposure and select the minimal p-value for each
for (i in 1:length(list)) {
  data = subset(hm1, exposure == list[i])
  new = data[which(data$pval.exposure == min(data$pval.exposure)), ]
  file = as.data.frame(rbind(file, new))
}

# Perform MR analysis using Wald ratio method
wald = mr(file, method_list = "mr_wald_ratio")
head(wald)

# Merge Wald results with file data
result = merge(file, wald, by = "exposure")

# Add empty columns for result formatting
result$id.exposure.y = c()
result$id.outcome.y = c()
result$eaf.exposure = c()
result$eaf.outcome = c()
result$id.outcome.x = c()
result$id.outcome.y = c()
result$data_source.outcome = c()
result$outcome.y = c()
result$samplesize.outcome = c()
result$palindromic = c()
result$remove = c()
result$ambiguous = c()
result$pval_origin.exposure = c()
head(result)

# Further result formatting
result$id.exposure.x = c()
result$data_source.exposure = c()
head(result)
result$mr_keep.outcome = c()
result$mr_keep.exposure = c()
result$outcome.x = c()

# Save the final MR results
write.table(result, "/
