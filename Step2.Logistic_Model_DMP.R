########################## Generate 101 Batch R Files for Logistic Model ##########################

# Load residual M-values for the current iteration (i)
row$i = read.table("${wd2}residual_mvals_block$i.txt", header = FALSE, sep = " ")
row.names = row$i[, 1]   # Set row names as the first column of data
data$i = t(row$i[, -1])   # Transpose the data to get the M-values
colnames(data$i) = row.names  # Set column names as row names

# Read in phenotype data (covariates)
covariance = read.table("${wd1}GDPH.A1.pheno.txt", sep = " ", header = TRUE)

# Read cell type information (selecting the first 6 columns)
cell = read.table("${wd1}GDPH.A1.Cell.Type.txt", sep = " ", header = TRUE)
cell = cell[, 1:6]

# Read ASA chip genotype information (PCA components)
ASA = read.table("${wd1}GDPH.A1.genotype_PCA.txt", header = TRUE)
genotypePC = ASA[, 2:11]  # Extract PCA columns 2 to 11

# Extract age and gender data from covariance
Age = covariance$Age
Gender = covariance$Gender

# Extract sample group, smoking status, and other relevant variables
Sample_group = covariance$Sample_Group
Smoking = covariance$Smoking

# Identify the rows corresponding to the disease condition for logistic regression
ID = which(pheno$$j %in% c(0, 1))  # Disease groups for current sample

# Initialize result variable to store logistic regression output
result = c()

# Loop through each CpG site (column in the M-values data)
for (i in 1:nrow(file)) {
  # Perform logistic regression for the current CpG site (i)
  fit = glm(Sample_group[ID] ~ data$i[ID, i] + Age[ID] + Smoking[ID] + 
            cell[ID, 1] + cell[ID, 2] + cell[ID, 3] + cell[ID, 4] + 
            cell[ID, 5] + cell[ID, 6], family = binomial)
  
  # Extract the coefficient estimates (Estimate, SE, P-value) for the logistic regression
  result = rbind(result, c(colnames(data$i)[i], coef(summary(fit))[2, c(1, 2, 4)]))
}

# Assign column names to the result data frame
colnames(result) = c("CpG", "Estimate", "SE", "Pvalue")

# Convert P-values and coefficients to numeric and format them properly
result$Pvalue = as.numeric(result$Pvalue)
result$SE = as.numeric(result$SE)
result$Estimate = as.numeric(result$Estimate)
result$Trait = "$j"  # Add the disease trait information

# Define a custom formatting function for scientific notation and rounding
custom_format <- function(x) {
  if (abs(x) < 1e-3) {
    return(format(x, scientific = TRUE, digits = 3))  # Format in scientific notation if very small
  } else if (abs(x) >= 0.001 && abs(x) < 0.005) {
    return(format(round(x, 3), nsmall = 3))  # Round to 3 decimals if in the range
  } else {
    return(format(round(x, 2), nsmall = 2))  # Round to 2 decimals otherwise
  }
}

# Apply custom formatting to P-values, Estimates, and SE
result$Pvalue = sapply(result$Pvalue, custom_format)
result$Estimate = sapply(as.numeric(result$Estimate), custom_format)
result$SE = sapply(as.numeric(result$SE), custom_format)

# Sort the results by P-value (ascending order)
result = result[order(result$Pvalue), ]

# Display the top rows of the result
head(result)

# Save the result to a file
write.table(result, file = "${diff2}residual_result_block$i", sep = " ", row.names = FALSE, quote = FALSE)
