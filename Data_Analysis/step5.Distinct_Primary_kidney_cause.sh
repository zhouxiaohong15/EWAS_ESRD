# Directories Creation

# Main working directory for differential CpG analysis
mkdir -p /share/home/liujianjun/ESRD_Project/Part2.Subtype_ESRD_EWAS

# Directory for MVAL files
mkdir -p /share/home/liujianjun/ESRD_Project/Part2.Subtype_ESRD_EWAS/step1.mval.files

# Disease-specific directories
mkdir -p /share/home/liujianjun/ESRD_Project/Part2.Subtype_ESRD_EWAS/${disease}

# Directory for shell scripts
mkdir -p /share/home/liujianjun/ESRD_Project/Part2.Subtype_ESRD_EWAS/Important.Shell.Code


# Set working directories for various steps
wd="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/"
wd1="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/01.Raw_to_QC/Important_files/"
wd2="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/02.Diffrential_CpG_analysis/step1.mval.files/"
Path="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/"

# Directories for logistic model R scripts and results
diff1="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/$j/step2.logistic.model.Rscript.Shell.files/"
diff2="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/$j/step2.logistic.model.result.files/"


# Create directories for each disease type
for i in IgAN SLEN PKD HRD; do
  mkdir -p /share/home/liujianjun/ESRD_Project/Part2.Subtype_ESRD_EWAS/$i
done


# Create subdirectories for R script shell files for each disease type
for i in IgAN SLEN PKD HRD; do
  mkdir -p /share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/$i/step2.logistic.model.Rscript.Shell.files
  mkdir -p /share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/$i/step2.logistic.model.result.files
done






# Set working directories for various steps
wd="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/"
wd1="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/01.Raw_to_QC/Important_files/"
wd2="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/02.Diffrential_CpG_analysis/step1.mval.files/"
Path="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/"


# Batch replace R script placeholders for each disease type

for j in SLEN IgAN PKD HRD; do
  # Define disease-specific paths for logistic model scripts and results
  diff1="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/$j/step2.logistic.model.Rscript.Shell.files/"
  diff2="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/$j/step2.logistic.model.result.files/"

  # Loop over the range 1 to 101 to generate logistic model scripts for each iteration
  for i in {1..101}; do
    # Ensure the file paths are correct, replace placeholders, and output the result to the target file
    sed -e "s|\${wd1}|${wd1}|g" \
        -e "s|\${wd2}|${wd2}|g" \
        -e "s|\${diff1}|${diff1}|g" \
        -e "s|\${diff2}|${diff2}|g" \
        -e "s|\$j|$j|g" \
        -e "s|\$i|$i|g" \
        ${Path}Raw.Code.sh > "${diff1}/block_residual_logistic_${i}.R"
  done
done



####################### Step 2: Generate Shell Scripts for R Logistic Models #############################

# Loop through disease conditions: SLEN, IgAN, PKD, HRD
for j in SLEN IgAN PKD HRD ; do
  diff1="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/$j/step2.logistic.model.Rscript.Shell.files/"

  # Loop through each block from 1 to 101 to generate individual R script shell files
  for i in $(seq 1 101); do
    echo "/share/apps/R/4.4.1/bin/Rscript ${diff1}block_residual_logistic_$i.R" > ${diff1}block_residual_logistic_$i.sh
  done
done

####################### Step 3: Submit Jobs to Run Shell Scripts on the Cluster #########################

# Loop through disease conditions: SLEN, IgAN, PKD, HRD
for j in SLEN IgAN PKD HRD; do
  diff1="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/$j/step2.logistic.model.Rscript.Shell.files/"

  # Submit jobs for blocks 21 to 40 using nohup to run in background
  for i in $(seq 1 101); do
    bsub   < ${diff1}block_residual_logistic_${i}.sh 
  done
done

####################### Step 4: Find Failed Scripts and Resubmit #########################

# Loop through disease conditions: SLEN, IgAN, PKD, HRD
for j in SLEN IgAN PKD HRD; do
  diff1="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/$j/step2.logistic.model.Rscript.Shell.files/"
  diff2="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/$j/step2.logistic.model.result.files/"

  # Create a list of all block numbers (1 to 101)
  seq 1 101 > ${diff2}all_numbers.txt

  # Get existing block numbers from results
  echo ${diff2}residual_result_block* | grep -oE '[0-9]+' | sort -n | uniq > ${diff2}existing_numbers.txt

  # Find the failed blocks by subtracting existing blocks from all numbers
  fail=$(comm -23 ${diff2}all_numbers.txt ${diff2}existing_numbers.txt | tr '\n' ' ')

  # Resubmit the failed blocks using nohup
  for i in ${fail}; do
    nohup sh ${diff1}block_residual_logistic_${i}.sh > ${diff1}output_${i}.log 2>&1 &
  done
done

####################### Step 5: Format the Results into a Unified Output File #########################

# Prepare directory for annotation results
for j in SLEN IgAN PKD HRD; do
  mkdir -p /share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/$j/step3.annotation.results.files/
done

# Format the logistic regression results
for j in SLEN IgAN PKD HRD ; do
  dir1="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/$j/step3.annotation.results.files/"
  diff2="/share/home/liujianjun/ESRD_Project/Back_up_data/tf11_data/PartA1.GDPH_ESRD_CKD/03.Subtype_Diffrential_CpG_analysis/$j/step2.logistic.model.result.files/"
  out="/share/home/liujianjun/ESRD_Project/Part2.Subtype_ESRD_EWAS/$j/step3.annotation.results.files/"

  # Initialize the linear result file
  echo "" > ${dir1}$j.linear_result_all.txt

  # Combine the residual results from all blocks into one file
  cat ${diff2}residual_result_block* >> ${dir1}$j.linear_result_all.txt

  # Remove unnecessary rows and clean up the file
  sed -i '/Estimate/d' ${dir1}$j.linear_result_all.txt  # Remove rows with the "Estimate" string
  sed -i "s/\"/ /g" ${dir1}$j.linear_result_all.txt  # Remove quotation marks
  sed -i "1iCpG   Estimate   SE   Pval   Trait" ${dir1}$j.linear_result_all.txt  # Add header
  sed -i "s/^ //g" ${dir1}$j.linear_result_all.txt  # Remove leading spaces
  sed -i "2d" ${dir1}$j.linear_result_all.txt  # Remove the second line (header from concatenation)
  cp ${dir1}$j.linear_result_all.txt ${out}$j.linear_result_all.txt
done

