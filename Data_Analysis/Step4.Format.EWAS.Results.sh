 
dir1="/share/home/liujianjun/ESRD_Project/Part2.Subtype_ESRD_EWAS/CKD/step3.annotation.results.files/"

diff2="/share/home/liujianjun/ESRD_Project/Part2.Subtype_ESRD_EWAS/CKD/step2.logistic.model.result.files/"

 # Initialize the linear result file
  echo "" > ${dir1}CKD.linear_result_all.txt

  # Combine the residual results from all blocks into one file
  cat ${diff2}residual_result_block* >> ${dir1}CKD.linear_result_all.txt

  # Remove unnecessary rows and clean up the file
  sed -i '/Estimate/d' ${dir1}CKD.linear_result_all.txt  # Remove rows with the "Estimate" string
  sed -i "s/\"/ /g" ${dir1}CKD.linear_result_all.txt  # Remove quotation marks
  sed -i "1iCpG   Estimate   SE   Pval   Trait" ${dir1}CKD.linear_result_all.txt  # Add header
  sed -i "s/^ //g" ${dir1}CKD.linear_result_all.txt  # Remove leading spaces
  sed -i "2d" ${dir1}CKD.linear_result_all.txt  # Remove the second line (header from concatenation)

