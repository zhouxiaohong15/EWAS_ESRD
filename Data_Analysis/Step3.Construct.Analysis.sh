cov="/share/home/liujianjun/ESRD_Project/Part0.GDPH_CKD_ESRD_EWAS/Important_files/"
path1="/share/home/liujianjun/ESRD_Project/Part2.Subtype_ESRD_EWAS/CKD/step1.mval.files/"



Path="/share/home/liujianjun/ESRD_Project/Part2.Subtype_ESRD_EWAS/"
CKD1="/share/home/liujianjun/ESRD_Project/Part2.Subtype_ESRD_EWAS/CKD/step2.logistic.model.Rscript.Shell.files/"
CKD2="/share/home/liujianjun/ESRD_Project/Part2.Subtype_ESRD_EWAS/CKD/step2.logistic.model.result.files/"

Script="/share/home/liujianjun/ESRD_Project/Part0.GDPH_CKD_ESRD_EWAS/Script/"



cd ${Path}

##########################分割为101个残差文件

#!usr/bin/bash
for i in $(seq 1 100)
do
sed -n "`echo $[7775*($i-1)+2]`,`echo  $[7775*($i-1)+1+7775]p`"  ${cov}GDPH.CKD_ESRD.residual_mvals.txt  >  ${path1}residual_mvals_block$i.txt
done;

tail -n 85 ${cov}GDPH.A1.residual_mvals.txt > ${path1}residual_mvals_block101.txt;




#################批量生成R
  for i in {1..101}; do
    # 确保每个文件路径正确，替换占位符并输出到目标文件
    sed -e "s|\${path1}|${path1}|g" \
        -e "s|\${cov}|${cov}|g" \
        -e "s|\${CKD1}|${CKD1}|g" \
        -e "s|\${CKD2}|${CKD2}|g" \
        -e "s|\$i|$i|g" \
        ${Script}Step2.Differences.Analysis.R > "${CKD1}block_residual_logistic_${i}.R"
  done

########################批量生成shell
 for i in $(seq 1 101); do
    echo "/share/apps/R/4.4.1/bin/Rscript ${CKD1}block_residual_logistic_$i.R" > ${CKD1}block_residual_logistic_$i.sh
  done

#####################批量运行
for i in $(seq 1 101); do
    bsub -e ${CKD1}error$i.log -o ${CKD1}output$i.log < ${CKD1}block_residual_logistic_${i}.sh 
  done

