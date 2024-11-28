######批量生成R对应的shell文件
group="adipose	blood	bone	brain	digestive	endocrine	endothelial	epithelial	esc	es	eye	heart	hsc	ipsc	kidney	liver	lung	lymphoblastoid	mesench	muscle	myosat	neurosph	pancreas	placenta	pns	reproductive	smmuscle	spleen	stromal	thymus	urinary";

for i in ${group} ; 
do 
echo "
/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/software/anaconda3/bin/Rscript /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/permutation.Enh.${i}.Bino.R " > /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/permutation.Enh.${i}.Bino.sh;
done;


##################提交所有shell文件到作业

wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/";
wd1="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/";

rm ${wd1}batch.qsub.enrichment.sh;

for i in ${group}

do 
echo "qsub -cwd -l vf=5G,p=1 -q st.q -P P18Z10200N0124  -binding linear:8 ${wd}permutation.Enh.${i}.Bino.sh" >> ${wd1}batch.qsub.enrichment.sh;
done;

sh ${wd1}batch.qsub.enrichment.sh;
