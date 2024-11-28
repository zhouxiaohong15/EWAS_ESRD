rm /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/All.Enrichment.Binomial.Pval.txt;

cat  /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/Final.Enrichment.Pval.*.txt  >> /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/All.Enrichment.Binomial.Pval.txt;

sed -i '/Pval/d' /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/All.Enrichment.Binomial.Pval.txt;

sed -i '/Sample.ID/d' /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/All.Enrichment.Binomial.Pval.txt;

sed -i '1i Sample.ID Obs.Frep Exp.Frep FC Pval' /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/All.Enrichment.Binomial.Pval.txt;

sed -i "s/\"//g" /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/All.Enrichment.Binomial.Pval.txt;

#awk '{print $2,$1}' /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/All.Enrichment.Binomial.Pval.txt > /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/All.Enrichment.Binomial.Pval.back.txt;

#rm /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/All.Enrichment.Binomial.Pval.txt;

#mv /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/All.Enrichment.Binomial.Pval.back.txt /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/04.chromatin_accessibility/Batch.permutation/All.Enrichment.Binomial.Pval.txt;
