#####################教程网址：https://hemtools.readthedocs.io/en/latest/content/Bioinformatics_tools/homer.html 以及 http://homer.ucsd.edu/homer/ngs/peakMotifs.html

##############该脚本的目的是对显著的重要的motif注释位置，TSS, intron, exon, etc.
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/01.homer_motif_analysis/Important_files/"
wd2="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/01.homer_motif_analysis/Important_files/Motif.Results/knownResults/"

/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/software/anaconda3/bin/annotatePeaks.pl ${wd}DMPpeaks.txt  -m  ${wd2}known1.motif ${wd2}known2.motif   ${wd2}known3.motif  -rmrevopp  hg19 -annStats   >  ${wd}DMP.Peak.KnownMotif.Anno.txt
