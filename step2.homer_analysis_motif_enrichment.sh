wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartG.Functional_analysis/01.homer_motif_analysis/Important_files/"

sed -i "s/\"//g" ${wd}DMPpeaks.txt;
/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/software/anaconda3/bin/findMotifsGenome.pl  ${wd}DMPpeaks.txt hg19 ${wd}Motif.Results -size 200 -mask -S 20 

#Number of motifs to find ("-S <#>", default 25)
#-mask The "-mask" is optional and tells the program to use the repeat-masked sequence.
#-size for Transcription Factor peaks, most of the motifs are found +/- 50-75 bp from the peak center, making it better to use a fixed size rather than depend on your peak size.
#网址：http://homer.ucsd.edu/homer/ngs/peakMotifs.html
