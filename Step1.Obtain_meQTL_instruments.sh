
#############################################Step1:extracting meQTL data of ESRD DMPs from a published Chinese dataset (Sijia Wang lab)


# Define working directory
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartD2.ESRD_MR_causal_analysis/Important.files/"

# Step 1: Filter rows where p-value < 5e-8 from the NSPT dataset and remove quotes
awk -F',' '$5 < 5e-8' /hwfssz1/pub/database/download.big.ac.cn/OMIX/OMIX004116/OMIX004116-01.csv > ${wd}NSPT.meQTL.csv
sed -i 's/\"//g' ${wd}NSPT.meQTL.csv

# Step 2: Extract CpG list from ESRD DMP file and save it to a new file
cd ${wd}
awk -F' ' '{print $1}' ${wd}Variable.Bacon.DMP.txt > ${wd}cpg_list.txt

# Step 3: Extract meQTL data for CpGs present in the list from the NSPT dataset
awk -F' ' 'NR==FNR {cpg[$1]; next} $5 in cpg' ${wd}cpg_list.txt ${wd}NSPT.meQTL.csv > ${wd}tmp.ESRD.DMP.meQTL.txt

# Step 4: Reformat the meQTL file to include separate columns for SNP, Chromosome, Position, and other relevant details
awk -F ' ' '{split($2, array, ":"); print $1 " " array[1] " " array[2] " " $3 " " $4 " " $5 " " $6 " " $7 " " $8}' ${wd}tmp.ESRD.DMP.meQTL.txt > ${wd}ESRD.DMP.meQTL.txt

# Step 5: Clean up potential carriage return characters and add column headers
sed -i 's/\r//' ${wd}ESRD.DMP.meQTL.txt
sed -i '1iSNP Chr BP EA OA DNAm Beta Tval Pval' ${wd}ESRD.DMP.meQTL.txt

# Step 6: Count the number of unique p-values in the file
awk '{print $5}' ${wd}ESRD.DMP.meQTL.txt | sort | uniq | wc -l

# Step 7: Clean up temporary files and rename the output file
rm ${wd}tmp.ESRD.DMP.meQTL.txt
rm tmp*
mv ${wd}ESRD.DMP.meQTL.txt ${wd}ESRD.1390.DMP.meQTL.txt




###################################################Step2. add RsID column to data

#############################################add rsID to SNP  in  GWAS data
# Set the working directories
wd1="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/GWAS.OF.ESRD/Important.file/"
wd0="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/DBSNP/"
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartD2.ESRD_MR_causal_analysis/Important.files/"

# Merge ESRD.GWAS.Results.txt and rsID.hg38 files based on the first two columns (rsid and chromosomal position)
awk '
    NR==FNR {                 # If reading the first file (ESRD.GWAS.Results.txt)
        key = $1 " " $2       # Create a key using the first two columns (rsid and chromosomal position)
        array[key] = $0       # Store the whole line in an array with the key as the identifier
        next                  # Skip to the next line
    }
    {                         # Process the second file (rsID.hg38)
        key = $1 " " $2       # Create the same key (rsid and chromosomal position)
        if (key in array) {   # If the key exists in the array from the first file
            print array[key], $3  # Print the matching line from the first file and the third column from the second file
        }
    }
' ${wd}ESRD.GWAS.Results.txt ${wd0}rsID.hg38 > ${wd}ESRD.GWAS.rsID.Results.txt   # Output the merged result to a new file


###################################################add rsID to SNP  in  meQTL data
# Set working directories
wd0="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/DBSNP/"
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/PartD2.ESRD_MR_causal_analysis/Important.files/"
wd1="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/GWAS.OF.ESRD/Important.file/"

# Merge ESRD.1390.DMP.meQTL.new.txt and dbsnp_151_hg37.vcf files based on the specified key (rsid and position)
awk '
    NR==FNR {                             # If reading the first file (ESRD.1390.DMP.meQTL.new.txt)
        key = $2 "\t" $3                  # Create a key using the second and third columns (chromosomal position and some identifier)
        array[key] = $0                   # Store the entire line in an array with the key
        next                               # Skip to the next line
    }
    {                                      # Process the second file (dbsnp_151_hg37.vcf)
        key = $1 "\t" $2                  # Create a key using the first two columns (rsid and chromosomal position)
        if (key in array) {               # If the key exists in the array (from the first file)
            print array[key], $3          # Print the matching line from the first file and the third column from the second file
        }
    }
' ${wd}ESRD.1390.DMP.meQTL.new.txt ${wd0}dbsnp_151_hg37.vcf > ${wd}ESRD.1390.DMP.meQTL.rsID.txt  # Output the merged results to a new file
