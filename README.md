# EntropyAnalysis
This repository is for the entropy based 3WII analysis of the variants prioritized by GWAS analysis. 
In the first step variants should be filtered from PLINK files. If BED files are provided in the first step then these files should be converted to pedmap format by using the following plink command:

plink --bfile inputFileName --recode12 --tab --out outputFileName

THen the prioritized variants can be filtered from pedmap files with this command:

plink --file pedmapFile --out outputFileName  --recode12 --extract prioritizedVariantListFile


