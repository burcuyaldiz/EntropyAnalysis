# EntropyAnalysis
This repository is for the entropy based 3WII analysis of the variants prioritized by GWAS analysis. 
In the first step variants should be filtered from PLINK files. If BED files are provided in the first step then these files should be converted to pedmap format by using the following plink command:

plink --bfile inputFileName --recode12 --tab --out outputFileName

Then the prioritized variants can be filtered from pedmap files with this command:

plink --file pedmapFile --out outputFileName  --recode12 --extract prioritizedVariantListFile

Afterwards case and control groups are filtered: 

plink --file filteredSNPs --filter-cases --recode12 --out filteredSNPs_Case 
plink --file filteredSNPs --filter-controls --recode12 --out filteredSNPs_Control

Lastly a list file which involves the list of members belongs to genotypes is generated with:

plink --file filteredSNPs_Case --list --out caseSNPList
plink --file filteredSNPs_Control --list --out controlSNPList

Then these list files are used as an input for genotype_freq_calculation.py script which calculates the frequencies of pairwise and triplet SNPs
These frequencies are used for implementing 3WII and 2WIG analysis (Fan, R, Zhong, M, Wang, 2011)
2way_TIG_function_2015.R and 3-way_function_2015.R functions are adapted from https://georgetown.app.box.com/s/ptf0niqqquc5m3zxstehdoormpgh7hq5

