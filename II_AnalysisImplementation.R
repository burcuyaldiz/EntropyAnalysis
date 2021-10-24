library(dplyr)
library("optparse")
source("3-way_function_2015.R")
source("2way_TIG_function_2015.R")
option_list = list(
  make_option(c("--caseFile"), type="character", default=NULL, 
              help="Frequencies of triplets or pairs of SNPs for case group ", metavar="character"),
    make_option(c("--controlFile"), type="character", default=NULL, 
              help="Frequencies of triplets or pairs of SNPs for control group ", metavar="character"),
		make_option(c("--intType"), type="character", default=NULL, 
              help="Type of interaction:3-WII or 2WI", metavar="character"),
			make_option(c("--outputFile"), type="character", default=NULL, 
              help="Interaction information results output path", metavar="character")	  ); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
casePath=opt$caseFile
controlPath=opt$controlFile
interactionType=opt$intType
outputPath=opt$outputFile


implementII <-function (caseFile, controlFile,interactionType ,outputFile){
	result_case <- read.csv(casePath)
	result_control <- read.csv(controlFile)

	if(interactionType=="3WII")
	{
		TIIG_Results<-data.frame(SNP_1<-as.character(),SNP_2<-as.character(),SNP_3<-as.character(),
							IIG<-as.numeric(),TIIG<-as.numeric(),TIIG_p_value<-as.numeric()) 
	
		for (i in 1:nrow(result_case))
		{
			
			X1<-c(result_control[i,1],result_control[i,2],result_control[i,3],  
				result_control[i,4],result_control[i,5],result_control[i,6],
				result_control[i,7],result_control[i,8],result_control[i,9],
				result_control[i,10],result_control[i,11],result_control[i,12],
				result_control[i,13],result_control[i,14],result_control[i,15],
				result_control[i,16],result_control[i,17],result_control[i,18],
				result_control[i,19],result_control[i,20],result_control[i,21],
				result_control[i,22],result_control[i,23],result_control[i,24],
				result_control[i,25],result_control[i,26],result_control[i,27])
				
			Y1<-c(result_case[i,1],result_case[i,2],result_case[i,3],  
				result_case[i,4],result_case[i,5],result_case[i,6],
				result_case[i,7],result_case[i,8],result_case[i,9],
				result_case[i,10],result_case[i,11],result_case[i,12],
				result_case[i,13],result_case[i,14],result_case[i,15],
				result_case[i,16],result_case[i,17],result_case[i,18],
				result_case[i,19],result_case[i,20],result_case[i,21],
				result_case[i,22],result_case[i,23],result_case[i,24],
				result_case[i,25],result_case[i,26],result_case[i,27])
		
			
		
			T_IIG_pvalue = p_value_func_T_TIIG(X1, Y1)
			temp_results<-data.frame(SNP_1<-result_case[i,28],SNP_2<-result_case[i,29],
									SNP_3<-result_case[i,30],IIG<-T_IIG_pvalue[1], TIIG<-T_IIG_pvalue[2],
									TIIG_p_value<-T_IIG_pvalue[3])
			TIIG_Results<-rbind(TIIG_Results,temp_results,stringsAsFactors = FALSE)		
			temp_results<-c()
		
		}
		colnames(TIIG_Results)=c("SNP1","SNP2","SNP3","IIG_Value","T_IIG_Value","T_IIG_p_value")
			
	}
	else
	{
		TIIG_Results<-data.frame(SNP_1<-as.character(),SNP_2<-as.character(),
                         TIG<-as.numeric(),TIG_p_value<-as.numeric())
		for (i in 1:nrow(result_case))
		{
			
			X1<-c(result_control[i,1],result_control[i,2],result_control[i,3],  
				result_control[i,4],result_control[i,5],result_control[i,6],
				result_control[i,7],result_control[i,8],result_control[i,9])
				
			Y1<-c(result_case[i,1],result_case[i,2],result_case[i,3],  
				result_case[i,4],result_case[i,5],result_case[i,6],
				result_case[i,7],result_case[i,8],result_case[i,9])

			T_IG_pvalue = p_value_func_T_IG(X1, Y1)
			temp_results<-data.frame(SNP_1<-result_case[i,10],SNP_2<-result_case[i,11],
            TIG<-T_IG_pvalue[1],TIG_p_value<-T_IG_pvalue[2])
			TIIG_Results<-rbind(TIIG_Results,temp_results,stringsAsFactors = FALSE)
			temp_results<-c()
		
		}
		colnames(TIIG_Results)=c("SNP1","SNP2","TIG_Value","TIG_p_value")
			
	}			
	write.csv(TIIG_Results,outputPath,row.names=FALSE)	
}
	
implementII(casePath,controlPath,interactionType,outputPath)

