import json
import configparser
import argparse
import pandas as pd

parser=argparse.ArgumentParser()
parser.add_argument('--snpListFile', type=str,help='Prioritized SNPs genotype list file path for case group')
parser.add_argument('--interactionType',type=str,help='Type of interaction information can be 3WII or 2WI')
parser.add_argument('--outputFilePath',type=str,help='Json file path for the frequencies of triplets or pairs of SNPs')
args=parser.parse_args()

###########################################################
#This function remove the family IDs and recode the genotypes in case control genotyping files
###########################################################

def convertGtypeCoding(snpListFilePath):
    f = open(snpListFilePath, "r")
    #g = open(controlFilePath, "r")
    snplist = f.readlines()
    f.close()
    for n, i in enumerate(snplist):
        snplist[n]=snplist[n].split()
        temp_list=snplist[n][3:len(snplist[n])]
        #Remove the familiy IDs with this step
        del temp_list[0::2]
        temp_list=list(dict.fromkeys(temp_list))
        snplist[n][3:len(snplist[n])]=temp_list

    for i in range (0,len(snplist)-1):
        if snplist[i][2]=='11':
            snplist[i][2]='0'
        elif snplist[i][2]=='12' or snplist[i][2]=='21' :
            snplist[i][2]='1'
        elif snplist[i][2] == '22' :
            snplist[i][2]='2'
    return(snplist)


##############################################################
#This function generates a csv file with frequencies of triplets or pairs
###############################################################
def convertgtypeII(snplist,interactionType,outputPath):
    freq="Freq_"
    nofcom={}
    temp_dic={}
    x=0
    if interactionType=="3WII":
        for i in range (0,len(snplist)-11,4) :
            for j in range (i+4,len(snplist)-7,4):
                for k in range (j+4,len(snplist)-3,4):
                    for l in range (i,i+4):       
                        for m in range (j,j+4):
                            for n in range (k,k+4):
                                if l%4 != 3 and m%4 != 3 and n%4 != 3:
                                    name_col=freq+snplist[l][2]+snplist[m][2]+snplist[n][2]
                                    set_1=set(snplist[l][3:len(snplist[l])])
                                    set_2=set(snplist[m][3:len(snplist[m])])
                                    set_3=set(snplist[n][3:len(snplist[n])])
                                    temp_dic.update({name_col:len(set_1 & set_2 & set_3)})
                    if l%4==3:
                        temp_dic.update({'SNP_1':snplist[l][1],'SNP_2':snplist[m][1],'SNP_3':snplist[n][1]})
                        nofcom.update({x:temp_dic})
                        x=x+1
                        temp_dic={}
    else:
        for j in range (0,len(snplist)-7,4):
            for k in range (j+4,len(snplist)-3,4):       
                for m in range (j,j+4):
                    for n in range (k,k+4):
                        if m%4 != 3 and n%4 != 3:
                            name_col=freq+snplist[m][2]+snplist[n][2]
                            set_2=set(snplist[m][3:len(snplist[m])])
                            set_3=set(snplist[n][3:len(snplist[n])])
                            temp_dic.update({name_col:len(set_2 & set_3)})
                if m%4==3:
                    temp_dic.update({'SNP_1':snplist[m][1],'SNP_2':snplist[n][1]})
                    nofcom.update({x:temp_dic})
                    x=x+1
                    temp_dic={}
        
    nofcom=pd.DataFrame.from_dict(nofcom,orient='index')
    nofcom.to_csv(outputPath,index=False)
    return 


interaction=args.interactionType
snpListFilePath=args.snpListFile
outputPath=args.outputFilePath

snpListEncoded=convertGtypeCoding(snpListFilePath)
convertgtypeII(snpListEncoded,interaction,outputPath)


