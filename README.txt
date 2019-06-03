

#########################################################################################################################################
# A package for probabilistic multi-omics data matching procedure (proMODMatcher)                                                       #
# Author: Dr. Eunjee Lee, Dr. Seungyeul Yoo, and  Dr. Jun zhu                                                                           #
# Requirement: R version (3.5.1 or later)                                                                                               #
# Other requirement: Require R package mnormt                                                                                           #
# Computational cost and CPU time : Alignment of 408 samples requires 503 MB memory, 802 seconds CPU time in CPU processor 3.50 GHz     #
#########################################################################################################################################

####################################################
# Method 1. use user specified cis_table           #
####################################################
Usage: Rscript Run_paired_alignment.R --arg1=type1 --arg2=type1_file --arg3=type2 --arg4=type2_file --arg5=method --arg6=cis_table   

Input data: The script requires two different types of data and matching table of cis-relationship between type1 and type2 data as input. 
1. type1 : Type of the first data profile. e.g. mRNA, Methylation, miRNA, CNV, and Protein.
2. type1_file : Name of type 1 profile. Type1 profile should be tab delimited txt file of matrix with rows with gene or probe names and columns are sample ID. The sample ID of type1 profile should match to type2 smaple ID. 
3. type2 : Type of the second data profile. e.g. mRNA, Methylation, miRNA, CNV, and Protein.
4. type2_file : Name of type2 profile. Type2 profile should be tab delimited txt file of matrix with rows with gene or probe names and columns are sample ID.
5. method: Sample alignment method. It should be MODMatcher or proMODMatcher.   
6. cis_table : Tab delimited txt file of table for cis-associations between probe names of type1 profile and probe name of type2 profile. The first column corresponds to probe ID of type1,the second colum corresponds to probe ID of type2. 

Output data: SA_type1_type2_MODMatcher.RData or SA_type1_type2_proMODMatcher.RData. Output data includes FinalMap, inMap, nCisObs. FinalMap indicates the mapping result between type1 and type2 profiles. inMap indicates the annotated sample ID between type1 and type2. nCisObs indicates the number of significant cis-association at each iteration.  
Sample code :  Rscript Run_paired_alignment.R --arg1=mRNA --arg2=./data/TCGA_BRCA_Array_tumor.txt --arg3=miRNA --arg4=./data/TCGA_BRCA_miRNASeq_tumor.txt --arg5=proMODMatcher --arg6=./data/Matching_array_miRNA.txt 


#############################################################################################################
# Method 2. use pre-defined cis-table for miRNA, CNV, Protein (RPPA only), Methylation (HM450 or HM27 only) #
# Please check ID from your data file and cis-table first.                                                  # 
#############################################################################################################

Usage: Rscript Run_paired_alignment.R --arg1=type1 --arg2=type1_file --arg3=type2 --arg4=type2_file --arg5=method --arg6=cis_table
Input data: The script requires two different types of data and matching table of cis-relationship between type1 and type2 data as input.
1. type1 : Type of the first data profile. It should be "mRNA". 
2. type1_file : Name of type 1 profile. Type1 profile should be tab delimited txt file of matrix with rows with gene or probe names and columns are sample ID. The sample ID of type1 profile should match to type2 smaple ID.
3. type2 : Type of the second data profile. It should be one of these: mRNA, Methylation, miRNA, CNV, and RPPA.
4. type2_file : Name of type2 profile. Type2 profile should be tab delimited txt file of matrix with rows with gene or probe names and columns are sample ID.
5. method: Sample alignment method. It should be MODMatcher or proMODMatcher.
6. cis_table : It should be "None".

Output data: SA_type1_type2_MODMatcher.RData or SA_type1_type2_proMODMatcher.RData. Output data includes FinalMap, inMap, nCisObs. FinalMap indicates the mapping result between type1 and type2 profiles. inMap indicates the annotated sample ID between type1 and type2. nCisObs indicates the number of significant cis-association at each iteration.
Sample code :  Rscript Run_paired_alignment.R --arg1=mRNA --arg2=./data/TCGA_BRCA_Array_tumor.txt --arg3=miRNA --arg4=./data/TCGA_BRCA_miRNASeq_tumor.txt --arg5=proMODMatcher --arg6=./data/Matching_array_miRNA.txt

