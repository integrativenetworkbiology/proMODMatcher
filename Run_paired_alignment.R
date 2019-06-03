#!/usr/bin/Rscript

args<-commandArgs(TRUE)


first<-args[grep("arg1",args)]
type1 <- sub("--arg1=","",first)
second<-args[grep("arg2",args)]
type1_file <- sub("--arg2=","",second)
third<-args[grep("arg3",args)]
type2<-sub("--arg3=","",third)
fourth<-args[grep("arg4",args)]
type2_file<-sub("--arg4=","",fourth)
fifth<-args[grep("arg5",args)]
method<-sub("--arg5=","",fifth)
sixth<-args[grep("arg6",args)]
cis_table<-sub("--arg6=","",sixth)



type1_tumr_pro <- read.delim(type1_file, row.names=1, check.names=F)
type2_tumr_pro <- read.delim(type2_file, row.names=1, check.names=F)
if(cis_table=="None") { 
	if(!type1 =="mRNA") { print("wrong type1") } 
	else if (type2=="mRNA") { 
		cis_table <- cbind(rownames(type1_tumr_pro), rownames(type1_tumr_pro))
	} else if(type2=="miRNA") { 
		cis_table <- read.delim("./data/Matching_array_miRNA.txt", header=T)
	} else if (type2 =="HM27") { 
		cis_table <- read.delim("./data/Matching_array_MethylationHM27.txt", header=T)
	} else if (type2 =="HM450") {
                cis_table <- read.delim("./data/Matching_array_MethylationHM450.txt", header=T)
	} else if (type2 =="CNV") { 
                cis_table <- cbind(rownames(type1_tumr_pro), rownames(type1_tumr_pro))
	} else if (type2 =="RPPA") { 
		cis_table <- read.delim("./data/Matching_array_protein.txt", header=T)			
	} else { 
		print ("wrong type2") 
	}
} else { 
cis_table <- read.delim(cis_table, header=T)
}
cismatch <- apply(cis_table, 1, function(x) { tt <- c(match(x[1], rownames(type1_tumr_pro)), match(x[2], rownames(type2_tumr_pro)))
					      sum(is.na(tt))})
cis_table_filter <- cis_table[which(cismatch==0),]


type1_tumr_pro <- type1_tumr_pro[match(cis_table_filter[,1], rownames(type1_tumr_pro)),]
type2_tumr_pro <- type2_tumr_pro[match(cis_table_filter[,2], rownames(type2_tumr_pro)),]
rownames(type1_tumr_pro) <- rownames(type2_tumr_pro) <- paste(rownames(type1_tumr_pro), rownames(type2_tumr_pro), sep=":")
tumr_gene_common <- intersect(rownames(type1_tumr_pro), rownames(type2_tumr_pro))


type1_tumr_pro_com<-type1_tumr_pro[match(tumr_gene_common, rownames(type1_tumr_pro)),]
type2_tumr_pro_com<-type2_tumr_pro[match(tumr_gene_common, rownames(type2_tumr_pro)),]

dd1 <- apply(type1_tumr_pro_com, 1, function(x) max(table(x))>(dim(type1_tumr_pro_com)[2]*0.95))
dd2 <- apply(type2_tumr_pro_com, 1, function(x) max(table(x))>(dim(type2_tumr_pro_com)[2]*0.95))
tmp <- union(which(dd1), which(dd2))
if(length(tmp)>0) {
type1_tumr_pro_com <- type1_tumr_pro_com[-tmp,]
type2_tumr_pro_com <- type2_tumr_pro_com[-tmp,]
}


### Load profile 1
type1_pro <- type1_tumr_pro_com
na_ind<-unique(which(is.na(type1_pro),arr.ind=T)[,1])
if(length(na_ind)>0) { type1_pro<-type1_pro[-(na_ind),] }

### Load profile_2
type2_pro <- type2_tumr_pro_com
na_ind<-unique(which(is.na(type2_pro),arr.ind=T)[,1])
if(length(na_ind)>0) { type2_pro<-type2_pro[-(na_ind),] }

nn <- intersect(rownames(type1_pro) , rownames(type2_pro))
type1_pro <- type1_pro[nn,]
type2_pro <- type2_pro[nn,]

type1_gene<-rownames(type1_pro)
type1_sam<-colnames(type1_pro)
type2_gene<-rownames(type2_pro)
type2_sam<-colnames(type2_pro)

## generate initial mapping result 
sam_common <- intersect(colnames(type1_pro), colnames(type2_pro))
align <- as.data.frame(cbind(sam_common, sam_common, rep(1, length(sam_common))))

## filtering genes with low expression levels 
dd <- apply(as.matrix(1:dim(type1_pro)[1]), 1, function(x) { tmp <- rbind(type2_pro[x,sam_common], type1_pro[x,sam_common])
                                                             tt <- apply(tmp, 2, function(y) length(which(y==0)))
                                                             if(length(which(tt>0))>length(sam_common)*0.75) { FALSE} else { TRUE}})
type1_pro <- type1_pro[which(dd),]
type2_pro <- type2_pro[which(dd),]
if(type2 =="Methylation") { 
type2_pro <- -type2_pro
}
if( type1=="Methylation"){ 
type1_pro <- -type1_pro
}

inMap<-align
nCis_detected <- c()

if(method=="MODMatcher") { 
source("paired_alignment_core_simple.R")
FinalMap <- map
nCisObs <- nCis_detected
save(FinalMap, inMap, nCisObs,  file=paste("SA_", type1,"_", type2, "_simple.RData", sep="") )
} else if(method=="proMODMatcher") { 
sprior<- 1/dim(type1_pro)[2]
source("paired_alignment_core.R")
FinalMap <- map
nCisObs <- nCis_detected
save(FinalMap, inMap, nCisObs,  file=paste("SA_", type1,"_", type2, "_mvprob.RData", sep="") )
} else { 
print("Wrong method")
}


tt <- apply(FinalMap, 1, function(x) x[1]==x[2])
cat(length(sam_common), ' out of ', dim(type1_pro)[2], ' are common', '\n')
cat(dim(FinalMap)[1], ' out of ', dim(type1_pro)[2], ' were mapped','\n')
cat(length(which(tt)), ' out of ', length(sam_common), ' were self-mapped','\n')
cat(length(which(!tt)), ' were cross-mapped','\n')

