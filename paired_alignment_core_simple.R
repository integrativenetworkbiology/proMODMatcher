#!/usr/bin/Rscript


#### to find cis gene pairs
source("calculateCisGeneCorrelation_simple.R")
map<-inMap;
outMap<-matrix(data='NONE', length(type1_sam), 3);
prev <- paste(map[,1], map[,2], sep="-")
new <- paste(outMap[,1], outMap[,2], sep="-")
iter=0;
while (length(setdiff(prev, new))+length(setdiff(new, prev))>0 & iter<10){  #  refine till no change in mapping
if (iter>0) {prev<-new}

cisGenes<-calculateCisGeneCorrelation_simple(type1_pro, type2_pro, map) 
gene_common <- cisGenes[[1]]
names(gene_common)<- cisGenes[[1]]
cis_rho<-cisGenes[[2]]
cis_pva<-cisGenes[[3]]
cis_qva<-cisGenes[[4]]
id_cis<-which(cis_rho>0 & cis_pva<0.01 & cis_qva<0.05)
cat('ngene=', length(gene_common), ' cis=', length(id_cis), '(used). \n')
## select only one pair based on rho 
gn <- apply(as.matrix(rownames(type2_pro)), 1, function(x) strsplit(x, split=":")[[1]][1])
pn <- apply(as.matrix(rownames(type2_pro)), 1, function(x) strsplit(x, split=":")[[1]][2])
gg<- unique(gn[duplicated(gn)])
gsetsd <- unique(unlist(apply(as.matrix(gg), 1, function(x) setdiff(which(gn==x), which(gn==x)[which.max(cis_rho[which(gn==x)])]))))
id_cis <- setdiff(id_cis, gsetsd)
test_id_cis <- id_cis
cat('ngene=', length(gene_common), ' cis=', length(id_cis), '(unique). \n')
nCis_detected <- c(nCis_detected,length(id_cis))


### Obtain common samples
sam_common <- intersect(type1_sam,type2_sam)
###### data matrices
m1 <- as.matrix(type1_pro[match(gene_common, rownames(type1_pro)), match(sam_common, colnames(type1_pro))])
m2 <- as.matrix(type2_pro[match(gene_common, rownames(type2_pro)), match(sam_common, colnames(type2_pro))])
ngene <- dim(m1)[1]
nsample <- dim(m1)[2]

#sample correlation
ngene_cis<- length(id_cis)
m1<-m1[id_cis,]
m2<-m2[id_cis,]
### standardization
m1<-t(apply (m1, 1, rank))
m2<-t(apply(m2, 1, rank))

### sample correlation
rs <- cor(m1, m2, method="spearman")

#check rank result for self-alignment
rkv1<-sapply (seq.int(nsample), function (i) nsample-rank(rs[i,])[i]+1)
rkv2<-sapply (seq.int(nsample), function (i) nsample-rank(rs[,i])[i]+1)
#check rank result for all possible pairs
rk1<-sapply (seq.int(nsample), function (i) nsample-rank(rs[i,])+1)
rk2<-sapply (seq.int(nsample), function (i) nsample-rank(rs[,i])+1)




#search each type1 profile against all type 2 profile
idx<-matrix(0,1, dim(type2_pro)[2])
idx[which(match(colnames(type2_pro), sam_common)>0)]<- 1
sam_type2_all<-c(sam_common, colnames(type2_pro)[which(idx==0)])
#search each type2 profile against all type 1 profile
idx<-matrix(0,1, dim(type1_pro)[2])
idx[which(match(colnames(type1_pro), sam_common)>0)]<- 1
sam_type1_all<-c(sam_common, colnames(type1_pro)[which(idx==0)])
sam_all <- c(sam_common, setdiff(sam_type1_all, sam_common), setdiff(sam_type2_all, sam_common))


m2_all<-as.matrix(type2_pro[match(names(gene_common[id_cis]), rownames(type2_pro)), match(sam_type2_all, colnames(type2_pro))])
m2_all<-t(apply(m2_all, 1, rank))
m1_all<-as.matrix(type1_pro[match(names(gene_common[id_cis]), rownames(type1_pro)), match(sam_type1_all, colnames(type1_pro))])
m1_all<-t(apply(m1_all, 1, rank))

rs_all <- cor(m1_all, m2_all, method="spearman")

#output alignment
outMap<-matrix(data='NONE', dim(rs_all)[1], 3);



#check for all possible alignment
RK1<-t(sapply (seq.int(dim(m1_all)[2]), function (i) dim(rs_all)[2]-rank(rs_all[i,])+1))
RK2<-(sapply (seq.int(dim(m2_all)[2]), function (i) dim(rs_all)[1]-rank(rs_all[,i])+1))

## self alignment threshold ; minimum of 5% or 20 , cross alignment threshold ; only  1st 
cutoff=0.05;
ntop  = nsample*cutoff
if(ntop>20) { ntop <- 20} 
for (i in 1:dim(inMap)[1]) {
  ix<-match(inMap[i,1], sam_common)
  iy<-match(inMap[i,2], sam_common)
  if(rk1[ix, iy]<=ntop & rk2[ix, iy]<=ntop) { 
    ## comparison with 
    outMap[ix, 1] <-sam_common[ix] 
    outMap[ix, 2] <-sam_common[iy]
    outMap[ix, 3] <-1
  }
}


### check all unmapped samples
ntop_cross = 1
for (i in 1:dim(RK1)[1]) {
  if (outMap[i,1] == 'NONE') {
    ix<-match(sam_type1_all[i], sam_type1_all)
    iy<-sort(RK1[ix,], index.return = TRUE)$ix[1]
    indic <- (iy > dim(inMap)[1])
    if(!indic) { indic2 <- is.na(match(inMap[iy, 2], outMap[,2]))} else { indic2 <- TRUE}
    if(indic2) {
    if(ix== sort(RK2[,iy], index.return = TRUE)$ix[1] & RK1[ix, iy]<=ntop_cross) { ### top at both directions
      ## comparison with
      outMap[ix, 1] <-sam_type1_all[ix]
      outMap[ix, 2] <-sam_type2_all[iy]
      outMap[ix, 3] <-1
    }
  }
  }
}

#### prepare for the next iteration
#cat( 'iteration #', iter,  'nT=', nT, ', nF=', nF,'\n')
map<-outMap[which(outMap[,3]==1),]
new <- paste(map[,1], map[,2], sep="-")
iter= iter+1;
}

###output map
map<-outMap[which(outMap[,3]==1),1:2]
colnames(map) <- c(type1, type2)
      
