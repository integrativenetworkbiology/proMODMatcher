#!/usr/bin/Rscript


#### to find cis gene pairs
source("MahaloanobisDistance_core.R")
source('calculateCisGeneCorrelation_simple.R')

#########
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
sam_common<-intersect(type1_sam,type2_sam)

#####search each type1 profile against all type 2 profile
idx<-matrix(0,1, dim(type2_pro)[2])
idx[which(match(colnames(type2_pro), sam_common)>0)]<- 1
sam_type2_all<-c(sam_common, colnames(type2_pro)[which(idx==0)])
#search each type2 profile against all type 1 profile
idx<-matrix(0,1, dim(type1_pro)[2])
idx[which(match(colnames(type1_pro), sam_common)>0)]<- 1
sam_type1_all<-c(sam_common, colnames(type1_pro)[which(idx==0)])
sam_all <- c(sam_common, setdiff(sam_type1_all, sam_common), setdiff(sam_type2_all, sam_common))

###### data matrices
m1 <- as.matrix(type1_pro[match(gene_common, rownames(type1_pro)), match(sam_type1_all, colnames(type1_pro))])
m2 <- as.matrix(type2_pro[match(gene_common, rownames(type2_pro)), match(sam_type2_all, colnames(type2_pro))])
ngene <- dim(m1)[1]
nsample <- dim(m1)[2]

#sample correlation
ngene_cis<- length(id_cis)
m1<-m1[id_cis,]
m2<-m2[id_cis,]
### standardization for real data only - inverse normal 
m1<-t(apply (m1, 1, function(x) qnorm((rank(x, na.last="keep")-0.5)/sum(!is.na(x)))))
m2<-t(apply (m2, 1, function(x) qnorm((rank(x, na.last="keep")-0.5)/sum(!is.na(x)))))
rs <- cor(m1, m2, method="spearman")
mean_cc <-  mean(diag(rs[inMap[,1], inMap[,2]]), na.rm=T)
sd_cc <- sd(diag(rs[inMap[,1], inMap[,2]]), na.rm=T)
corr_cutoff <- mean_cc-sd_cc
corr_cutoff_self <- mean_cc-2.576*sd_cc

v1 <- sapply(1:dim(m1)[2],  function(x) sd(m1[,x]))
v2 <- sapply(1:dim(m2)[2],  function(x) sd(m2[,x]))


#more random samples by permutation
moreS <-1000
if(moreS >0 ) {
  ##generate permuted samples
  m1s<-matrix(data=0, nrow=dim(m1)[1], ncol=moreS )
  m2s<-matrix(data=0, nrow=dim(m2)[1], ncol=moreS ) 
  for (i in 1:moreS) {
    pi = sample(1:dim(m1)[2])[1];
    permi = sample(1:dim(m1)[1]);
    m1s[,i]<-m1[permi,pi];
    pi = sample(1:dim(m2)[2])[1];
    permi = sample(1:dim(m2)[1]);
    m2s[,i]<-m2[permi,pi];
  }


rs_s <- cor(cbind(m1, m1s), cbind(m2, m2s), method="spearman")
tmp_prob <- t(sapply(seq.int(nsample), function(i) mnorm_Others(i, rs_s, dim(rs_s)[2])$prob))
rk1_prob <-tmp_prob[, 1:dim(m2)[2]]
rk1_p <- 1-rk1_prob
rkv1_prob <- diag(rk1_prob)
} else { 
rkv1_prob <- sapply(seq.int(nsample), function(i) mnorm_self_sim(i,rs)$prob)
rk1_prob <- t(sapply(seq.int(nsample), function(i) mnorm_Others(i, rs, nsample)$prob))
rk1_p <- 1-rk1_prob
}

rk1_p [which (rk1_p>0.05)] <-nsample
p_cutoff <-min(0.001, 1/nsample^1.25)

### adjust prior by input MAP
for (i in 1:dim(inMap)[1]) {
  ix<-match(inMap[i,1], sam_common)
  iy<-match(inMap[i,2], sam_common)
  rk1_p[ix, iy] <- rk1_p[ix,iy]*sprior
}

#output alignment
outMap<-matrix(data='NONE', dim(rk1_p)[1], 3);

### check input MAP
if(length(id_cis)>1000) {ntop = 1 } else { ntop=2} 
for (i in 1:dim(inMap)[1]) {
  ix<-match(inMap[i,1], sam_common)
  iy<-match(inMap[i,2], sam_common)
  mx<-sort(rk1_p[ix,])
  my<-sort(rk1_p[,iy])
  if(rk1_p[ix, iy]<=mx[ntop] & rk1_p[ix, iy]<=my[ntop]) { ### top at both directions
  if(rs[ix, iy] >=corr_cutoff_self) { ### correlation cutoff 
    ## comparison with 
    outMap[ix, 1] <-sam_common[ix]
    outMap[ix, 2] <-sam_common[iy]
    outMap[ix, 3] <-1
  }
  }
}
### check all unmapped samples
for (i in 1:dim(rk1_p)[1]) {
  if (outMap[i,1] == 'NONE') {
    ix<-match(sam_type1_all[i], sam_type1_all)
    iy<-sort(rk1_p[ix,], index.return = TRUE)$ix[1]
    indic <- (iy > dim(inMap)[1])  
    if(!indic) { indic2 <- is.na(match(inMap[iy, 2], outMap[,2]))} else { indic2 <- TRUE}
    if(indic2) { 
    if(ix== sort(rk1_p[,iy], index.return = TRUE)$ix[1] & rk1_p[ix, iy]<=p_cutoff & rs[ix, iy] >= corr_cutoff) { ### top at both directions
      ## comparison with
      outMap[ix, 1] <-sam_type1_all[ix]
      outMap[ix, 2] <-sam_type2_all[iy]
      outMap[ix, 3] <-1
    }
  }
  }
}

### RECHECK input MAP
ntop=5;
for (i in 1:dim(inMap)[1]) {
  if (is.na(match(inMap[i,1], outMap[,1])) & is.na(match(inMap[i,2], outMap[,2]))) {
#      cat('i=', i, '\n')
      ix<-match(inMap[i,1], sam_common)
      iy<-match(inMap[i,2], sam_common)
      mx<-sort(rk1_p[ix,])
      my<-sort(rk1_p[,iy])
      if(rk1_p[ix, iy]<=mx[ntop] & rk1_p[ix, iy]<=my[ntop]) { ### top at both directions
        ## comparison with 
        outMap[ix, 1] <-sam_common[ix]
        outMap[ix, 2] <-sam_common[iy]
        outMap[ix, 3] <-1
      }
  }
}

###missed pairs

#### prepare for the next iteration
#cat( 'iteration #', iter,  'nT=', nT, ', nF=', nF,'\n')
map<-outMap[which(outMap[,3]==1),]
new <- paste(map[,1], map[,2], sep="-")
iter= iter+1;
}

###output map
map<-outMap[which(outMap[,3]==1),1:2]
colnames(map) <- c(type1, type2)




