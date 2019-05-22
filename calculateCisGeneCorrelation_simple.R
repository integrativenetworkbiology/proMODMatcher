#### to calculate correlation between pairs of cis genes
calculateCisGeneCorrelation_simple <- function  (prof1, prof2, map) {


### common genes
gene_common <- intersect(rownames(prof1), rownames(prof2))

### sample weight
wt<-as.numeric(map[,3])

###### data matrices
m1 <- as.matrix(prof1[match(gene_common, rownames(prof1)), match(map[,1], colnames(prof1))])
m2 <- as.matrix(prof2[match(gene_common, rownames(prof2)), match(map[,2], colnames(prof2))])
ngene <- dim(m1)[1]
nsample <- dim(m1)[2]

#### to find cis gene pairs
wtf<-TRUE
if (!wtf) {
  cat('unweighted version\n')
cis<-sapply(seq.int(ngene),function(i) cor.test(as.numeric(m1[i,]),as.numeric(m2[i,]),method="spearman",exact=T))
#cis<-sapply(seq.int(ngene),function(i) cor.test(as.numeric(m1[i,]),as.numeric(m2[i,]),exact=T))
cis_rho<-sapply(seq.int(ngene), function(i) cis[,i]$estimate)
cis_pva<-sapply(seq.int(ngene), function(i) cis[,i]$p.value)
cis_qva<-p.adjust(cis_pva,method="BH",n=length(cis_pva))
}
else {
##### use weighted correlation
M1<-t(sapply(seq.int(dim(m1)[1]), function(i) rank(m1[i,])))
M2<-t(sapply(seq.int(dim(m2)[1]), function(i) rank(m2[i,])))
cis_rho<-sapply (seq.int(ngene), function(i) cov.wt(t(rbind(M1[i,], M2[i,])), wt=wt, cor=TRUE)$cor[1,2])

d<-nsample
r<-cis_rho
cis_pva<- 2*pt(-abs(r)*sqrt((d-2)/(1-r^2)), d-2)
cis_qva<-p.adjust(cis_pva,method="BH",n=length(cis_pva))
}
return (list(gene_common, cis_rho, cis_pva, cis_qva))
}



calculateCisGeneCorrelation_simple_wtf <- function  (prof1, prof2, map, wtf) {


### common genes
gene_common <- intersect(rownames(prof1), rownames(prof2))

### sample weight
wt<-as.numeric(map[,3])

###### data matrices
m1 <- as.matrix(prof1[match(gene_common, rownames(prof1)), match(map[,1], colnames(prof1))])
m2 <- as.matrix(prof2[match(gene_common, rownames(prof2)), match(map[,2], colnames(prof2))])
ngene <- dim(m1)[1]
nsample <- dim(m1)[2]

#### to find cis gene pairs
#wtf<-TRUE
if (!wtf) {
  cat('unweighted version\n')
cis<-sapply(seq.int(ngene),function(i) cor.test(as.numeric(m1[i,]),as.numeric(m2[i,]),method="pearson"))
cis_rho<-sapply(seq.int(ngene), function(i) cis[,i]$estimate)
cis_pva<-sapply(seq.int(ngene), function(i) cis[,i]$p.value)
cis_qva<-p.adjust(cis_pva,method="BH",n=length(cis_pva))
}
else {
##### use weighted correlation
M1<-t(sapply(seq.int(dim(m1)[1]), function(i) rank(m1[i,])))
M2<-t(sapply(seq.int(dim(m2)[1]), function(i) rank(m2[i,])))
cis_rho<-sapply (seq.int(ngene), function(i) cov.wt(t(rbind(M1[i,], M2[i,])), wt=wt, cor=TRUE)$cor[1,2])

d<-nsample
r<-cis_rho
cis_pva<- 2*pt(-abs(r)*sqrt((d-2)/(1-r^2)), d-2)
cis_qva<-p.adjust(cis_pva,method="BH",n=length(cis_pva))
}
return (list(gene_common, cis_rho, cis_pva, cis_qva))
}

