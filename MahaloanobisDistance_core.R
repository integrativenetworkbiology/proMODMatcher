## row= type1, col=type2
## distance with sign : if the correlation itself is smaller than mean -> takes negative 

distance_self <- function(xx,corr_mat){
	re <- list()
        tmp <- cbind(corr_mat[xx,], corr_mat[,xx])
        tt <- mahalanobis(tmp, colMeans(tmp), cov(tmp))
        sign_tt <- apply(as.matrix(1:length(tt)), 1, function(x){ if( corr_mat[xx,x] < colMeans(tmp)[1] | corr_mat[x,xx]< colMeans(tmp)[2]) { -tt[x]} else { tt[x]}})
        re$distance <- sign_tt[xx]
	re$order <- length(sign_tt)-rank(sign_tt)[xx]+1
	re$p.value <- 1-pchisq(tt[xx], df=1)
        re
}


distance_Others <- function(xx, corr_mat, nsample){
        re <- list()
        if(dim(corr_mat)[1]==nsample){
                others <- seq.int(nsample)
                dis_others <- apply(as.matrix(others), 1, function(pp) {tmp <- cbind(corr_mat[xx,], corr_mat[, pp])
                                                                  tmp <- tmp[-unique(c(xx,pp)),]
                                                                  tmp <- rbind(tmp, c(corr_mat[xx,pp], corr_mat[xx, pp]))
                                                                  tt <- mahalanobis(tmp, colMeans(tmp), cov(tmp))
                                                                  if(corr_mat[xx,pp] < colMeans(tmp)[1]| corr_mat[xx,pp] < colMeans(tmp)[2]) { dis <- -tt[length(tt)] } else { dis <- tt[length(tt)]}
                                                                  dis})
                re$distance <- dis_others
                re$p.value <- 1-pchisq(abs(dis_others), df=1)
                re
	} else {
		print("non-common samples")
	}
}


## with mRNA of xx and other's miRNA 
distance_Others_type1 <- function(xx, corr_mat, nsample){
	re <- list()
        others <- setdiff(seq.int(dim(corr_mat)[2]), seq.int(nsample))
        dis_others <- apply(as.matrix(others), 1, function(pp) { tmp <- cbind(corr_mat[xx, 1:nsample], corr_mat[1:nsample,pp])
                                                                 tmp <- rbind(tmp, c(corr_mat[xx,pp], corr_mat[xx,pp]))
                                                                 tt <- mahalanobis(tmp, colMeans(tmp), cov(tmp))
                                                                 if(corr_mat[xx,pp] < colMeans(tmp)[1]| corr_mat[xx,pp] < colMeans(tmp)[2]) { dis <- -tt[length(tt)] } else { dis <- tt[length(tt)]}
                                                                 dis})
                re$distance <- dis_others
                re$p.value <- 1-pchisq(abs(dis_others), df=1)
                re

}

## with miRNA of xx and other's mRNA 
distance_Others_type2 <- function(xx, corr_mat, nsample){
## with miRNA of xx and other's mRNA
#        others <- setdiff(sam, xx)
#        others <- setdiff(seq.int(nsample), xx)
	re <- list()
        others <- setdiff(seq.int(dim(corr_mat)[1]), seq.int(nsample))
        dis_others <- apply(as.matrix(others), 1, function(pp) {tmp <- cbind(corr_mat[1:nsample,xx], corr_mat[pp,1:nsample])
                                                          tmp <- rbind(tmp, c(corr_mat[pp,xx], corr_mat[pp,xx]))
                                                          tt <- mahalanobis(tmp, colMeans(tmp), cov(tmp))
                                                          if(corr_mat[pp,xx] < colMeans(tmp)[1]|corr_mat[pp,xx] < colMeans(tmp)[2]) { dis <- -tt[length(tt)] } else { dis <- tt[length(tt)]}
                                                          dis})
	re$distance <- dis_others
        re$p.value <- 1-pchisq(abs(dis_others), df=1)
        re
}


#####multivariate normal distribution 
require(mnormt)
mnorm_self <- function(xx,corr_mat){
        re <- list()
        tmp <- cbind(corr_mat[xx,], corr_mat[,xx])
        tt <- pmnorm(tmp[xx,], colMeans(tmp), cov(tmp))
	tt_others <- apply(tmp, 1, function(x)  pmnorm(x, colMeans(tmp), cov(tmp)))
        re$prob <- tt[[1]]
        re$order <- length(tt_others)-rank(tt_others)[xx]+1
        re
}

mnorm_self_sim <- function(xx,corr_mat){
        re <- list()
        tmp <- cbind(corr_mat[xx,], corr_mat[,xx])
        tt <- pmnorm(tmp[xx,], colMeans(tmp), cov(tmp))
        re$prob <- tt[[1]]
        re
}

mnorm_Others <- function(xx, corr_mat, nsample){
        re <- list()
#        if(dim(corr_mat)[1]==nsample){
                others <- seq.int(nsample)
                prob_others <- apply(as.matrix(others), 1, function(pp) {mm <- min(dim(corr_mat)[1], dim(corr_mat)[2])
								  tmp <- cbind(corr_mat[xx,1:mm], corr_mat[1:mm, pp])
                                                                  tmp <- tmp[-unique(c(xx,pp)),]
                                                                  tmp <- rbind(tmp, c(corr_mat[xx,pp], corr_mat[xx, pp]))
                                                                  tt <- pmnorm(tmp[dim(tmp)[1],], colMeans(tmp), cov(tmp))[[1]]
                                                                  tt})
                re$prob <- prob_others
                re
#        } else {
#                print("non-common samples")
#        }
}



mnorm_Others_moreSamples <- function(xx, corr_mat, nsample, corr_mat_permute){
        re <- list()
        if(dim(corr_mat)[1]==nsample){
                others <- seq.int(nsample)
                prob_others <- apply(as.matrix(others), 1, function(pp) {tmp <- cbind(corr_mat[xx,], corr_mat[, pp])
                                                                  tmp <- tmp[-unique(c(xx,pp)),]
                                                                  tmp <- rbind(tmp, c(corr_mat[xx,pp], corr_mat[xx, pp]))
								  tmp <- rbind(cbind(corr_mat_permute[xx,], corr_mat_permute[,pp]), tmp)
                                                                  tt <- pmnorm(tmp[dim(tmp)[1],], colMeans(tmp), cov(tmp))[[1]]
                                                                  tt})
                re$prob <- prob_others
                re
        } else {
                print("non-common samples")
        }
}






mnorm_Others_type1 <- function(xx, corr_mat, nsample){
# with mRNA of xx and other's miRNA
#        others <- setdiff(sam, xx)
#        others <- setdiff(seq.int(nsample), xx)
        re <- list()
        others <- setdiff(seq.int(dim(corr_mat)[2]), seq.int(nsample))
        prob_others <- apply(as.matrix(others), 1, function(pp) { tmp <- cbind(corr_mat[xx, 1:nsample], corr_mat[1:nsample,pp])
                                                                 tmp <- rbind(tmp, c(corr_mat[xx,pp], corr_mat[xx,pp]))
                                                                 tt <- pmnorm(tmp[dim(tmp)[1],], colMeans(tmp), cov(tmp))[[1]]
                                                                 tt})
        re$prob <- prob_others
        re
}

## with miRNA of xx and other's mRNA 
mnorm_Others_type2 <- function(xx, corr_mat, nsample){
## with miRNA of xx and other's mRNA
#        others <- setdiff(sam, xx)
#        others <- setdiff(seq.int(nsample), xx)
        re <- list()
        others <- setdiff(seq.int(dim(corr_mat)[1]), seq.int(nsample))
        prob_others <- apply(as.matrix(others), 1, function(pp) {tmp <- cbind(corr_mat[1:nsample,xx], corr_mat[pp,1:nsample])
                                                          tmp <- rbind(tmp, c(corr_mat[pp,xx], corr_mat[pp,xx]))
                                                          tt <- pmnorm(tmp[dim(tmp)[1],], colMeans(tmp), cov(tmp))[[1]]
                                                          tt})
        re$prob <- prob_others
        re
}


mnorm_Others_type1_moreSamples <- function(xx, corr_mat, nsample, corr_mat_permute){
# with mRNA of xx and other's miRNA
#        others <- setdiff(sam, xx)
#        others <- setdiff(seq.int(nsample), xx)
        re <- list()
        others <- setdiff(seq.int(dim(corr_mat)[2]), seq.int(nsample))
        prob_others <- apply(as.matrix(others), 1, function(pp) { tmp <- cbind(corr_mat[xx, 1:nsample], corr_mat[1:nsample,pp])
                                                                 tmp <- rbind(tmp, c(corr_mat[xx,pp], corr_mat[xx,pp]))
								 tmp <- rbind(cbind(corr_mat_permute[xx,], corr_mat_permute[,pp]), tmp)
                                                                 tt <- pmnorm(tmp[dim(tmp)[1],], colMeans(tmp), cov(tmp))[[1]]
                                                                 tt})
        re$prob <- prob_others
        re
}

mnorm_Others_type2_moreSamples <- function(xx, corr_mat, nsample, corr_mat_permute){
## with miRNA of xx and other's mRNA
#        others <- setdiff(sam, xx)
#        others <- setdiff(seq.int(nsample), xx)
        re <- list()
        others <- setdiff(seq.int(dim(corr_mat)[1]), seq.int(nsample))
        prob_others <- apply(as.matrix(others), 1, function(pp) {tmp <- cbind(corr_mat[1:nsample,xx], corr_mat[pp,1:nsample])
                                                          tmp <- rbind(tmp, c(corr_mat[pp,xx], corr_mat[pp,xx]))
							  tmp <- rbind(cbind(corr_mat_permute[,xx], corr_mat_permute[pp,]), tmp)
                                                          tt <- pmnorm(tmp[dim(tmp)[1],], colMeans(tmp), cov(tmp))[[1]]
                                                          tt})
        re$prob <- prob_others
        re
}



