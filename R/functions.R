
if (!require(data.table)) {
  install.packages("data.table")
  require(data.table)
}
if (!require("colorspace")) {
  install.packages("colorspace")
  require(colorspace)
}
## recursively merge a set of reads that are indexed by tag
mergeall <- function(mm) {
  merge1 <- function(ll) {
    n <- length(ll)
    res <- list()
    for (ii in seq(2,n,2)) {
      res[[ii%/%2]] <- merge(ll[[ii-1]], ll[[ii]], all=TRUE)
    }  
    if (length(ll)%%2==1) {
      res[[length(res)+1]] <- ll[[n]]
    }
    res
  }
  while(length(mm)>1) {
    mm <- merge1(mm)
  }
  mm[[1]]
}

merge_by_tag <- function(d1, d2) {
  da <- merge(d1, d2, by="tag", all=TRUE)
  da[is.na(da)] <- 0
  da
}

setup <- function(experimentnames) {
  l <- lapply(experimentnames,
          function(x){
            load(paste("rda/",x,".rda",sep=""));
            get(x);
          }
  )
  m <- mergeall(l)
  m[is.na(m)] <- 0
  m
}

#############################################################################################
rawestimatep <- function(p, p_baseline) {
  mean((p-p_baseline)^2)/mean(p_baseline*(1-p_baseline))
} 
#############################################################################################
D2 <- function(p1, p2) {
##  http://evolution.genetics.washington.edu/phylip/doc/gendist.html
  (sum((p1-p2)^2))/(sqrt(sum(p1^2))*sqrt(sum(p2^2)))
}

FF <- function(p1, p2) {
  (sum((p1-p2)^2))/(2*(1-sum(p1*p2)))
} 


estimateD2 <- function(count1, count2, min_count=0) {
  u <- count1 >= min_count | count2 > min_count
  p1 <- count1[u]/sum(count1[u])
  p2 <- count2[u]/sum(count2[u])
  
  D2(p1, p2)
}


estimateFF <- function(count1, count2, min_count=0) {
  u <- count1 >= min_count | count2 > min_count
  p1 <- count1[u]/sum(count1[u])
  p2 <- count2[u]/sum(count2[u])
  
  FF(p1, p2)
}
#############################
estimatecount <- function(count2, count_baseline, min_baseline=0, min_p=0) {
  pbase <- count_baseline/sum(count_baseline)
  u <- count_baseline>=min_baseline & pbase >= min_p
  p <- count2[u]/sum(count2[u])
  p_b <- count_baseline[u]/sum(count_baseline[u])
  rawestimatep(p, p_b)
}
#############################
estimatecount_equalfreq <- function(count2) {
  p <- 1/length(count2)
  p_b <- count2/sum(count2)
  mean((p_b-p)^2)*length(count2)
}
#############################
estimatecount2 <- function(count2, count_baseline, min_baseline=0, min_p=0) {
  pbase <- count_baseline/sum(count_baseline)
  u <- count_baseline>=min_baseline & pbase >= min_p
  p <- count2[u]/sum(count2[u])
  p_b <- count_baseline[u]/sum(count_baseline[u])
  mean((p-p_b)^2)/mean(p_b*(1-p_b)) - 1/sum(count2) - 1/sum(count_baseline)
}
## Code run to get the experiment names
experiments <- scan("all.files", sep="\n",what=character(0))
experiments <- gsub(pattern = "Mouse ","",experiments)
experiments <- gsub(pattern = ".xlsx","",experiments)
ex <- experiments




plot_tagfreq <- function(m, n_colours=12, minp=0.005) {
  require(colorspace)
  ptmp <- rowSums(m)/sum(m)
  u <- ptmp > minp
  p <- sweep(m[u,], 2, colSums(m[u,]), "/")
  o <- order(p[,1], decreasing = TRUE)
  p <- p[o,][1:n_colours,]
  p <- rbind(p, 1-colSums(p))
  
  colours <-c(rainbow_hcl(n_colours, c=90, l=70), "lightgrey")
  barplot(as.matrix(p), col=colours, beside=FALSE)
}


