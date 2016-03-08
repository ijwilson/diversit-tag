## Install required packages
if (!require(data.table)) {
  install.packages("data.table")
  require(data.table)
}
if (!require("colorspace")) {
  install.packages("colorspace")
  require(colorspace)
}
if (!require("stringdist")) {
  install.packages("stringdist")
  require(stringdist)
}

merge_by_tag <- function(d1, d2) {
  da <- merge(d1, d2, by="tag", all=TRUE)
  da[is.na(da)] <- 0
  da
}

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

#############################
estimate_G <- function(count, count_baseline, min_baseline=0, min_p=0.0) {
  rawestimatep <- function(p, p_baseline) {
    mean((p-p_baseline)^2)/mean(p_baseline*(1-p_baseline))
  } 
  if (missing(count_baseline)) {
    p <- 1/length(count)
    p_b <- count/sum(count)
    return(mean((p_b-p)^2)*length(count))
  } else {
    pbase <- count_baseline/sum(count_baseline)
    u <- count_baseline>=min_baseline & pbase >= min_p
    pc <- count[u]/sum(count[u])
    p_b <- count_baseline[u]/sum(count_baseline[u])
    return(rawestimatep(pc, p_b))
  }
}
##############################################
tagdistplot <- function(tags, f, howmany=5) {
  maxtags <- tags[order(f, decreasing = TRUE)][1:howmany]
  tag1closest <- rep("none", length(tags))
  labels=c(c("1st", "2nd", "3rd"), paste(4:20,"th", sep=""))
  for (i in seq(howmany,1,-1)) {
    d <- stringdist(maxtags[i], tags)
    tag1closest[d==1] <- labels[i]
  }
  require(ggplot2)
  u <- d>0 & f>0
  
  pl <- ggplot(data.frame(x=d[u], y=f[u], closest=factor(tag1closest[u])),aes(x=x,y=y, col=closest)) 
  pl <- pl + geom_jitter(width=0.05, alpha=0.8) + scale_y_log10("Frequency of Tag") 
  pl <- pl + scale_x_continuous("Distance from Most Frequent Tag", breaks=c(1,4,8,12))
  pl + scale_colour_brewer(palette="Dark2")
}
##################################################
simulate_drifttags <- function(n_tags, total_depth, drift_G, alpha=0.02) {
  rdirichlet <- function(k, beta) {
    a <- rgamma(k, beta, 1.0)
    a/sum(a)
  }
  background_p <- rdirichlet(n_tags, alpha)
  background_count <- rpois(n_tags, total_depth*background_p)
  driftdistribution <- background_p*(1/drift_G-1)
  p <- rdirichlet(n_tags, driftdistribution)
  count <- rpois(n_tags, total_depth*p)
  data.frame(count=count, background_count=background_count, p=p, background_p = background_p)
}
##################################################
