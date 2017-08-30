

#--------------------------------
# coastal-mixture-distributions.R
#
# estimate seropositivity cutoffs
# from continuous Ab levels using
# gaussian mixture models
#
# select the number of groups
# data-adaptively using AIC loss
# for distributions that are
# really wide -- malaria
#--------------------------------

#--------------------------------
# preamble
#--------------------------------
rm(list=ls())
library(mixtools)


#--------------------------------
# fn to determin optimal k
# using AIC 
# (not used)
#--------------------------------

optmix <- function(x,lambda,k) {
  # estimate mixtures for values of k
  nk <- length(k)
  estk <- list()
  for(i in 1:nk) {
    cat("\n k = ",k[i],"\n")
    estk[[i]] <- normalmixEM(x=x[!is.na(x)],lambda=lambda,k=k[i])
    summary(estk[[i]])
  }
  
  # compute AIC
  LLs <-  sapply(estk, function(y) y$loglik)
  AICs <- sapply(estk, function(y) {
    return(2*(length(c(y$lambda,y$mu,y$sigma))-1) - 2*y$loglik)
    })
  cat("\nLogLik and AICs for each value k:\n")
  print(cbind(k,LLs,AICs))
  
  # return optimal fit
  optfit <- estk[which(AICs==min(AICs))][[1]]
  cat("\nOptimal k in terms of AIC:",length(optfit$mu),"\n\n")
  return(optfit)
}


#--------------------------------
# load the formatted dataset
#--------------------------------
load("~/dropbox/coastalkenya/data/final/coastal_kenya.RData")
d <- coastal_kenya

#--------------------------------
# Malaria
#--------------------------------
set.seed(1234)
cspmix <- normalmixEM(d$csp,lambda=0.5,k=2)
  summary(cspmix)
pfmix <- normalmixEM(d$msp1pf,lambda=0.5,k=2)
  summary(pfmix)
pmmix <- normalmixEM(d$msp1pm,lambda=0.5,k=2)
  summary(pmmix)

cspmixcut <- (cspmix$mu+3*cspmix$sigma)[1]
  cspmixcut
pfmixcut <- (pfmix$mu+3*pfmix$sigma)[1]
  pfmixcut
pmmixcut <- (pmmix$mu+3*pmmix$sigma)[1]
  pmmixcut

#--------------------------------
# LF
#--------------------------------

set.seed(1234)
wb123mix <- normalmixEM(d$wb123,lambda=0.5,k=2)
  summary(wb123mix)
bm14mix <- normalmixEM(d$bm14,lambda=0.5,k=2)
  summary(bm14mix)
bm33mix <- normalmixEM(d$bm33,lambda=0.5,k=2)
  summary(bm33mix)
  
wb123mixcut <- (wb123mix$mu+3*wb123mix$sigma)[1]
  wb123mixcut
bm14mixcut <- (bm14mix$mu+3*bm14mix$sigma)[1]
  bm14mixcut
bm33mixcut <- (bm33mix$mu+3*bm33mix$sigma)[1]
  bm33mixcut

#--------------------------------
# Schisto
#--------------------------------
  
set.seed(1234)
seamix  <- normalmixEM(d$sea,lambda=0.5,k=2)
  summary(seamix)
sm25mix <- normalmixEM(d$sm25,lambda=0.5,k=2)
  summary(sm25mix)
  
seamixcut <- (seamix$mu+3*seamix$sigma)[1]
  seamixcut
sm25mixcut <- (sm25mix$mu+3*sm25mix$sigma)[1]
  sm25mixcut
  
#--------------------------------
# Strongy + Ascaris
#--------------------------------
set.seed(1234)
niemix  <- normalmixEM(d$nie,lambda=0.5,k=2)
  summary(niemix)
ascarismix <- normalmixEM(d$ascaris,lambda=0.5,k=2)
  summary(ascarismix)
  
niemixcut <- (niemix$mu+3*niemix$sigma)[1]
  niemixcut
ascarismixcut <- (ascarismix$mu+3*ascarismix$sigma)[1]
  ascarismixcut
  
#--------------------------------
# save results
#--------------------------------
save.image(file="~/dropbox/coastalkenya/results/raw/coastal-mixtures.RData")


