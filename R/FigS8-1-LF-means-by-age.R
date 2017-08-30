

#--------------------------------
# xx-coastal-LF-means-by-age.R
#
# summarize LF Ab by community
# and by age group
#--------------------------------

#--------------------------------
# preamble
#--------------------------------
rm(list=ls())
library(SuperLearner)
library(tmle)
library(tmleAb)

#--------------------------------
# load the gaussian mixture model
# results for seropositive cutoffs
#--------------------------------
load(file="~/dropbox/coastalkenya/results/raw/coastal-mixtures.RData")


#--------------------------------
# load the formatted dataset
#--------------------------------
load("~/dropbox/coastalkenya/data/final/coastal_kenya.RData")
d <- coastal_kenya

#--------------------------------
# identify seropositive individuals
#--------------------------------
d$pwb123 <- ifelse(d$wb123>wb123mixcut,1,0)
d$pbm14 <- ifelse(d$bm14>bm14mixcut,1,0)
d$pbm33 <- ifelse(d$bm33>bm33mixcut,1,0)

#--------------------------------
# set negative and zero values to
# 1 before the log transform
#--------------------------------
d["wb123"][d["wb123"]<=0] <-1
d["bm14"][d["bm14"]<=0] <-1
d["bm33"][d["bm33"]<=0] <-1


#--------------------------------
# fn to estimate prev + exact 95% cis
# given the very low seroprevalence
# in some communities, particularly
# for ICT
#--------------------------------
exactprev <- function(x) {
  # x : a binary indicator of seropositive(1) vs. seronegative(0)
  tabx <- table(x)
  if(length(tabx)<2) {
    if(names(tabx)=="1") {
      tabx <- c(0,tabx)
    } else{
      tabx <- c(tabx,0)
    }
  } 
  estx <- binom.test(x=tabx[2],n=sum(tabx))
  res <- c(estx$parameter,estx$statistic,estx$estimate,estx$conf.int)
  names(res) <- c("N","n","prev","min95","max95")
  return(res)
}

d$agecat <- cut(d$age,breaks=c(0,5,10,15,20,30,40,50,100))


#--------------------------------
# Wb123 mean Ab and
# seroprevalence by community
#--------------------------------
set.seed(1234)
wb123mus <- list()
for(i in 1:10) {
  wb123mus[[i]] <- sapply(levels(d$agecat),
                           FUN = function(x) {
                             tmleAb(Y=log10(d$wb123[d$community==i & d$agecat==x]),W=d[d$community==i & d$agecat==x,"age"],SL.library=c("SL.mean","SL.glm","SL.gam"),gamdf=2:6)[c("psi","se","lb","ub")]}
  )
}

pwb123mus <- list()
for(i in 1:10) {
  pwb123mus[[i]] <- sapply(levels(d$agecat),
                         FUN = function(x) {
                           exactprev(x=d$pwb123[d$community==i & d$agecat==x])}
  )
}


#--------------------------------
# Bm14 mean Ab and
# seroprevalence by community
#--------------------------------
set.seed(1234)
bm14mus <- list()
for(i in 1:10) {
  bm14mus[[i]] <- sapply(levels(d$agecat),
                          FUN = function(x) {
                            tmleAb(Y=log10(d$bm14[d$community==i & d$agecat==x]),W=d[d$community==i & d$agecat==x,"age"],SL.library=c("SL.mean","SL.glm","SL.gam"),gamdf=2:6)[c("psi","se","lb","ub")]}
  )
}

pbm14mus <- list()
for(i in 1:10) {
  pbm14mus[[i]] <- sapply(levels(d$agecat),
                           FUN = function(x) {
                             exactprev(x=d$pbm14[d$community==i & d$agecat==x])}
  )
}



#--------------------------------
# Bm33 mean Ab and
# seroprevalence by community
#--------------------------------
set.seed(1234)
bm33mus <- list()
for(i in 1:10) {
  bm33mus[[i]] <- sapply(levels(d$agecat),
                         FUN = function(x) {
                           tmleAb(Y=log10(d$bm33[d$community==i & d$agecat==x]),W=d[d$community==i & d$agecat==x,"age"],SL.library=c("SL.mean","SL.glm","SL.gam"),gamdf=2:6)[c("psi","se","lb","ub")]}
  )
}

pbm33mus <- list()
for(i in 1:10) {
  pbm33mus[[i]] <- sapply(levels(d$agecat),
                          FUN = function(x) {
                            exactprev(x=d$pbm33[d$community==i & d$agecat==x])}
  )
}



#--------------------------------
# ICT prevalence by community
#--------------------------------
ictmus <- list()
for(i in 1:10) {
  ictmus[[i]] <- sapply(levels(d$agecat),
                          FUN = function(x) {
                            exactprev(x=d$ict[d$community==i & d$agecat==x])}
                        )
}


#--------------------------------
# pull out means, lb and ub
# from the TMLE fitted objects
#--------------------------------

# Wb123
wb123_psi <- matrix(NA,nrow=10,ncol=8)
wb123_lb <- wb123_psi
wb123_ub <- wb123_psi
for(i in 1:10) {
  wb123_psi[i,] <- mapply(function(x1,x2) wb123mus[[x1]][,x2]$psi,x1=i,x2=1:8)
  wb123_lb[i,] <- mapply(function(x1,x2) wb123mus[[x1]][,x2]$lb,x1=i,x2=1:8)
  wb123_ub[i,] <- mapply(function(x1,x2) wb123mus[[x1]][,x2]$ub,x1=i,x2=1:8)
}
rownames(wb123_psi) <- rownames(wb123_lb) <- rownames(wb123_ub) <- 1:10
colnames(wb123_psi) <- colnames(wb123_lb) <- colnames(wb123_ub) <- levels(d$agecat)


# Bm14
bm14_psi <- matrix(NA,nrow=10,ncol=8)
bm14_lb <- bm14_psi
bm14_ub <- bm14_psi
for(i in 1:10) {
  bm14_psi[i,] <- mapply(function(x1,x2) bm14mus[[x1]][,x2]$psi,x1=i,x2=1:8)
  bm14_lb[i,] <- mapply(function(x1,x2) bm14mus[[x1]][,x2]$lb,x1=i,x2=1:8)
  bm14_ub[i,] <- mapply(function(x1,x2) bm14mus[[x1]][,x2]$ub,x1=i,x2=1:8)
}
rownames(bm14_psi) <- rownames(bm14_lb) <- rownames(bm14_ub) <- 1:10
colnames(bm14_psi) <- colnames(bm14_lb) <- colnames(bm14_ub) <- levels(d$agecat)

# Bm33
bm33_psi <- matrix(NA,nrow=10,ncol=8)
bm33_lb <- bm33_psi
bm33_ub <- bm33_psi
for(i in 1:10) {
  bm33_psi[i,] <- mapply(function(x1,x2) bm33mus[[x1]][,x2]$psi,x1=i,x2=1:8)
  bm33_lb[i,] <- mapply(function(x1,x2) bm33mus[[x1]][,x2]$lb,x1=i,x2=1:8)
  bm33_ub[i,] <- mapply(function(x1,x2) bm33mus[[x1]][,x2]$ub,x1=i,x2=1:8)
}
rownames(bm33_psi) <- rownames(bm33_lb) <- rownames(bm33_ub) <- 1:10
colnames(bm33_psi) <- colnames(bm33_lb) <- colnames(bm33_ub) <- levels(d$agecat)

#--------------------------------
# pull out prev, lb and ub
# from the exactprev objects
#--------------------------------

# Wb123
pwb123_psi <- matrix(NA,nrow=10,ncol=8)
pwb123_lb <- wb123_psi
pwb123_ub <- wb123_psi
for(i in 1:10) {
  pwb123_psi[i,] <- mapply(function(x1,x2) pwb123mus[[x1]][3,x2],x1=i,x2=1:8)
  pwb123_lb[i,] <- mapply(function(x1,x2) pwb123mus[[x1]][4,x2],x1=i,x2=1:8)
  pwb123_ub[i,] <- mapply(function(x1,x2) pwb123mus[[x1]][5,x2],x1=i,x2=1:8)
}
rownames(pwb123_psi) <- rownames(pwb123_lb) <- rownames(pwb123_ub) <- 1:10
colnames(pwb123_psi) <- colnames(pwb123_lb) <- colnames(pwb123_ub) <- levels(d$agecat)

# Bm14
pbm14_psi <- matrix(NA,nrow=10,ncol=8)
pbm14_lb <- wb123_psi
pbm14_ub <- wb123_psi
for(i in 1:10) {
  pbm14_psi[i,] <- mapply(function(x1,x2) pbm14mus[[x1]][3,x2],x1=i,x2=1:8)
  pbm14_lb[i,] <- mapply(function(x1,x2) pbm14mus[[x1]][4,x2],x1=i,x2=1:8)
  pbm14_ub[i,] <- mapply(function(x1,x2) pbm14mus[[x1]][5,x2],x1=i,x2=1:8)
}
rownames(pbm14_psi) <- rownames(pbm14_lb) <- rownames(pbm14_ub) <- 1:10
colnames(pbm14_psi) <- colnames(pbm14_lb) <- colnames(pbm14_ub) <- levels(d$agecat)

# Bm33
pbm33_psi <- matrix(NA,nrow=10,ncol=8)
pbm33_lb <- wb123_psi
pbm33_ub <- wb123_psi
for(i in 1:10) {
  pbm33_psi[i,] <- mapply(function(x1,x2) pbm33mus[[x1]][3,x2],x1=i,x2=1:8)
  pbm33_lb[i,] <- mapply(function(x1,x2) pbm33mus[[x1]][4,x2],x1=i,x2=1:8)
  pbm33_ub[i,] <- mapply(function(x1,x2) pbm33mus[[x1]][5,x2],x1=i,x2=1:8)
}
rownames(pbm33_psi) <- rownames(pbm33_lb) <- rownames(pbm33_ub) <- 1:10
colnames(pbm33_psi) <- colnames(pbm33_lb) <- colnames(pbm33_ub) <- levels(d$agecat)

# ICT
ict_psi <- matrix(NA,nrow=10,ncol=8)
ict_lb <- ict_psi
ict_ub <- ict_psi
for(i in 1:10) {
  ict_psi[i,] <- mapply(function(x1,x2) ictmus[[x1]][3,x2],x1=i,x2=1:8)
  ict_lb[i,]  <- mapply(function(x1,x2) ictmus[[x1]][4,x2],x1=i,x2=1:8)
  ict_ub[i,]  <- mapply(function(x1,x2) ictmus[[x1]][5,x2],x1=i,x2=1:8)
}
rownames(ict_psi) <- rownames(ict_lb) <- rownames(ict_ub)  <- 1:10
colnames(ict_psi) <- colnames(ict_lb) <- colnames(ict_ub) <- levels(d$agecat)

#--------------------------------
# save results
#--------------------------------
save.image(file="~/dropbox/coastalkenya/results/raw/coastal-LF-by-age.RData")


