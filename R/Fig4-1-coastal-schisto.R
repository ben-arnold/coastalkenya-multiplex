

#--------------------------------
# coastal-schisto.R
#
# summarize schisto Ab results
# by community
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
# results for alternate LF cutoffs
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
d$psea  <- ifelse(d$sea>seamixcut,1,0)
d$psm25 <- ifelse(d$sm25>sm25mixcut,1,0)

#--------------------------------
# set negative and zero values to
# 1 before the log transform
#--------------------------------
d["sea"][d["sea"]<=0] <-1
d["sm25"][d["sm25"]<=0] <-1

#--------------------------------
# summarize antibody curves
# by community
#--------------------------------

# ensemble library
SL.library <-  c("SL.mean", "SL.glm", "SL.gam", "SL.loess")


set.seed(1234)
sea_curves <- lapply(unique(d$community),
                     FUN = function(x) {
                       agecurveAb(Y=log10(d$sea[d$community==x]),
                                  Age=d$age[d$community==x],
                                  SL.library=SL.library,
                                  gamdf=2:6
                       )
                     }
)

set.seed(1234)
sm25_curves <- lapply(unique(d$community),
                        FUN = function(x) {
                          agecurveAb(Y=log10(d$sm25[d$community==x]),
                                     Age=d$age[d$community==x],
                                     SL.library=SL.library,
                                     gamdf=2:6
                          )
                        }
)


#--------------------------------
# Estmiate village-level mean 
# Ab levels EY_x
#--------------------------------

# SEA
set.seed(1234)
seaEYx <- lapply(unique(d$community),FUN = function(x) {
  tmleAb(Y=log10(d$sea[d$community==x]),
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library,gamdf=2:6)
})

# Sm25
set.seed(1234)
sm25EYx <- lapply(unique(d$community),FUN = function(x) {
  tmleAb(Y=log10(d$sm25[d$community==x]),
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library,gamdf=2:6)
})


#--------------------------------
# pull out means, lb and ub
# from the TMLE fitted objects
#--------------------------------
getpsi <- function(x) {
  ests <- sapply(x, function(xx) c(xx$psi,xx$lb,xx$ub))
  rownames(ests) <- c("psi","lb","ub") 
  colnames(ests) <- 1:10
  return(ests)
}

seaEYxs <- getpsi(seaEYx)
sm25EYxs  <- getpsi(sm25EYx)


#--------------------------------
# seroprevalence curves
# by community
#--------------------------------

set.seed(1234)
psea_curves <- lapply(unique(d$community),
                       FUN = function(x) {
                         agecurveAb(Y=d$psea[d$community==x],
                                    Age=d$age[d$community==x],
                                    SL.library=SL.library,
                                    gamdf=2:6
                         )
                       }
)

set.seed(1234)
psm25_curves <- lapply(unique(d$community),
                      FUN = function(x) {
                        agecurveAb(Y=d$psm25[d$community==x],
                                   Age=d$age[d$community==x],
                                   SL.library=SL.library,
                                   gamdf=2:6
                        )
                      }
)


#--------------------------------
# fn to estimate prev + exact 95% cis
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

#--------------------------------
# SEA seroprevalence by community
#--------------------------------
seap <- tapply(d$psea,d$community,exactprev)
seap <- matrix(unlist(seap),nrow=10,ncol=5,byrow=T)
sea_prev <- data.frame(community=1:10,seap)
names(sea_prev) <- c("community","N","n","prev","min95","max95")
round(sea_prev,3)

#--------------------------------
# Sm25 seroprevalence by community
#--------------------------------
sm25p <- tapply(d$psm25,d$community,exactprev)
sm25p <- matrix(unlist(sm25p),nrow=10,ncol=5,byrow=T)
sm25_prev <- data.frame(community=1:10,sm25p)
names(sm25_prev) <- c("community","N","n","prev","min95","max95")
round(sm25_prev,3)

#--------------------------------
# save results
#--------------------------------
save.image(file="~/dropbox/coastalkenya/results/raw/coastal-schisto.RData")


