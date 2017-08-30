

#--------------------------------
# coastal-sth.R
#
# summarize STH Ab results
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
d$pnie  <- ifelse(d$nie>niemixcut,1,0)
d$pascaris <- ifelse(d$ascaris>ascarismixcut,1,0)

#--------------------------------
# set negative and zero values to
# 1 before the log transform
#--------------------------------
d["nie"][d["nie"]<=0] <-1
d["ascaris"][d["ascaris"]<=0] <-1

#--------------------------------
# summarize antibody curves
# by community
#--------------------------------

# ensemble library
SL.library <-  c("SL.mean", "SL.glm", "SL.gam", "SL.loess")

# Strongy NIE
set.seed(1234)
nie_curves <- lapply(unique(d$community),
                     FUN = function(x) {
                       agecurveAb(Y=log10(d$nie[d$community==x]),
                                  Age=d$age[d$community==x],
                                  SL.library=SL.library,
                                  gamdf=2:6
                       )
                     }
)

# Ascaris AsHb
set.seed(1234)
ascaris_curves <- lapply(unique(d$community),
                        FUN = function(x) {
                          agecurveAb(Y=log10(d$ascaris[d$community==x]),
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

# Strongy NIE
set.seed(1234)
nieEYx <- lapply(unique(d$community),FUN = function(x) {
  tmleAb(Y=log10(d$nie[d$community==x]),
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library,gamdf=2:6)
})

# Ascaris AsHb
set.seed(1234)
ascarisEYx <- lapply(unique(d$community),FUN = function(x) {
  tmleAb(Y=log10(d$ascaris[d$community==x]),
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

nieEYxs <- getpsi(nieEYx)
ascarisEYxs  <- getpsi(ascarisEYx)


#--------------------------------
# seroprevalence curves
# by community
#--------------------------------

set.seed(1234)
pnie_curves <- lapply(unique(d$community),
                       FUN = function(x) {
                         agecurveAb(Y=d$pnie[d$community==x],
                                    Age=d$age[d$community==x],
                                    SL.library=SL.library,
                                    gamdf=2:6
                         )
                       }
)

set.seed(1234)
pascaris_curves <- lapply(unique(d$community),
                      FUN = function(x) {
                        agecurveAb(Y=d$pascaris[d$community==x],
                                   Age=d$age[d$community==x],
                                   SL.library=SL.library,
                                   gamdf=2:6
                        )
                      }
)


#--------------------------------
# estimate mean seroprevalence, 
# stratified by community
#--------------------------------

# StrongyNIE
set.seed(1234)
pnieEYx <- lapply(1:10,FUN = function(x) {
  tmleAb(Y=d$pnie[d$community==x],
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library,gamdf=2:6)
})

# Ascaris AsHb
set.seed(1234)
pascarisEYx <- lapply(1:10,FUN = function(x) {
  tmleAb(Y=d$pascaris[d$community==x],
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library,gamdf=2:6)
})

#--------------------------------
# pull out means, lb and ub
# from the TMLE fitted objects
#--------------------------------
pnieEYxs   <- getpsi(pnieEYx)
pascarisEYxs <- getpsi(pascarisEYx)


#--------------------------------
# save results
#--------------------------------
save.image(file="~/dropbox/coastalkenya/results/raw/coastal-sth.RData")


