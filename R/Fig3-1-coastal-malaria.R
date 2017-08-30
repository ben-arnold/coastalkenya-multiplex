

#--------------------------------
# coastal-malaria.R
#
# summarize malaria Ab levels
# and seroprevalence
# by age E(Y_ax) and community
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
d$pcsp    <- ifelse(d$csp>cspmixcut,1,0)
d$pmsp1pf <- ifelse(d$msp1pf>pfmixcut,1,0)
d$pmsp1pm <- ifelse(d$msp1pm>pmmixcut,1,0)

#--------------------------------
# set negative and zero values to
# 1 before the log transform
#--------------------------------
# csp (all >0)
table(d$msp1pf<=0)
d["msp1pf"][d["msp1pf"]<=0] <-1
# msp1pm (all >0)

#--------------------------------
# summarize antibody curves
# by community
#--------------------------------

# ensemble library
SL.library <-  c("SL.mean", "SL.glm", "SL.gam", "SL.loess")


set.seed(1234)
csp_curves <- lapply(1:10,
                       FUN = function(x) {
                         agecurveAb(Y=log10(d$csp[d$community==x]),
                                    Age=d$age[d$community==x],
                                    SL.library=SL.library,
                                    gamdf=2:6
                                    )
                       }
                       )

set.seed(1234)
msp1pf_curves <- lapply(1:10,
                      FUN = function(x) {
                        agecurveAb(Y=log10(d$msp1pf[d$community==x]),
                                   Age=d$age[d$community==x],
                                   SL.library=SL.library,
                                   gamdf=2:6
                        )
                      }
                      )

set.seed(1234)
msp1pm_curves <- lapply(1:10,
                      FUN = function(x) {
                        agecurveAb(Y=log10(d$msp1pm[d$community==x]),
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

# CSP
set.seed(1234)
cspEYx <- lapply(1:10,FUN = function(x) {
  tmleAb(Y=log10(d$csp[d$community==x]),
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library,gamdf=2:6)
})

# Pf MSP-1
set.seed(1234)
msp1pfEYx <- lapply(1:10,FUN = function(x) {
  tmleAb(Y=log10(d$msp1pf[d$community==x]),
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library,gamdf=2:6)
})

# Pm MSP-1
set.seed(1234)
msp1pmEYx <- lapply(1:10,FUN = function(x) {
  tmleAb(Y=log10(d$msp1pm[d$community==x]),
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

cspEYxs     <- getpsi(cspEYx)
msp1pfEYxs  <- getpsi(msp1pfEYx)
msp1pmEYxs  <- getpsi(msp1pmEYx)



#--------------------------------
# seroprevalence curves
# by community
#--------------------------------


set.seed(1234)
pcsp_curves <- lapply(1:10,
                     FUN = function(x) {
                       agecurveAb(Y=d$pcsp[d$community==x],
                                  Age=d$age[d$community==x],
                                  SL.library=SL.library,
                                  gamdf=2:6
                       )
                     }
                     )

set.seed(1234)
pmsp1pf_curves <- lapply(1:10,
                        FUN = function(x) {
                          agecurveAb(Y=d$pmsp1pf[d$community==x],
                                     Age=d$age[d$community==x],
                                     SL.library=SL.library,
                                     gamdf=2:6
                          )
                        }
                        )

set.seed(1234)
pmsp1pm_curves <- lapply(1:10,
                        FUN = function(x) {
                          agecurveAb(Y=d$pmsp1pm[d$community==x],
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

# CSP
set.seed(1234)
pcspEYx <- lapply(1:10,FUN = function(x) {
  tmleAb(Y=d$pcsp[d$community==x],
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library,gamdf=2:6)
         })

# MSP-1 Pf
set.seed(1234)
pmsp1pfEYx <- lapply(1:10,FUN = function(x) {
  tmleAb(Y=d$pmsp1pf[d$community==x],
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library,gamdf=2:6)
  })

# MSP-1 Pm
set.seed(1234)
pmsp1pmEYx <- lapply(1:10,FUN = function(x) {
  tmleAb(Y=d$pmsp1pm[d$community==x],
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library,gamdf=2:6)
  })


#--------------------------------
# calculate differences in means
# between Jaribuni and other
# communities in Kalifi county
# and Ndau and Kipini compared with
# others for P-value reporting
# differences are all highly significant
# based on estimates and SEs, but this
# is just a way to get a formal P-value
# (requested by co-authors)
#--------------------------------

# Jaribuni
jaribuni01 <- ifelse(d$cname=='Jaribuni',1,0)
jdiffcsppf <- tmleAb(Y=d$pcsp[d$county=='Kilifi'],
                      X=jaribuni01[d$county=='Kilifi'],
                      W=data.frame(Age=d[d$county=='Kilifi',c("age")]),
                      SL.library=SL.library,gamdf=2:6)

jdiffmsp1pf <- tmleAb(Y=d$pmsp1pf[d$county=='Kilifi'],
                     X=jaribuni01[d$county=='Kilifi'],
                     W=data.frame(Age=d[d$county=='Kilifi',c("age")]),
                     SL.library=SL.library,gamdf=2:6)

jdiffmsp1pf <- tmleAb(Y=d$pmsp1pm[d$county=='Kilifi'],
                      X=jaribuni01[d$county=='Kilifi'],
                      W=data.frame(Age=d[d$county=='Kilifi',c("age")]),
                      SL.library=SL.library,gamdf=2:6)

# Ndau
ndau01 <- ifelse(d$cname=='Ndau',1,0)
ndiffcsppf <- tmleAb(Y=d$pcsp,
                     X=ndau01,
                     W=data.frame(Age=d[c("age")]),
                     SL.library=SL.library,gamdf=2:6)

ndiffmsp1pf <- tmleAb(Y=d$pmsp1pf,
                      X=ndau01,
                      W=data.frame(Age=d[c("age")]),
                      SL.library=SL.library,gamdf=2:6)

ndiffmsp1pf <- tmleAb(Y=d$pmsp1pm,
                      X=ndau01,
                      W=data.frame(Age=d[c("age")]),
                      SL.library=SL.library,gamdf=2:6)
# Kipini
kipini01 <- ifelse(d$cname=='Kipini',1,0)
kdiffcsppf <- tmleAb(Y=d$pcsp,
                     X=kipini01,
                     W=data.frame(Age=d[c("age")]),
                     SL.library=SL.library,gamdf=2:6)

kdiffmsp1pf <- tmleAb(Y=d$pmsp1pf,
                      X=kipini01,
                      W=data.frame(Age=d[c("age")]),
                      SL.library=SL.library,gamdf=2:6)

kdiffmsp1pf <- tmleAb(Y=d$pmsp1pm,
                      X=kipini01,
                      W=data.frame(Age=d[c("age")]),
                      SL.library=SL.library,gamdf=2:6)


#--------------------------------
# pull out means, lb and ub
# from the TMLE fitted objects
#--------------------------------
pcspEYxs   <- getpsi(pcspEYx)
pmsp1pfEYxs <- getpsi(pmsp1pfEYx)
pmsp1pmEYxs <- getpsi(pmsp1pmEYx)


#--------------------------------
# save results
#--------------------------------
save.image(file="~/dropbox/coastalkenya/results/raw/coastal-malaria.RData")


