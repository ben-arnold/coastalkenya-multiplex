

#--------------------------------
# coastal-LF.R
#
# summarize LF Ab and ICT results
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
# summarize antibody curves
# by community
#--------------------------------

# ensemble library
SL.library <-  c("SL.mean", "SL.glm", "SL.gam", "SL.loess")

set.seed(1234)
wb123_curves <- lapply(unique(d$community),
                     FUN = function(x) {
                       agecurveAb(Y=log10(d$wb123[d$community==x]),
                                  Age=d$age[d$community==x],
                                  SL.library=SL.library,
                                  gamdf=2:6
                       )
                     }
)

set.seed(1234)
bm14_curves <- lapply(unique(d$community),
                        FUN = function(x) {
                          agecurveAb(Y=log10(d$bm14[d$community==x]),
                                     Age=d$age[d$community==x],
                                     SL.library=SL.library,
                                     gamdf=2:6
                          )
                        }
)

set.seed(1234)
bm33_curves <- lapply(unique(d$community),
                        FUN = function(x) {
                          agecurveAb(Y=log10(d$bm33[d$community==x]),
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

# Wb123
set.seed(1234)
wb123EYx <- lapply(unique(d$community),FUN = function(x) {
  tmleAb(Y=log10(d$wb123[d$community==x]),
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library,gamdf=2:6)
})

# Bm14
set.seed(1234)
bm14EYx <- lapply(unique(d$community),FUN = function(x) {
  tmleAb(Y=log10(d$bm14[d$community==x]),
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library,gamdf=2:6)
})

# Bm33
set.seed(1234)
bm33EYx <- lapply(unique(d$community),FUN = function(x) {
  tmleAb(Y=log10(d$bm33[d$community==x]),
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

wb123EYxs <- getpsi(wb123EYx)
bm14EYxs  <- getpsi(bm14EYx)
bm33EYxs  <- getpsi(bm33EYx)


#--------------------------------
# seroprevalence curves
# by community
#--------------------------------

set.seed(1234)
pwb123_curves <- lapply(unique(d$community),
                       FUN = function(x) {
                         agecurveAb(Y=d$pwb123[d$community==x],
                                    Age=d$age[d$community==x],
                                    SL.library=SL.library,
                                    gamdf=2:6
                         )
                       }
)

set.seed(1234)
pbm14_curves <- lapply(unique(d$community),
                      FUN = function(x) {
                        agecurveAb(Y=d$pbm14[d$community==x],
                                   Age=d$age[d$community==x],
                                   SL.library=SL.library,
                                   gamdf=2:6
                        )
                      }
)

set.seed(1234)
pbm33_curves <- lapply(unique(d$community),
                      FUN = function(x) {
                        agecurveAb(Y=d$pbm33[d$community==x],
                                   Age=d$age[d$community==x],
                                   SL.library=SL.library,
                                   gamdf=2:6
                        )
                      }
)


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


#--------------------------------
# ICT prevalence by community
#--------------------------------
ictp <- tapply(d$ict,d$community,exactprev)
ictp <- matrix(unlist(ictp),nrow=10,ncol=5,byrow=T)
ict_prev <- data.frame(community=1:10,ictp)
names(ict_prev) <- c("community","N","n","prev","min95","max95")

round(ict_prev,3)

#--------------------------------
# Wb123 seroprevalence by community
#--------------------------------
wb123p <- tapply(d$pwb123,d$community,exactprev)
wb123p <- matrix(unlist(wb123p),nrow=10,ncol=5,byrow=T)
wb123_prev <- data.frame(community=1:10,wb123p)
names(wb123_prev) <- c("community","N","n","prev","min95","max95")
round(wb123_prev,3)


#--------------------------------
# Bm14 seroprevalence by community
#--------------------------------
bm14p <- tapply(d$pbm14,d$community,exactprev)
bm14p <- matrix(unlist(bm14p),nrow=10,ncol=5,byrow=T)
bm14_prev <- data.frame(community=1:10,bm14p)
names(bm14_prev) <- c("community","N","n","prev","min95","max95")
round(bm14_prev,3)


#--------------------------------
# Bm33 seroprevalence by community
#--------------------------------
bm33p <- tapply(d$pbm33,d$community,exactprev)
bm33p <- matrix(unlist(bm33p),nrow=10,ncol=5,byrow=T)
bm33_prev <- data.frame(community=1:10,bm33p)
names(bm33_prev) <- c("community","N","n","prev","min95","max95")
round(bm33_prev,3)


#--------------------------------
# calculate differences in means
# between Jaribuni and other
# communities in Kalifi county
# for P-value reporting
# differences are all highly significant
# based on estimates and SEs, but this
# is just a way to get a formal P-value
# (requested by co-authors)
#--------------------------------

# Jaribuni
jaribuni01 <- ifelse(d$cname=='Jaribuni',1,0)
jdiffwb123 <- tmleAb(Y=d$pwb123[d$county=='Kilifi'],
                     X=jaribuni01[d$county=='Kilifi'],
                     W=data.frame(Age=d[d$county=='Kilifi',c("age")]),
                     SL.library=SL.library,gamdf=2:6)

jdiffbm14 <- tmleAb(Y=d$pbm14[d$county=='Kilifi'],
                      X=jaribuni01[d$county=='Kilifi'],
                      W=data.frame(Age=d[d$county=='Kilifi',c("age")]),
                      SL.library=SL.library,gamdf=2:6)

jdiffbm33 <- tmleAb(Y=d$pbm33[d$county=='Kilifi'],
                      X=jaribuni01[d$county=='Kilifi'],
                      W=data.frame(Age=d[d$county=='Kilifi',c("age")]),
                      SL.library=SL.library,gamdf=2:6)


#--------------------------------
# save results
#--------------------------------
save.image(file="~/dropbox/coastalkenya/results/raw/coastal-LF.RData")


