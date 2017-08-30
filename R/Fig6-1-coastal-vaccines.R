

#--------------------------------
# coastal-EYax-vaccines-vil.R
#
# summarize vaccine Ab levels
# by age E(Y_ax) and village
#--------------------------------

#--------------------------------
# preamble
#--------------------------------
rm(list=ls())
library(SuperLearner)
library(tmle)
library(tmleAb)

#--------------------------------
# load the formatted dataset
#--------------------------------
load("~/dropbox/coastalkenya/data/final/coastal_kenya.RData")
d <- coastal_kenya


#--------------------------------
# identify immunoprotected using
# cutoff values sent from jeff priest
#--------------------------------
d$pmeasles <- ifelse(d$measles > 178,1,0)
d$pdiptheria <- ifelse(d$diptheria > 4393,1,0)
d$ptetanus <- ifelse(d$tetanus > 118, 1, 0)

#--------------------------------
# estimate antibody mean curves
# by community
#--------------------------------

# ensemble library
SL.library <-  c("SL.mean", "SL.glm", "SL.gam", "SL.loess")

set.seed(1234)
mea_curves <- lapply(unique(d$community),
                       FUN = function(x) {
                         agecurveAb(Y=log10(d$measles[d$community==x]),
                                    Age=d$age[d$community==x],
                                    SL.library=SL.library,
                                    gamdf=2:6
                                    )
                       }
                       )

set.seed(1234)
dip_curves <- lapply(unique(d$community),
                      FUN = function(x) {
                        agecurveAb(Y=log10(d$diptheria[d$community==x]),
                                   Age=d$age[d$community==x],
                                   SL.library=SL.library,
                                   gamdf=2:6
                        )
                      }
                      )

set.seed(1234)
tet_curves <- lapply(unique(d$community),
                      FUN = function(x) {
                        agecurveAb(Y=log10(d$tetanus[d$community==x]),
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

# measles
set.seed(1234)
meaEYx <- lapply(unique(d$community),FUN = function(x) {
  tmleAb(Y=log10(d$measles[d$community==x]),
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library)
})

# diptheria
set.seed(1234)
dipEYx <- lapply(unique(d$community),FUN = function(x) {
  tmleAb(Y=log10(d$diptheria[d$community==x]),
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library)
})

# tetanus
set.seed(1234)
tetEYx <- lapply(unique(d$community),FUN = function(x) {
  tmleAb(Y=log10(d$tetanus[d$community==x]),
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library)
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
meaEYxs <- getpsi(meaEYx)
dipEYxs <- getpsi(dipEYx)
tetEYxs <- getpsi(tetEYx)


#--------------------------------
# estimate immunoprotection curves
# by community
#--------------------------------

set.seed(1234)
pmea_curves <- lapply(unique(d$community),
                      FUN = function(x) {
                        agecurveAb(Y=d$pmeasles[d$community==x],
                                   Age=d$age[d$community==x],
                                   SL.library=SL.library
                        )
                      }
)

set.seed(1234)
pdip_curves <- lapply(unique(d$community),
                      FUN = function(x) {
                        agecurveAb(Y=d$pdiptheria[d$community==x],
                                   Age=d$age[d$community==x],
                                   SL.library=SL.library
                        )
                      }
)

set.seed(1234)
ptet_curves <- lapply(unique(d$community),
                      FUN = function(x) {
                        agecurveAb(Y=d$ptetanus[d$community==x],
                                   Age=d$age[d$community==x],
                                   SL.library=SL.library
                        )
                      }
)

#--------------------------------
# estimate mean immunoprotection, 
# stratified by community
#--------------------------------

# measles
set.seed(1234)
pmeaEYx <- lapply(unique(d$community),FUN = function(x) {
  tmleAb(Y=d$pmeasles[d$community==x],
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library)
})

# diptheria
set.seed(1234)
pdipEYx <- lapply(unique(d$community),FUN = function(x) {
  tmleAb(Y=d$pdiptheria[d$community==x],
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library)
})

# tetanus
set.seed(1234)
ptetEYx <- lapply(unique(d$community),FUN = function(x) {
  tmleAb(Y=d$ptetanus[d$community==x],
         W=data.frame(Age=d[d$community==x,c("age")]),
         SL.library=SL.library)
})

#--------------------------------
# pull out means, lb and ub
# from the TMLE fitted objects
#--------------------------------
pmeaEYxs   <- getpsi(pmeaEYx)
pdipEYxs <- getpsi(pdipEYx)
ptetEYxs <- getpsi(ptetEYx)


#--------------------------------
# save results
#--------------------------------
save.image(file="~/dropbox/coastalkenya/results/raw/coastal-vaccines.RData")




