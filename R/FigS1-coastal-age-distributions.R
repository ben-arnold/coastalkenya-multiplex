#--------------------------------
# coastal-age-distributions.R
#
# summarize age distribution
# by community
#--------------------------------

#--------------------------------
# preamble
#--------------------------------
rm(list=ls())

#--------------------------------
# load the formatted dataset
#--------------------------------
load("~/dropbox/coastalkenya/data/final/coastal_kenya.RData")
d <- coastal_kenya


pdf("~/dropbox/coastalkenya/results/figs/coastal-age-dists.pdf",width=5,height=10)
lo <- layout(mat=matrix(1:12,nrow=6,ncol=2,byrow=T))
op <- par(mar=c(4,4,2,1)+0.1)

# each village
vnames <- unique(d$cname[order(d$community)])
for(i in 1:10) {
  hist(d$age[d$community==i],
       breaks=0:100,
       main="",
       xlim=c(0,100),ylim=c(0,20),
       col="gray70",border=NA,
       xlab="",ylab="",
       las=1)
  mtext(paste(vnames[i]," (n= ",length(d$age[d$community==i]),")",sep=""),side=3,line=0.5,at=0,adj=0)
  if(i %in% c(1,3,5,7,9)) {
    mtext("N",side=2,line=3,las=1,cex=0.75) 
  }
  if(i %in% c(10)) {
    mtext("Age, years",side=1,line=2.5,las=1,cex=0.75) 
  }

}

# all villages
hist(d$age,
     breaks=0:100,
     main="",
     xlim=c(0,100),ylim=c(0,140),
     col="gray70",border=NA,
     xlab="",ylab="",
     las=1)
mtext(paste("All communities (n= ",length(d$age),")",sep=""),side=3,line=0.5,at=0,adj=0)
mtext("N",side=2,line=3,las=1,cex=0.75) 
mtext("Age, years",side=1,line=2.5,las=1,cex=0.75)

par(op)
layout(matrix(1))
dev.off()
