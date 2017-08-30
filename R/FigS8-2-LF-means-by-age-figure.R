

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
library(scales)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# bright color blind palette:  https://personal.sron.nl/~pault/ 
cblack <- "#000004FF"
cblue <- "#3366AA"
cteal <- "#11AA99"
cgreen <- "#66AA55"
cchartr <- "#CCCC55"
cmagent <- "#992288"
cred <- "#EE3333"
corange <- "#EEA722"
cyellow <- "#FFEE33"
cgrey <- "#777777"

cols <- c(cblack,"gray40",cred,corange,cchartr,cgreen,cteal,cblue,cmagent, cbPalette[8])



#--------------------------------
# load saved results
#--------------------------------
# load(file="~/dropbox/coastalkenya/results/raw/coastal-malaria.RData")
load(file="~/dropbox/coastalkenya/results/raw/coastal-LF-by-age.RData")

# order all villages in this figure
# by CSP seroprevalence
# vord <- order(pcspEYxs[1,],decreasing=TRUE)
# vordn <- unique(d$cname[order(d$community)])[vord]

# order all villages by county and from south to north
vordn <- c("Kimorigo","Makwenyeni","Mirihini","Mwadimu","Kinarani","Jaribuni","Masindeni", "Mikinduni","Kipini","Ndau")
vord <- c(2,7,6,5,4,1,9,8,3,10)

# summarize the sample size in each community and age stratum
childNs <- sapply(ictmus,FUN=function(x) x[1,])
min(childNs[1:2,])
max(childNs[1:2,])
mean(childNs[1:2,])
quantile(childNs[1:2,],probs=c(0.25,0.5,0.75))

#--------------------------------
# open a PDF graphics window
# and set layout
#--------------------------------

pdf("~/dropbox/coastalkenya/results/figs/coastal-LF-by-age.pdf",height=12,width=9)
op <- par(mar=c(5,5,5,1)+0.1,xpd=TRUE)
lo <- layout(mat=matrix(c(1,2,3,4,9, 5,6,7,8,9),nrow=5,ncol=2),
             widths=rep(1,2),heights=c(1,1,1,1,0.3))

xtics <- 1:21
ytics <- seq(0,0.6,by=0.1)

#----------------------------
# ICT prev
#----------------------------

plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=range(ytics),bty="n")
axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100),las=1,cex.axis=1.25)
mtext("Seroprevalence",side=3,line=0.5,adj=0)
mtext("Seroprevalence (%)",side=2,line=3)
mtext("ICT",side=3,line=2,font=2,adj=0,cex=1.25)

# data
arrows(x0=1:10,y0=ict_lb[vord,1],y1=ict_ub[vord,1],angle=90,col=cols,code=3,length=0.05)
points(1:10,ict_psi[vord,1],pch=19,col=cols)

arrows(x0=12:21,y0=ict_lb[vord,2],y1=ict_ub[vord,2],angle=90,col=cols,code=3,length=0.05)
points(12:21,ict_psi[vord,2],pch=19,col=cols)

# labels
segments(x0=11,y0=-0.1,y1=0.6)
mtext(1:10,at=c(1:10,12:21),col=cols,side=1,line=0,cex=0.8)
mtext(c("Ages 2-5","Ages 6-10"),at=c(mean(1:10),mean(12:21)),side=1,line=2)

#----------------------------
# Wb123 seroprev
#----------------------------

  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
       xlim=range(xtics),ylim=range(ytics),bty="n")
  axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100),las=1,cex.axis=1.25)
  mtext("Seroprevalence",side=3,line=0.5,adj=0)
  mtext("Seroprevalence (%)",side=2,line=3)
  mtext("Wb123",side=3,line=2,font=2,adj=0,cex=1.25)
  
  # data
  arrows(x0=1:10,y0=pwb123_lb[vord,1],y1=pwb123_ub[vord,1],angle=90,col=cols,code=3,length=0.05)
  points(1:10,pwb123_psi[vord,1],pch=19,col=cols)
  
  arrows(x0=12:21,y0=pwb123_lb[vord,2],y1=pwb123_ub[vord,2],angle=90,col=cols,code=3,length=0.05)
  points(12:21,pwb123_psi[vord,2],pch=19,col=cols)
  
  # labels
  segments(x0=11,y0=-0.1,y1=0.6)
  mtext(1:10,at=c(1:10,12:21),col=cols,side=1,line=0,cex=0.8)
  mtext(c("Ages 2-5","Ages 6-10"),at=c(mean(1:10),mean(12:21)),side=1,line=2)
  
#----------------------------
# Bm14 seroprev
#----------------------------

  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
       xlim=range(xtics),ylim=range(ytics),bty="n")
  axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100),las=1,cex.axis=1.25)
  mtext("Seroprevalence",side=3,line=0.5,adj=0)
  mtext("Seroprevalence (%)",side=2,line=3)
  mtext("Bm14",side=3,line=2,font=2,adj=0,cex=1.25)
  
  # data
  arrows(x0=1:10,y0=pbm14_lb[vord,1],y1=pbm14_ub[vord,1],angle=90,col=cols,code=3,length=0.05)
  points(1:10,pbm14_psi[vord,1],pch=19,col=cols)
  
  arrows(x0=12:21,y0=pbm14_lb[vord,2],y1=pbm14_ub[vord,2],angle=90,col=cols,code=3,length=0.05)
  points(12:21,pbm14_psi[vord,2],pch=19,col=cols)
  
  # labels
  segments(x0=11,y0=-0.1,y1=0.6)
  mtext(1:10,at=c(1:10,12:21),col=cols,side=1,line=0,cex=0.8)
  mtext(c("Ages 2-5","Ages 6-10"),at=c(mean(1:10),mean(12:21)),side=1,line=2)
  
  #----------------------------
  # Bm33 seroprev
  #----------------------------
  
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
       xlim=range(xtics),ylim=range(ytics),bty="n")
  axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100),las=1,cex.axis=1.25)
  mtext("Seroprevalence",side=3,line=0.5,adj=0)
  mtext("Seroprevalence (%)",side=2,line=3)
  mtext("Bm33",side=3,line=2,font=2,adj=0,cex=1.25)
  
  # data
  arrows(x0=1:10,y0=pbm33_lb[vord,1],y1=pbm33_ub[vord,1],angle=90,col=cols,code=3,length=0.05)
  points(1:10,pbm33_psi[vord,1],pch=19,col=cols)
  
  arrows(x0=12:21,y0=pbm33_lb[vord,2],y1=pbm33_ub[vord,2],angle=90,col=cols,code=3,length=0.05)
  points(12:21,pbm33_psi[vord,2],pch=19,col=cols)
  
  # labels
  segments(x0=11,y0=-0.1,y1=0.6)
  mtext(1:10,at=c(1:10,12:21),col=cols,side=1,line=0,cex=0.8)
  mtext(c("Ages 2-5","Ages 6-10"),at=c(mean(1:10),mean(12:21)),side=1,line=2)


  #--------------------------------
  # plot age- and community- stratified means
  #--------------------------------
  ytics <- 1:3
  
  #----------------------------
  # ICT (blank plot)
  #----------------------------
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
       xlim=range(xtics),ylim=range(ytics),bty="n")
  # mtext("ICT",side=3,line=2,font=2,adj=0,cex=1.25)
  
  
  #----------------------------
  # Wb123 means
  #----------------------------
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
       xlim=range(xtics),ylim=range(ytics),bty="n")
  axis(2,at=ytics,labels=c(
    expression(10^1),
    expression(10^2),
    expression(10^3)), las=1,cex.axis=1.25)
  mtext("Mean response",side=3,line=0.5,adj=0)
  mtext("Antibody response (MFI-bg)",side=2,line=3)
  # mtext("Wb123",side=3,line=2,font=2,adj=0,cex=1.25)
  
  # data
  arrows(x0=1:10,y0=wb123_lb[vord,1],y1=wb123_ub[vord,1],angle=90,col=cols,code=3,length=0.05)
  points(1:10,wb123_psi[vord,1],pch=19,col=cols)
  
  arrows(x0=12:21,y0=wb123_lb[vord,2],y1=wb123_ub[vord,2],angle=90,col=cols,code=3,length=0.05)
  points(12:21,wb123_psi[vord,2],pch=19,col=cols)
  
  # labels
  segments(x0=11,y0=0.6,y1=3)
  mtext(1:10,at=c(1:10,12:21),col=cols,side=1,line=0,cex=0.8)
  mtext(c("Ages 2-5","Ages 6-10"),at=c(mean(1:10),mean(12:21)),side=1,line=2)
  
  #----------------------------  
  # Bm14 means
  #----------------------------
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
       xlim=range(xtics),ylim=range(ytics),bty="n")
  axis(2,at=ytics,labels=c(
    expression(10^1),
    expression(10^2),
    expression(10^3)), las=1,cex.axis=1.25)
  mtext("Mean response",side=3,line=0.5,adj=0)
  mtext("Antibody response (MFI-bg)",side=2,line=3)
  # mtext("Bm14",side=3,line=2,font=2,adj=0,cex=1.25)
  
  # data
  arrows(x0=1:10,y0=bm14_lb[vord,1],y1=bm14_ub[vord,1],angle=90,col=cols,code=3,length=0.05)
  points(1:10,bm14_psi[vord,1],pch=19,col=cols)
  
  arrows(x0=12:21,y0=bm14_lb[vord,2],y1=bm14_ub[vord,2],angle=90,col=cols,code=3,length=0.05)
  points(12:21,bm14_psi[vord,2],pch=19,col=cols)
  
  # labels
  segments(x0=11,y0=0.6,y1=3)
  mtext(1:10,at=c(1:10,12:21),col=cols,side=1,line=0,cex=0.8)
  mtext(c("Ages 2-5","Ages 6-10"),at=c(mean(1:10),mean(12:21)),side=1,line=2)
  
  #----------------------------
  # Bm33 means
  #----------------------------
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
       xlim=range(xtics),ylim=range(ytics),bty="n")
  axis(2,at=ytics,labels=c(
    expression(10^1),
    expression(10^2),
    expression(10^3)), las=1,cex.axis=1.25)
  mtext("Mean response",side=3,line=0.5,adj=0)
  mtext("Antibody response (MFI-bg)",side=2,line=3)
  # mtext("Bm33",side=3,line=2,font=2,adj=0,cex=1.25)
  
  # data
  arrows(x0=1:10,y0=bm33_lb[vord,1],y1=bm33_ub[vord,1],angle=90,col=cols,code=3,length=0.05)
  points(1:10,bm33_psi[vord,1],pch=19,col=cols)
  
  arrows(x0=12:21,y0=bm33_lb[vord,2],y1=bm33_ub[vord,2],angle=90,col=cols,code=3,length=0.05)
  points(12:21,bm33_psi[vord,2],pch=19,col=cols)
  
  # labels
  segments(x0=11,y0=0.6,y1=3)
  mtext(1:10,at=c(1:10,12:21),col=cols,side=1,line=0,cex=0.8)
  mtext(c("Ages 2-5","Ages 6-10"),at=c(mean(1:10),mean(12:21)),side=1,line=2)
  
  

  
  
  #--------------------------------
  # add a community name key along the bottom
  #--------------------------------
  # resetplot <- function() {
  #   layout(mat=matrix(1))
  #   par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  #   plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
  # }
  # resetplot()
  # legend('bottom',legend=paste(1:10,vordn),text.col=cols,ncol=5,bty="n",xpd=NA)
  
  
  #--------------------------------
  # add a community name key along the bottom
  # alternative format, by county
  #--------------------------------
  resetplot <- function() {
    layout(mat=matrix(1))
    par(mfrow=c(1, 1), oma=rep(0, 4), mar=c(3,0,0,0), new=TRUE)
    plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
  }
  resetplot()
  legxs <- c(0,0.2,0.4,0.6,0.8)+0.05
  mtext(c("Taita Taveta","Kwale","Kilifi","Tana River","Lamu"),side=1,line=-1,at=legxs,adj=0,font=2)
  mtext(paste(1,vordn[1]),side=1,line=0,adj=0,at=legxs[1],col=cols[1])
  
  mtext(paste(2,vordn[2]),side=1,line=0,adj=0,at=legxs[2],col=cols[2])
  mtext(paste(3,vordn[3]),side=1,line=1,adj=0,at=legxs[2],col=cols[3])
  mtext(paste(4,vordn[4]),side=1,line=2,adj=0,at=legxs[2],col=cols[4])
  
  mtext(paste(5,vordn[5]),side=1,line=0,adj=0,at=legxs[3],col=cols[5])
  mtext(paste(6,vordn[6]),side=1,line=1,adj=0,at=legxs[3],col=cols[6])
  mtext(paste(7,vordn[7]),side=1,line=2,adj=0,at=legxs[3],col=cols[7])
  
  mtext(paste(8,vordn[8]),side=1,line=0,adj=0,at=legxs[4],col=cols[8])
  mtext(paste(9,vordn[9]),side=1,line=1,adj=0,at=legxs[4],col=cols[9])
  
  mtext(paste(10,vordn[10]),side=1,line=0,adj=0,at=legxs[5],col=cols[10])
  
  par(op)
  layout(mat=matrix(1))
  dev.off()
  

  