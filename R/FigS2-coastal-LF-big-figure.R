

#--------------------------------
# coastal-LF-figure-big-figure.R
#
# plot LF Ab levels and seroprev
# by age E(Y_ax) and community
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
load(file="~/dropbox/coastalkenya/results/raw/coastal-LF.RData")

# order all villages by county and from south to north
vordn <- c("Kimorigo","Makwenyeni","Mirihini","Mwadimu","Kinarani","Jaribuni","Masindeni", "Mikinduni","Kipini","Ndau")
vord <- c(2,7,6,5,4,1,9,8,3,10)

#--------------------------------
# open a PDF graphics window
# and set layout
#--------------------------------
pdf("~/dropbox/coastalkenya/results/figs/coastal-means-prev-LF.pdf",height=13,width=18)
op <- par(mar=c(5,5,5,1)+0.1,xpd=TRUE)
lo <- layout(mat=matrix(c(1,2,3,14, 4,5,6,13, 7,8,9,15, 10,11,12,16),nrow=4,ncol=4),widths=rep(1,4),heights=c(1,1,1,0.2))

#--------------------------------
# plot community-stratified curves
# antibody titers
#--------------------------------
xtics <- seq(0,60,by=10)
ytics <- seq(1,4,by=1)

# Wb123
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=c(1,4.25),bty="n")
axis(1,at=xtics,las=1,cex.axis=1.5)
axis(2,at=ytics,labels=c(
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)), las=1,cex.axis=1.25)
mtext("Age, years",side=1,line=3,cex=1.1)
mtext("Antibody response (MFI-bg)",side=2,line=3,cex=1.1)
mtext("Age-dependent mean response, by community",side=3,line=2.5,adj=0.5,cex=1.1)
mtext('Lymphatic filariasis Wb123',side=3,line=0.5,adj=0,cex=1.1,font=1)
mtext("A",side=3,line=2.5,adj=0,at=-10,cex=1.5,font=2)


j <- 1
for(i in vord){ 
  lines(wb123_curves[[i]]$Age,wb123_curves[[i]]$pY,col=cols[j])
  j <- j+1
}


# Bm14
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=c(1,4.25),bty="n")
axis(1,at=xtics,las=1,cex.axis=1.5)
axis(2,at=ytics,labels=c(
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)), las=1,cex.axis=1.25)
# mtext("Age, years",side=1,line=3)
mtext("Antibody response (MFI-bg)",side=2,line=3)
# mtext("Age-dependent mean response, by community",side=3,line=0.5,adj=0,cex=1)
mtext('Lymphatic filariasis Bm14',side=3,line=0.5,adj=0,cex=1.1,font=1)

j <- 1
for(i in vord){ 
  lines(bm14_curves[[i]]$Age,bm14_curves[[i]]$pY,col=cols[j])
  j <- j+1
}


# Bm33
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=c(1,4.25),bty="n")
axis(1,at=xtics,las=1,cex.axis=1.5)
axis(2,at=ytics,labels=c(
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)), las=1,cex.axis=1.25)
mtext("Age, years",side=1,line=3,cex=1.1)
mtext("Antibody response (MFI-bg)",side=2,line=3,cex=1.1)
# mtext("Age-dependent mean response, by community",side=3,line=0.5,adj=0,cex=1)
mtext('Lymphatic filariasis Bm33',side=3,line=0.5,adj=0,cex=1.1,font=1)

j <- 1
for(i in vord){ 
  lines(bm33_curves[[i]]$Age,bm33_curves[[i]]$pY,col=cols[j])
  j <- j+1
}

#--------------------------------
# plot community-level means
#--------------------------------

xtics <- 1:10

# Wb123
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=c(1,4.25),bty="n")
axis(2,at=ytics,labels=c(
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)), las=1,cex.axis=1.25)
mtext(1:10,side=1,line=0,at=1:10,col=cols)
mtext("Mean response, by community",side=3,line=2.5,adj=0.5,cex=1.1)
# mtext("Community",side=1,line=3)


arrows(x0=1:10,y0=wb123EYxs[2,vord],y1=wb123EYxs[3,vord],angle=90,col=cols,code=3,length=0.05)
points(1:10,wb123EYxs[1,vord],pch=19,col=cols)


# Bm14
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=c(1,4.25),bty="n")
axis(2,at=ytics,labels=c(
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)), las=1,cex.axis=1.25)
mtext(1:10,side=1,line=0,at=1:10,col=cols)
# mtext("Mean response, by community",side=3,line=0.5,adj=0)
# mtext("Community",side=1,line=3)

arrows(x0=1:10,y0=bm14EYxs[2,vord],y1=bm14EYxs[3,vord],angle=90,col=cols,code=3,length=0.05)
points(1:10,bm14EYxs[1,vord],pch=19,col=cols)


# Bm33
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=c(1,4.25),bty="n")
axis(2,at=ytics,labels=c(
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)), las=1,cex.axis=1.25)
mtext(1:10,side=1,line=0,at=1:10,col=cols)
# mtext("Mean response, by community",side=3,line=0.5,adj=0)
mtext("Community",side=1,line=3,cex=1.1)

arrows(x0=1:10,y0=bm33EYxs[2,vord],y1=bm33EYxs[3,vord],angle=90,col=cols,code=3,length=0.05)
points(1:10,bm33EYxs[1,vord],pch=19,col=cols)

#--------------------------------
# add community labels
#--------------------------------
# plot(1,1,type="n",yaxt="n",xaxt="n",ylab="",xlab="",bty="n")
# plot(1,1,type="n",yaxt="n",xaxt="n",ylab="",xlab="",bty="n")
# op <- par(xpd=TRUE)
# plot(1,1,type="n",yaxt="n",xaxt="n",ylab="",xlab="",bty="n",xlim=c(1,3),ylim=c(0,1.1))  
# text(x=rep(1.1,10),y=seq(1.1,0,length=10),1:10,col=cols,adj=1,cex=1.5)
# text(x=rep(1.3,10),y=seq(1.1,0,length=10),paste("community name",1:10),adj=0,cex=1.5)


#--------------------------------
# plot seroprevalence curves means
#--------------------------------
xtics <- seq(0,60,by=10)
ytics <- seq(0,1,by=0.1)

# Wb123
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=range(ytics),bty="n")
axis(1,at=xtics,las=1,cex.axis=1.5)
axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100), las=1,cex.axis=1.25)
mtext("Age-dependent seroprevalence, by community",side=3,line=2.5,adj=0.5,cex=1.1)
# mtext("Age, years",side=1,line=3,cex=1.1)
mtext("Seroprevalence (%)",side=2,line=3,cex=1.1)
mtext('Lymphatic filariasis Wb123',side=3,line=0.5,adj=0,cex=1.1,font=1)
mtext("B",side=3,line=2.5,adj=0,at=-15,cex=1.5,font=2)



j <- 1
for(i in vord){ 
  lines(pwb123_curves[[i]]$Age,pwb123_curves[[i]]$pY,col=cols[j])
  j <- j+1
}


# Bm14
xtics <- seq(0,60,by=10)
ytics <- seq(0,1,by=0.1)
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=range(ytics),bty="n")
axis(1,at=xtics,las=1,cex.axis=1.5)
axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100), las=1,cex.axis=1.25)
# mtext("Age-dependent seroprevalence, by community",side=3,line=0.5,adj=0,cex=1)
# mtext("Age, years",side=1,line=3)
mtext("Seroprevalence (%)",side=2,line=3,cex=1.1)
mtext('Lymphatic filariasis Bm14',side=3,line=0.5,adj=0,cex=1.1,font=1)

j <- 1
for(i in vord){ 
  lines(pbm14_curves[[i]]$Age,pbm14_curves[[i]]$pY,col=cols[j])
  j <- j+1
}


# Bm33
xtics <- seq(0,60,by=10)
ytics <- seq(0,1,by=0.1)
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=range(ytics),bty="n")
axis(1,at=xtics,las=1,cex.axis=1.5)
axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100), las=1,cex.axis=1.25)
# mtext("Age-dependent seroprevalence, by community",side=3,line=0.5,adj=0,cex=1)
mtext("Age, years",side=1,line=3,cex=1.1)
mtext("Seroprevalence (%)",side=2,line=3,cex=1.1)
mtext('Lymphatic filariasis Bm33',side=3,line=0.5,adj=0,cex=1.1,font=1)

j <- 1
for(i in vord){ 
  lines(pbm33_curves[[i]]$Age,pbm33_curves[[i]]$pY,col=cols[j])
  j <- j+1
}

#--------------------------------
# plot community-level means
#--------------------------------
xtics <- 1:10

# Wb123
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=range(ytics),bty="n")
axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100), las=1,cex.axis=1.25)
mtext(1:10,side=1,line=0,at=1:10,col=cols)
mtext("Seroprevalence, by community",side=3,line=2.5,adj=0.5,cex=1.1)
# mtext("Community",side=1,line=3)

arrows(x0=1:10,y0=wb123_prev[vord,5],y1=wb123_prev[vord,6],angle=90,col=cols,code=3,length=0.05)
points(1:10,wb123_prev[vord,4],pch=19,col=cols)


# Bm14
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=range(ytics),bty="n")
axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100), las=1,cex.axis=1.25)
mtext(1:10,side=1,line=0,at=1:10,col=cols)
# mtext("Seroprevalence, by community",side=3,line=0.5,adj=0)
# mtext("Community",side=1,line=3)

arrows(x0=1:10,y0=bm14_prev[vord,5],y1=bm14_prev[vord,6],angle=90,col=cols,code=3,length=0.05)
points(1:10,bm14_prev[vord,4],pch=19,col=cols)


# P. malariae
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=range(ytics),bty="n")
axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100), las=1,cex.axis=1.25)
mtext(1:10,side=1,line=0,at=1:10,col=cols)
# mtext("Seroprevalence, by community",side=3,line=0.5,adj=0)
mtext("Community",side=1,line=3,cex=1.1)

arrows(x0=1:10,y0=bm33_prev[vord,5],y1=bm33_prev[vord,6],angle=90,col=cols,code=3,length=0.05)
points(1:10,bm33_prev[vord,4],pch=19,col=cols)

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




