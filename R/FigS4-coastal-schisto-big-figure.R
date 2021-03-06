

#--------------------------------
# coastal-schisto-big-figure.R
#
# plot schisto Ab levels and seroprevalence
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
load(file="~/dropbox/coastalkenya/results/raw/coastal-schisto.RData")

# order all villages by county and from south to north
vordn <- c("Kimorigo","Makwenyeni","Mirihini","Mwadimu","Kinarani","Jaribuni","Masindeni", "Mikinduni","Kipini","Ndau")
vord <- c(2,7,6,5,4,1,9,8,3,10)

#--------------------------------
# open a PDF graphics window
# and set layout
#--------------------------------
pdf("~/dropbox/coastalkenya/results/figs/coastal-means-prev-schisto.pdf",height=9,width=18)
op <- par(mar=c(5,5,5,1)+0.1,xpd=TRUE)
lo <- layout(mat=matrix(c(1,2,10, 3,4,9, 5,6,11, 7,8,12),nrow=3,ncol=4),widths=rep(1,4),heights=c(1,1,0.2))

#--------------------------------
# plot community-stratified curves
# antibody titers
#--------------------------------
xtics <- seq(0,60,by=10)
ytics <- 1:4

# SEA
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=c(1,4.25),bty="n")
axis(1,at=xtics,las=1,cex.axis=1.5)
axis(2,at=ytics,labels=c(
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4)), las=1,cex.axis=1.25)
# mtext("Age, years",side=1,line=3)
mtext("Antibody response (MFI-bg)",side=2,line=3,cex=1.1)
mtext("Age-dependent mean response, by community",side=3,line=2.5,adj=0.5,cex=1.1)
mtext(expression(paste(italic("S. mansoni")," SEA")),side=3,line=0.5,adj=0,cex=1.1,font=1)
mtext("A",side=3,line=2.5,adj=0,at=-10,cex=1.5,font=2)

j <- 1
for(i in vord){ 
  lines(sea_curves[[i]]$Age,sea_curves[[i]]$pY,col=cols[j])
  j <- j+1
}


# Sm25
ytics <- 0:3
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=c(0,3.25),bty="n")
axis(1,at=xtics,las=1,cex.axis=1.5)
axis(2,at=ytics,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3)), las=1,cex.axis=1.25)
mtext("Age, years",side=1,line=3,cex=1.1)
mtext("Antibody response (MFI-bg)",side=2,line=3,cex=1.1)
# mtext("Age-dependent mean response, by community",side=3,line=0.5,adj=0,cex=1)
mtext(expression(paste(italic("S. mansoni")," Sm25")),side=3,line=0.5,adj=0,cex=1.1,font=1)

j <- 1
for(i in vord){ 
  lines(sm25_curves[[i]]$Age,sm25_curves[[i]]$pY,col=cols[j])
  j <- j+1
}


#--------------------------------
# plot community-level means
#--------------------------------

xtics <- 1:10
ytics <- 1:4
# SEA
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

arrows(x0=1:10,y0=seaEYxs[2,vord],y1=seaEYxs[3,vord],angle=90,col=cols,code=3,length=0.05)
points(1:10,seaEYxs[1,vord],pch=19,col=cols)


# Sm25
ytics <- 0:3
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=c(0,3.25),bty="n")
axis(2,at=ytics,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3)), las=1,cex.axis=1.25)
mtext(1:10,side=1,line=0,at=1:10,col=cols)
# mtext("Mean response, by community",side=3,line=0.5,adj=0)
mtext("Community",side=1,line=3,cex=1.1)

arrows(x0=1:10,y0=sm25EYxs[2,vord],y1=sm25EYxs[3,vord],angle=90,col=cols,code=3,length=0.05)
points(1:10,sm25EYxs[1,vord],pch=19,col=cols)


#--------------------------------
# plot seroprevalence curves means
#--------------------------------
xtics <- seq(0,60,by=10)
ytics <- seq(0,1,by=0.1)

# SEA
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=range(ytics),bty="n")
axis(1,at=xtics,las=1,cex.axis=1.5)
axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100), las=1,cex.axis=1.25)
mtext("Age-dependent seroprevalence, by community",side=3,line=2.5,adj=0.5,cex=1.1)
# mtext("Age, years",side=1,line=3)
mtext("Seroprevalence (%)",side=2,line=3,cex=1.1)
mtext(expression(paste(italic("S. mansoni")," SEA")),side=3,line=0.5,adj=0,cex=1.1,font=1)
mtext("B",side=3,line=2.5,adj=0,at=-15,cex=1.5,font=2)

j <- 1
for(i in vord){ 
  lines(psea_curves[[i]]$Age,psea_curves[[i]]$pY,col=cols[j])
  j <- j+1
}


# Sm25
xtics <- seq(0,60,by=10)
ytics <- seq(0,1,by=0.1)
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=range(ytics),bty="n")
axis(1,at=xtics,las=1,cex.axis=1.5)
axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100), las=1,cex.axis=1.25)
# mtext("Age-dependent seroprevalence, by community",side=3,line=0.5,adj=0,cex=1)
mtext("Age, years",side=1,line=3,cex=1.1)
mtext("Seroprevalence (%)",side=2,line=3,cex=1.1)
mtext(expression(paste(italic("S. mansoni")," Sm25")),side=3,line=0.5,adj=0,cex=1.1,font=1)

j <- 1
for(i in vord){ 
  lines(psm25_curves[[i]]$Age,psm25_curves[[i]]$pY,col=cols[j])
  j <- j+1
}

#--------------------------------
# plot community-level means
#--------------------------------
xtics <- 1:10

# SEA
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=range(ytics),bty="n")
axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100), las=1,cex.axis=1.25)
mtext(1:10,side=1,line=0,at=1:10,col=cols)
mtext("Seroprevalence, by community",side=3,line=2.5,adj=0.5,cex=1.1)
# mtext("Community",side=1,line=3)

arrows(x0=1:10,y0=sea_prev[vord,5],y1=sea_prev[vord,6],angle=90,col=cols,code=3,length=0.05)
points(1:10,sea_prev[vord,4],pch=19,col=cols)


# Sm25
plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=range(xtics),ylim=range(ytics),bty="n")
axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100), las=1,cex.axis=1.25)
mtext(1:10,side=1,line=0,at=1:10,col=cols)
# mtext("Seroprevalence, by community",side=3,line=0.5,adj=0)
mtext("Community",side=1,line=3,cex=1.1)

arrows(x0=1:10,y0=sm25_prev[vord,5],y1=sm25_prev[vord,6],angle=90,col=cols,code=3,length=0.05)
points(1:10,sm25_prev[vord,4],pch=19,col=cols)

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




