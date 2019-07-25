

#--------------------------------
# coastal-LF-boxplot.R
#
# plot LF Ab distributions by
# ICT status
#--------------------------------

#--------------------------------
# preamble
#--------------------------------
rm(list=ls())
library(scales)

#--------------------------------
# load saved results
#--------------------------------
load(file="~/dropbox/coastalkenya/results/raw/coastal-LF.RData")

#--------------------------------
# color pallete
#--------------------------------
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# brighter color blind palette:  https://personal.sron.nl/~pault/ 
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
cols=c(cblack,"gray50",cblue,cteal,cgreen,cchartr,corange,cred,cbPalette[8],cmagent)


#--------------------------------
# conduct mann-whitnet U-test
# for differences by ICT status
#--------------------------------

wb123u <- wilcox.test(log10(d$wb123),d$ict)
  wb123u
bm14u <- wilcox.test(log10(d$bm14),d$ict)
  bm14u
bm33u <- wilcox.test(log10(d$bm33),d$ict)
  bm33u

#--------------------------------
# make a boxplot of 
# Ab response stratified by
# ICT status
#--------------------------------

pdf("~/dropbox/coastalkenya/results/figs/coastal-LF-ICT-boxplot.pdf",height=5,width=10)
op <- par(mar=c(3,5,4,0)+0.1)
lo <- layout(mat=matrix(1:3,nrow=1,ncol=3))
ytics <- 0:5
xlims <- c(0,1)
xs <- c(0.25,0.75)

# Wb123
wb123n <- table(d$ict[!is.na(d$wb123)])
plot(1,1,ylim=range(ytics),xlim=xlims,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
bpt <- boxplot(log10(wb123)~ict,data=d,
               at=xs,
               lty=1, cex=0.4,staplewex=0.25,
               boxwex=0.25,
               col=alpha(c(cteal,corange),alpha=0.8),
               outcol=alpha(c(cteal,corange),alpha=0.4),pch=19,
               axes=FALSE,
               add=TRUE)
axis(2,at=ytics,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4),
  expression(10^5)), las=1,cex.axis=1.5)
mtext("Antibody response (MFI-bg)",side=2,line=3)
mtext(c("ICT neg","ICT pos"),side=1,line=0,at=xs,cex=1.25,col=c(cteal,corange))
mtext(paste("n =",format(wb123n,big.mark = ",")),side=1,line=1.5,at=xs,cex=0.9,col="gray20")
mtext("Wb123",side=3,line=0,adj=0,cex=1.25)


# Bm14
bm14n <- table(d$ict[!is.na(d$bm14)])
plot(1,1,ylim=range(ytics),xlim=xlims,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
bpt <- boxplot(log10(bm14)~ict,data=d,
               at=xs,
               lty=1,cex=0.4,staplewex=0.25,
               boxwex=0.25,
               col=alpha(c(cteal,corange),alpha=0.8),
               outcol=alpha(c(cteal,corange),alpha=0.4),pch=19,
               axes=FALSE,
               add=TRUE)
axis(2,at=ytics,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4),
  expression(10^5)), las=1,cex.axis=1.5)
# mtext("Antibody response (MFI-bg)",side=2,line=3)
mtext(c("ICT neg","ICT pos"),side=1,line=0,at=xs,cex=1.25,col=c(cteal,corange))
mtext(paste("n =",format(bm14n,big.mark = ",")),side=1,line=1.5,at=xs,cex=0.9,col="gray20")
mtext("Bm14",side=3,line=0,adj=0,cex=1.25)



# Bm33
bm33n <- table(d$ict[!is.na(d$bm33)])
plot(1,1,ylim=range(ytics),xlim=xlims,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
bpt <- boxplot(log10(bm33)~ict,data=d,
               at=xs,
               lty=1,cex=0.4,staplewex=0.25,
               boxwex=0.25,
               col=alpha(c(cteal,corange),alpha=0.8),
               outcol=alpha(c(cteal,corange),alpha=0.4),pch=19,
               axes=FALSE,
               add=TRUE)
axis(2,at=ytics,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4),
  expression(10^5)), las=1,cex.axis=1.5)
# mtext("Antibody response (MFI-bg)",side=2,line=3)
mtext(c("ICT neg","ICT pos"),side=1,line=0,at=xs,cex=1.25,col=c(cteal,corange))
mtext(paste("n =",format(bm33n,big.mark = ",")),side=1,line=1.5,at=xs,cex=0.9,col="gray20")
mtext("Bm33",side=3,line=0,adj=0,cex=1.25)

layout(mat=matrix(1))
par(op)
dev.off()






