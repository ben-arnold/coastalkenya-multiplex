

#--------------------------------
# coastal-Ab-distribution-figure.R
#
# summarize coastal kenya Ab
# antibody distributions
#--------------------------------


#--------------------------------
# preamble
#--------------------------------
rm(list=ls())
library(scales)

#--------------------------------
# load gaussian mixture model
# results
#--------------------------------
load(file="~/dropbox/coastalkenya/results/raw/coastal-mixtures.RData")

#--------------------------------
# load the formatted dataset
#--------------------------------
load("~/dropbox/coastalkenya/data/final/coastal_kenya.RData")
d <- coastal_kenya

#--------------------------------
# set negative and zero values to
# 1 before the log transform
#--------------------------------
table(d$wb123<=0)
d["wb123"][d["wb123"]<=0] <-1
# bm14 (all >0)
# bm33 (all >0)
# csp (all >0)
table(d$msp1pf<=0)
d["msp1pf"][d["msp1pf"]<=0] <-1
# msp1pm (all >0)
# strongy NIE (all >0)
table(d$ascaris<=0)
d["ascaris"][d["ascaris"]<=0] <- 1

# schisto SEA (all >0)
table(d$sm25<=0)
# for sm25 -- set all to missing (N=351 seems high)
d["sm25"][d["sm25"]<=0] <-NA


#--------------------------------
# kernel density plotting function
#--------------------------------

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
cols=c(cblack,"gray50",cblue,cteal,cgreen,cchartr,corange,cred,cmagent)


denplot <- function(x,col,main,header=NULL,xlab=NULL,cutoff=NULL,cutoffmix=NULL,cutlab=NULL) {
  dest <- density(x[!is.na(x)])
  xtics <- 0:5
  ytics <- seq(0,1.2,by=0.1)
  plot(dest,type="n",
       main="",
       ylim=range(ytics),
       xlim=range(xtics),xaxt="n",xlab="",
       las=1,bty="n"
  )
  axis(1,at=0:5,labels=c(
    expression(10^0),
    expression(10^1),
    expression(10^2),
    expression(10^3),
    expression(10^4),
    expression(10^5)
  ), las=1,cex.axis=1.25
  )
  mtext(header,side=3,line=3,adj=0,at=-0.5,cex=1.25,font=2)
  mtext(main,side=3,line=1,adj=0,cex=1.1)
  mtext(xlab,side=1,line=2.5,cex=1.1)
  polygon(dest,col=alpha(col,alpha=0.75),border=col)
  
  if(!is.null(cutoff)){
    segments(x0=cutoff,y0=min(ytics),y1=max(ytics),col="gray20")
  }
  if(!is.null(cutoffmix)){
    segments(x0=cutoffmix,y0=min(ytics),y1=max(ytics),col="gray20",lty=2)
  }
  if(!is.null(cutoff) | !is.null(cutoffmix)){
    text(x=max(cutoff,cutoffmix),y=1.1,cutlab,col="gray20",pos=4,cex=1.2)
  }

}


#--------------------------------
# composite plot with cutoffs
#--------------------------------

pdf("~/dropbox/coastalkenya/results/figs/coastal-Ab-distributions-cutoffs.pdf",height=15,width=12)
op <- par(mar=c(4,4,5,1)+0.1)
lo <- layout(mat=matrix(1:15,nrow=5,ncol=3,byrow=T))
cols <- c(cblue,cteal,cgreen,corange,cmagent)

### vaccine preventable

# measles
denplot(log10(d$measles),col=cols[1],main="Measles MV-N",header="Vaccine Preventable",
        cutoff=log10(178),cutlab="Seroprotected")

# diphtheria
denplot(log10(d$diphtheria),col=cols[1],main="Diphtheria toxoid",
        cutoff=log10(4393),cutlab="Seroprotected")

# tetanus
denplot(log10(d$tetanus),col=cols[1],main="Tetanus toxoid",
        cutoff=log10(118),cutlab="Seroprotected")

### Malaria
# CSP
denplot(log10(d$csp),col=cols[2],main=expression(paste(italic('Plasmodium falciparum')," CSP")),header="Malaria",
        cutoffmix=log10(766),cutlab="Seropositive")

# MSP-1 falciparum
denplot(log10(d$msp1pf),col=cols[2],main=expression(paste(italic('Plasmodium falciparum '),MSP-1[19])),
        cutoff=log10(256),cutoffmix=log10(pfmixcut),cutlab="Seropositive")

# MSP-1 malariae
denplot(log10(d$msp1pm),col=cols[2],main=expression(paste(italic('Plasmodium malariae '),MSP-1[19])),
        cutoff=log10(403),cutoffmix=log10(pmmixcut),cutlab="Seropositive")

### LF
# Wb123
denplot(log10(d$wb123),col=cols[3],main="LF Wb123",header="Lymphatic filariasis",
        cutoff=log10(342),cutoffmix=log10(wb123mixcut),cutlab="Seropositive")

# Bm14
denplot(log10(d$bm14),col=cols[3],main="LF Bm14",
        cutoff=log10(444),cutoffmix=log10(bm14mixcut),cutlab="Seropositive")

# Bm33
denplot(log10(d$bm33),col=cols[3],main="LF Bm33",xlab="Antibody response (MFI-bg)",
        cutoff=log10(369),cutoffmix=log10(bm33mixcut),cutlab="Seropositive")


### Schistosomiasis
# Schisto SEA
denplot(log10(d$sea),col=cols[4],main=expression(paste(italic('Schistosoma mansoni')," SEA")),header="Schistosomiasis",
        cutoff=log10(965),cutoffmix=log10(seamixcut),cutlab="Seropositive")

# Schisto Sm25
denplot(log10(d$sm25),col=cols[4],main=expression(paste(italic('Schistosoma mansoni')," Sm25")),
        cutoff=log10(38),cutoffmix=log10(sm25mixcut),cutlab="Seropositive")

# empty
plot(1,1,type="n",yaxt="n",xaxt="n",ylab="",xlab="",bty="n")

### Other helminths
# Strongy NIE
denplot(log10(d$nie),col=cols[5],main=expression(paste(italic('Strongyloides stercoralis')," NIE")),
        header="Other helminths",xlab="Antibody response (MFI-bg)",
        cutoff=log10(628),cutoffmix=log10(niemixcut),cutlab="Seropositive")

denplot(log10(d$ascaris),col=cols[5],main=expression(paste(italic('Ascaris')," spp. AsHb")),
        xlab="Antibody response (MFI-bg)",
        cutoff=log10(386),cutoffmix=log10(ascarismixcut),cutlab="Seropositive")

# empty
plot(1,1,type="n",yaxt="n",xaxt="n",ylab="",xlab="",bty="n",ylim=c(0,1),xlim=c(0,1))
legend(0,0.75,legend=c("Cutoff from standard curve\n(vaccine prev), or unexposed\nUS adults (others)","Cutoff from Gaussian mixture\nmodel"),
       lty=c(1,2),col="gray20",bty="n",cex=1.6)


par(op)
layout(mat=matrix(1))
dev.off()














