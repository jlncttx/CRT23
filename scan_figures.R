library(RColorBrewer)
library(sf)
library(graticule)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(cowplot)

#---------------------------------------------------------------------------------
# Color palettes and levels definitions
#---------------------------------------------------------------------------------
# Unique color scale used in CR23 - 13 colors for 14 levels (incl min/max)
COLSM=c("lightblue",brewer.pal(9,"YlOrRd"),"magenta4","magenta","plum1")
#COLSM=c("lightblue",rev(heat.colors(length(BRKSM)-6)),"red4","magenta4","magenta","pink")

# Levels for 1d-heatmap (in percentile levels) (the 2nd one is zoomed at the end)
BRKS1=c(0,50,80,90,95,98,99,99.5,99.8,99.9,99.95,99.99,99.999,100)
BRKS2=c(BRKS1[-c(3:6,length(BRKS1))],99.9996,99.9997,99.9998,99.9999,100)

# Levels for 2d-maps
BRKSM=c(0,50,70,80,90,95,97,98,99,99.5,99.7,99.8,99.9,100)
#BRKSM=100*(1-1/c(1,2,5,8,10,20,50,80,100,200,500,800,1000))#c(0,50,80,90,95,98,99,99.5,99.8,99.9,100)

# Levels for chronology
BRKS0=c(0,75,90,95,96,97,98,99,99.5,99.8,99.9,99.99,99.999,100)

# Other color scales
COLS0=c("lightblue",rev(heat.colors(length(BRKS0)-7)),"red3","red4","magenta4","magenta","pink")
COLS1=c("lightblue",rev(heat.colors(length(BRKS1)-6)),"red3","red4","magenta4","violet")
COLS2=c(COLS1[-c(3:5,length(COLS1))],"magenta3","magenta","violet","white")

#---------------------------------------------------------------------------------
# Proj4strings for crs (map projections)
#---------------------------------------------------------------------------------
crs.robinson=list(crs="+proj=robin +lon_0=0 +x_0=0 +y_0=0 +units=m",xlim=NULL,ylim=NULL)

crs.lonlat=function(xlim=NULL,ylim=NULL)
  list(crs="+proj=longlat",xlim=xlim,ylim=ylim)

crs.lambert=function(lon0,lat0,xkm,ykm)
  list(crs=paste0("+proj=laea +lon_0=",lon0," +lat_0=",lat0," +x_0=0 +y_0=0 +units=m"),
       xlim=c(-1,1)*xkm*1e3,ylim=c(-1,1)*ykm*1e3)

# France
crs.France.lonlat=crs.lonlat(xlim=c(-4.8,10.2),ylim=c(41,51.6))
crs.France.lambert=crs.lambert(2.7,46.3,550,550)

#---------------------------------------------------------------------------------
# A few things for the legends
#---------------------------------------------------------------------------------
# Colors / names for the 3 methods
method.color=function(m){
  out=1
  if (m %in% c("calend","calend2")) out="darkmagenta"
  if (m %in% c("annmax","annmax1","annmax2","annmax3")) out=2
  if (m %in% c("annmin","annmin1","annmin2","annmin3")) out=4
  if (m %in% c("locmax","locmin")) out=3
  out
}

method.color2=function(m){
  out=1
  if (m %in% c("calend","calend2")) out="lavenderblush"
  if (m %in% c("annmax","annmax1","annmax2","annmax3")) out="antiquewhite"
  if (m %in% c("annmin","annmin1","annmin2","annmin3")) out="aliceblue"
  if (m %in% c("locmax","locmin")) out="palegreen"
  out
}

method.name=function(m){
  out="METHOD"
  if (m=="calend") out="Calendar"
  if (m=="calend2") out="Calendar SN"
  if (m=="annmax") out="Annual-maxima"
  if (m=="annmin") out="Annual-minima"
  if (m=="annmax1") out="Annual-maxima GEVcon"
  if (m=="annmin1") out="Annual-minima GEVcon"
  if (m=="annmax2") out="Annual-maxima Gumbel"
  if (m=="annmin2") out="Annual-minima Gumbel"
  if (m=="annmax3") out="Annual-maxima Normal"
  if (m=="annmin3") out="Annual-minima Normal"
  if (m=="locmax") out="Local-maxima"
  if (m=="locmin") out="Local-minima"
  out
}

method.sample=function(m,w,id1,id2){
  out="SAMPLE"
  if (m %in% c("calend","calend2")) out=paste("on",idays2leg(id1,id2))
  if (m %in% c("annmax","annmax1","annmax2","annmax3")) out=paste0(w,"-day maxima")
  if (m %in% c("annmin","annmin1","annmin2","annmin3")) out=paste0(w,"-day minima")
  if (m %in% c("locmax","locmin")) out=paste("near",idays2leg(id1,id2))
  out
}

distrib.name=function(d){
  out="DISTRIB"
  if (d=="gauss") out="Normal"
  if (d=="snorm") out="Skew-N"
  if (d %in% c("gev","gevmin")) out="GEV"
  if (d %in% c("gevmle","gevminmle")) out="GEV-const_ksi"
  if (d %in% c("gumbel","gumbelmin")) out="Gumbel"
  out
}

var.name=function(v){
  if (v %in% c("T","TM","TX","TN")) out="Temperature"
  paste0(out," (",var.unit(v),")")
}

var.unit=function(v){
  if (v %in% c("T","TM","TX","TN")) out="K"
  out
}

var.leg=function(v) paste(v,"in",var.unit(v))

# Months names
MONTHS=c("January","February","March","April","May","June","July",
         "August","September","October","November","December")
MTHS=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
MTHS1=c("J","F","M","A","M","J","J","A","S","O","N","D")

# Calendar axis
my.calendar.axis=function(side=1,mlab=MTHS,dlab=c(1,15),lab=T,lty=3,lwd=0.5,tcl=0.5){
  d0=makedates(y=1:3,avg="day",caltype="noleap")
  if (side %in% c(1,3)) abline(v=which(d0$d==1)-0.5,lty=lty,col='gray')
  if (side %in% c(2,4)) abline(h=which(d0$d==1)-0.5,lty=lty,col='gray')
  dd=which(d0$d %in% seq(5,25,5) & ! d0$d %in% dlab)
  axis(side,at=dd,lab=F,tcl=tcl/2)
  dd=which(d0$d %in% dlab)
  if (lab)  axis(side,at=dd,lab=d0$d[dd],cex.axis=0.66,mgp=c(0,-0.2,0),tcl=tcl)
  if (!lab) axis(side,at=dd,lab=F,tcl=tcl)
  if (lab) {
    dd=which(d0$d==15)
    for (id in dd) axis(1,at=id,lab=mlab[d0$m[id]],tcl=0,mgp=c(0,0.6,0))
  }
}

# Title right axis
my.right.axis.title=function(title="title"){
  lims=par("usr"); par(xpd = TRUE)
  text(lims[2]+0.08*(lims[2]-lims[1]),mean(lims[3:4]),title,srt = 270)
  par(xpd=FALSE)
}

# Functions that convert percentile levels into probabilities (with scientific writing)
perc2prob=function(lev){
  lev=1-lev/100
  out=c(); for (l in lev){
    if (l>=0.00099) out=c(out,signif(l,3))
    else out=c(out,format(l,sci=T))
  }; return(out)
}

prob2leg=function(probs){
  out=c(); for (p in probs){
    dum=round(p,digits=2); if (dum==0) dum=signif(p,digits=1)
    out=c(out,dum)
  }; return(out)
}

# Time window in titles/legends
idays2leg=function(id1,id2,short=F){
  dd=makedates(y=1,avg="day",caltype="noleap")[c(id1,id2),]
  mm=MONTHS ; if (short) mm=MTHS
  out=paste(mm[dd$m[1]],dd$d[1])
  if (id2!=id1){
    if (dd$m[1]==dd$m[2]) out=paste(out,"-",dd$d[2])
    else out=paste(out,"-",mm[dd$m[2]],dd$d[2])
  }; return(out)
}

# Dates in titles/legends
dates2leg=function(date1,date2,short=F){
  dd=yyyymmdd2mdy(c(date1,date2))
  mm=MONTHS ; if (short) mm=MTHS
  out=paste(mm[dd$m[1]],dd$d[1])
  if (date2!=date1){
    if (dd$m[1]==dd$m[2]) out=paste(out,"-",dd$d[2])
    else out=paste(out,"-",mm[dd$m[2]],dd$d[2])
  }; return(out)
}

#---------------------------------------------------------------------------------
# Figure to illustrate the detrending
# Warning: s must be a calendar scan!
#---------------------------------------------------------------------------------
myfig.detrend=function(s=NULL,var="TM",loc="France",set="mf",
                       y=2022,f.df=12,h.df=6,
                       title=c("","","","")){
  # Data
  if (is.null(s)) s=get.scan.1d(var,"calend",loc,set)
  var=s$var; meth=s$method; loc=s$location; set=s$set
  # Periods data vs. cmip
  iy=which(CMIP.PERIOD %in% s$period)
  # Multi-model mean
  mmm=s$cmip.mmm
  # Forced response + centered over data period
  dum=forced_response_by_gam_fit(mmm,Enat,10)
  nat=dum$nat; ant=dum$ant; fr=dum$all
  cfr=fr[iy]-mean(fr[iy])
  # Yearly time series from scan (including incomplete years)
  yts=apply(s$data[,1,],2,mean,na.rm=T)
  ytr=apply(s$trend[,1,],2,mean,na.rm=T)
  # Calendar scaling
  f_cal=t(apply(s$data[,1,],1,function(x) lm(x~cfr)$coeff))
  sf_cal=cbind(mysmoothy(f_cal[,1],df=f.df,per=T),mysmoothy(f_cal[,2],df=h.df,per=T))
  cyc1=sf_cal[,1]+sf_cal[,2]*cfr[which(s$period==s$yref)]
  cyc2=sf_cal[,1]+sf_cal[,2]*cfr[which(s$period==y)]
  sf_tyr=mean(sf_cal[,2])
  # Annmax/annmin scaling
  sx=get.scan.1d(var,"annmax",loc,set)
  sn=get.scan.1d(var,"annmin",loc,set)
  f_tsx=t(apply(sx$data[1,,],1,function(x) lm(x~cfr)$coef))[,2]
  f_tsn=t(apply(sn$data[1,,],1,function(x) lm(x~cfr)$coef))[,2]
  sf_tsx=mysmoothy(f_tsx,sdf=2)
  sf_tsn=mysmoothy(f_tsn,sdf=2)
  # Plot -------------------------------------------------------
  par(mar=c(3,3,1.5,3),mgp=c(1.6,0.2,0),tcl=0.5,las=1)
  layout(cbind(c(1,1,2,2,7),c(3,4,5,6,8)),height=c(1,0.3,0.3,0.5,1.2))
  #  - annual cycle (1)
  par(mar=c(0,3,1.5,3)); cols=c(1,"darkmagenta",4,2)
  ylim=range(c(f_cal[,1],cyc1,cyc2)); ylim=ylim+c(-0.05,0.05)*(ylim[2]-ylim[1])
  matplot(cbind(f_cal[,1],sf_cal[,1],cyc1,cyc2),type="l",lwd=c(1,2,1,1),col=cols,lty=c(1,1,2,2),
          axes=F,frame.plot=T,ylim=ylim,main="",xlab="",ylab=var.leg(var))
  for (i in c(1,3)) my.calendar.axis(i,lab=F)
  for (i in c(2,4)) axis(i,at=seq(-100,100,5),lab=(i==2))
  for (i in c(2,4)) axis(i,at=seq(-100,100,1),lab=F,tcl=0.25)
  legend("topleft",bty="n",inset=0.02,
         leg=c(expression(paste(T[cyc]," raw")),expression(paste(T[cyc]," smoothed")),
               paste(s$yref,"normals"),paste(y,"normals")),
         lwd=c(1,2,1,1),col=cols,lty=c(1,1,2,2),text.col=cols)
  title(paste(title[1],var,loc,"annual cycle"),adj=0)
  #  - delta annual cycle (2)
  par(mar=c(3,3,0,3))
  ylim=range(cyc2-cyc1); ylim=ylim+c(-0.1,0.1)*(ylim[2]-ylim[1])
  plot(cyc2-cyc1,type="l",lwd=2,col=cols[2],axes=F,frame.plot=T,
       ylim=ylim,main="",xlab="Calendar days",ylab=var.unit(var))
  abline(h=mean(cyc2-cyc1),col=cols[2],lty=2)
  for (i in c(1,3)) my.calendar.axis(i,MTHS1,lab=(i==1))
  for (i in c(2,4)) axis(i,lab=(i==2))
  legend("topleft",bty="n",inset=0.02,leg=paste(y,"vs.",s$yref),
         lwd=2,col=cols[2],text.col=cols[2])
  #  - yearly trend (3)
  par(mar=c(0,3,1.5,3))
  xlim=range(s$period)+c(-10,10)
  ylim=range(yts,na.rm=T); ylim=ylim+c(-0.05,0.05)*(ylim[2]-ylim[1])
  matplot(s$period,cbind(yts,ytr),type="l",lty=1,lwd=c(1,2),col=c(1,4),
          axes=F,frame.plot=T,xlim=xlim,ylim=ylim,
          main="",xlab="",ylab=var.leg(var))
  if (any(is.na(s$data[,1,length(s$period)-1])))
    text(s$period[length(yts)-1],yts[length(yts)-1],"!",pos=4,font=2,col=2)
  points(y,yts[which(s$period==y)],pch=8,cex=1.5)
  abline(h=mean(ytr),lty=3)
  abline(v=range(s$years),lty=3)
  for (i in c(1,3)) axis(i,at=seq(1800,2100,10),lab=F)
  for (i in c(1,3)) axis(i,at=seq(1800,2100,2),lab=F,tcl=0.25)
  axis(2,seq(-100,100,1)); axis(2,at=seq(-100,100,0.2),lab=F,tcl=0.25)
  axis(4,at=mean(ytr)+seq(-100,100,1),lab=seq(-100,100,1))
  axis(4,at=mean(ytr)+seq(-100,100,0.2),lab=F,tcl=0.25)
  my.right.axis.title("Centered")
  legend("top",bty="n",inset=0.02,leg=c("T",expression(paste(beta[T]," . F"))),
         lty=1,lwd=c(1,2),col=c(1,4),text.col=c(1,4))
  title(paste(title[2],var,loc,"yearly trend"),adj=0)
  #  - mmm and F (4)
  par(mar=c(0,3,0,3))
  ylim=range(mmm[iy]); ylim=ylim+c(-1,1)*0.5*(ylim[2]-ylim[1])
  matplot(CMIP.PERIOD,cbind(mmm,fr),type="l",lty=1,lwd=c(1,2),col=c(1,4),
          axes=F,frame.plot=T,xlim=xlim,ylim=ylim,
          main="",xlab="",ylab=var.unit(var))
  abline(h=mean(fr[iy]),lty=3)
  abline(v=range(s$years),lty=3)
  for (i in c(1,3)) axis(i,at=seq(1800,2100,10),lab=F)
  for (i in c(1,3)) axis(i,at=seq(1800,2100,2),lab=F,tcl=0.25)
  axis(2,at=seq(-100,100,1))
  axis(4,at=mean(fr[iy])+seq(-100,100,1),lab=seq(-100,100,1))
  axis(4,at=mean(fr[iy])+seq(-100,100,0.5),lab=F,tcl=0.25)
  my.right.axis.title("Centered")
  legend("top",bty="n",inset=0.02,leg=c("MMM","F"),lty=1,lwd=c(1,2),col=c(1,4),text.col=c(1,4),ncol=2)
  #  - nat (5)
  par(mar=c(0,3,0,3))
  plot(CMIP.PERIOD,nat,type="l",lwd=2,col=3,axes=F,frame.plot=T,
       xlim=xlim,ylim=ylim+mean(nat[iy])-mean(mmm[iy]),
       main="",xlab="",ylab=var.unit(var))
  abline(v=range(s$years),lty=3)
  for (i in c(1,3)) axis(i,at=seq(1800,2100,10),lab=F)
  for (i in c(1,3)) axis(i,at=seq(1800,2100,2),lab=F,tcl=0.25)
  for (i in c(2,4)) axis(i,at=seq(-100,100,1),lab=(i==2))
  legend("top",bty="n",inset=0.02,leg="NAT",lty=1,lwd=2,col=3,text.col=3)
  #  - ant (6)
  par(mar=c(3,3,0,3))
  plot(CMIP.PERIOD,ant,type="l",lwd=2,col=2,axes=F,frame.plot=T,
       xlim=xlim,ylim=ylim+mean(ant[iy])-mean(mmm[iy]),
       main="",xlab="",ylab=var.unit(var))
  abline(v=range(s$years),lty=3)
  for (i in c(1,3)) axis(i,at=seq(1800,2100,10),lab=(i==1))
  for (i in c(1,3)) axis(i,at=seq(1800,2100,2),lab=F,tcl=0.25)
  for (i in c(2,4)) axis(i,at=seq(-100,100,1),lab=(i==2))
  legend("top",bty="n",inset=0.02,leg="ANT",lty=1,lwd=2,col=2,text.col=2)
  title(xlab="Years",mgp=c(1.2,0.1,0))
  #  - calendar scaling (7)
  par(mar=c(3,3,1.5,3))
  ylim=range(f_cal[,2]); ylim=ylim+c(-0.1,0.1)*(ylim[2]-ylim[1])
  matplot(cbind(f_cal[,2],sf_cal[,2]),type="l",lty=1,lwd=c(1,2),col=c(1,"darkmagenta"),
          ylim=ylim,axes=F,frame.plot=T,main="",xlab="Calendar days",ylab="Scaling factor")
  abline(h=1,lty=3); abline(h=sf_tyr,lty=2,col="darkmagenta")
  for (i in c(1,3)) my.calendar.axis(i,MTHS1,lab=(i==1))
  axis(2,seq(-10,10,1)); axis(2,seq(-10,10,0.2),lab=F,tcl=0.25)
  axis(4,at=seq(0,10,0.4)*sf_tyr,label=round(seq(0,10,0.4),1))
  my.right.axis.title("Relative to yearly scaling")
  legend("topleft",bty="n",inset=0.02,
         leg=c(expression(paste(beta[T]," raw")),expression(paste(beta[T]," smoothed"))),
         lty=1,lwd=c(1,2),col=c(1,"darkmagenta"),text.col=c(1,"darkmagenta"))
  title(paste(title[3],var,loc,"calendar trends"),adj=0)
  #  - annmax/annmin scaling (8)
  ylim=range(c(sf_tsx,sf_tsn)); ylim=ylim+c(-0.1,0.1)*(ylim[2]-ylim[1])
  matplot(cbind(sf_tsx,sf_tsn),type="l",lty=1,lwd=2,col=c(2,4),
          ylim=ylim,axes=F,frame.plot=T,main="",xlab="n in days",ylab="Scaling factor")
  matplot(cbind(sf_tsx,sf_tsn),pch=16,col=c(2,4),add=T)
  abline(h=1,lty=3); abline(h=sf_tyr,lty=2,col="darkmagenta")
  axis(1,at=1:length(s$windows),lab=s$windows,las=2)
  axis(3,at=1:length(s$windows),lab=F)
  axis(2,seq(-10,10,0.2)); axis(2,seq(-10,10,0.04),lab=F,tcl=0.25)
  axis(4,at=seq(0,10,0.2)*sf_tyr,label=round(seq(0,10,0.2),1))
  my.right.axis.title("Relative to yearly scaling")
  legend("topright",bty="n",inset=0.02,
         leg=c(expression(beta[Tx]),expression(beta[Tn])),
         lty=1,lwd=2,pch=16,col=c(2,4),text.col=c(2,4))
  title(paste(title[4],var,loc,"annual max/min trends"),adj=0)
}

#---------------------------------------------------------------------------------
# Figure to illustrate the methodology
#---------------------------------------------------------------------------------
myfig.meth=function(s=NULL,var="TM",meth="calend",loc="France",set="mf",
                    y=2022,id=168,iw=5,stat="p1",reverse.prob=F,
                    title=c("","")){
  # Data
  if (is.null(s)) s=get.scan.1d(var,meth,loc,set)
  var=s$var; meth=s$method; loc=s$location; set=s$set
  sy=s$scan[[paste0("Y",y)]]
  # Indices of days (id1,id2) corresponding to (id,iw)
  id1=mymod(id-floor((s$windows[iw]-1)/2),365)
  id2=mymod(id+floor(s$windows[iw]/2),365)
  # Value of the year of interest
  iy=which(s$period==y)
  dumy=sy$value[id,iw]
  # Year of reference
  iy0=which(s$period==s$yref)
  # Get data, fit, pdfs and percentiles
  if (s$params$annual.maxima) dum=s$data[1,iw,] else dum=s$data[id,iw,]
  dumtr=s$trend[id,iw,]
  dumtr1=dumtr[iy]
  dumtr0=dumtr[iy0]
  dum1=dum-dumtr+dumtr[iy]
  dum0=dum-dumtr+dumtr[iy0]
  pdf1=mypdffun(s$params$distrib,sy$fit1[id,iw,])
  pdf0=mypdffun(s$params$distrib,sy$fit0[id,iw,])
  p1=myperclevel(dumy,s$params$distrib,sy$fit1[id,iw,],exceedance=s$params$exceedance)
  p0=myperclevel(dumy,s$params$distrib,sy$fit0[id,iw,],exceedance=s$params$exceedance)
  if (reverse.prob) {p1=1-p1 ; p0=1-p0}
  pr=(1-p1)/(1-p0)
  # Verif
  pp1=sy$p1[id,iw] ; pp0=sy$p0[id,iw]
  if (reverse.prob) {pp1=1-pp1 ; pp0=1-pp0}
  ppr=(1-pp1)/(1-pp0)
  # Exchange 0 and 1 if stat=p0
  if (stat=="p0"){
    x=dumtr1; dumtr1=dumtr0; dumtr0=x
    x=dum1; dum1=dum0; dum0=x
    x=pdf1; pdf1=pdf0; pdf0=x
  }
  # Panel plot
  layout(matrix(1:3,1,3),widths=c(1,0.1,1))
  coly="darkmagenta"
  colm=method.color(meth)
  # Plot 1: time series
  par(mar=c(3,3,1.5,0.1),mgp=c(1.6,0.2,0),tcl=0.5,las=1)
  ylim=range(c(dum,dum1),na.rm=T); ylim=ylim+c(-1,1)*0.1*(ylim[2]-ylim[1])
  matplot(s$period,cbind(dum,dum1,dumtr,rep(dumtr1,length(s$period))),type="l",
          lty=c(1,2,1,2),lwd=c(1.5,1,3,1),col=c(1,1,colm,colm),
          axes=F,frame.plot=T,ylim=ylim,main="",xlab="Years",ylab=var.name(var))
  # Value of interest
  abline(v=y,col=coly); abline(h=dumy,col=coly)
  points(y,dumy,col=coly,pch=16,cex=1.5)
#  text(y,ylim[1],"Event",srt=90,adj=c(0,-0.5),font=2,col=coly)
  text(y,ylim[1],paste0(" Event = ",idays2leg(id1,id2),", ",y," "),adj=c((y>mean(s$period)),0),font=2,col=coly)
  # Axes, title and legends
  for (i in c(1,3)) axis(i,at=seq(1800,2100,10),lab=(i==1))
  for (i in c(1,3)) axis(i,at=seq(1800,2100,2),lab=F,tcl=0.25)
  for (i in c(2,4)) axis(i,at=seq(-100,100,1),lab=(i==2))
  for (i in c(2,4)) axis(i,at=seq(-100,100,0.2),lab=F,tcl=0.25)
  legend("topleft",bty="n",inset=0.01,leg=c("T raw","T detrended","Long-term trend",paste(y,"climate")),
         lty=c(1,2,1,2),lwd=c(1.5,1,3,1),col=c(1,1,colm,colm),ncol=2)
  title(paste(title[1],var,loc,method.sample(meth,s$windows[iw],id1,id2)),adj=0)
  title(method.name(meth),adj=1,col.main=colm)
  # Plot 2: vertical barcode
  par(mar=c(3,0.1,1.5,0.5))
  plot0(ylim=ylim)
  abline(h=dum1)
  abline(h=dumy,col=coly)
  # Plot 3: histogram + fit
  par(mar=c(3,3,1.5,0.5))
  xlim=ylim; ylim=c(-0.3,1.5)*pdf1(mean(dum1,na.rm=T))
  # Densities
  plot(pdf1,xlim[1],xlim[2],ylim=ylim,lwd=3,col=colm,axes=F,xlab=var.name(var),ylab="Probability",frame.plot=T,type="n")
  hist(dum1,prob=T,add=T,col=method.color2(meth),breaks = 10); abline(h=0,lwd=1.5)
  plot(pdf1,xlim[1],xlim[2],ylim=ylim,lwd=3,col=colm,add=T)
  plot(pdf0,xlim[1],xlim[2],ylim=ylim,lwd=2,col="gray",add=T)
  abline(v=dumy,lwd=2,col=coly); text(dumy,ylim[2],"Event",srt=90,adj=c(1,-0.5),font=2,col=coly)
  # Horizontal barcodes
  for (d in dum1) segments(d,0.9*ylim[1],d,0.2*ylim[1],lwd=0.5)
  # KS test + trend value
  leg.ks=paste("KS-test p-value =",round(mykstest(dum1,s$params$distrib,sy$fit1[id,iw,]),2))
  leg.tr=paste("Estimated long-term trend =",round(dumtr1-dumtr0,2),"K")
  legend("topleft",bty="n",inset=0.01,c(leg.ks,leg.tr))
  # Axes, title and legends
  for (i in c(1,3)) axis(i,at=seq(-100,100,1),lab=(i==1))
  for (i in c(1,3)) axis(i,at=seq(-100,100,0.2),lab=F,tcl=0.25)
  for (i in c(2,4)) axis(i,at=c(-999,seq(0,1,0.1)),lab=(i==2))
  title(paste(title[2],distrib.name(s$params$distrib),"fit"),adj=0)
  title(paste("p1 =",prob2leg(1-p1),"; p0 =",prob2leg(1-p0),"; pr=",round(pr)),adj=1,col.main=colm)
#  legend("topright",paste("p1 =",prob2leg(1-pp1),"; p0 =",prob2leg(1-pp0),"; pr=",round(ppr,2)),bty="n")
}

#---------------------------------------------------------------------------------
# Figure that shows the event in the calendar context
# Warning: s must be a calendar scan!
#---------------------------------------------------------------------------------
myfig.eve=function(s=NULL,var="TM",loc="France",set="mf",
                   y=2022,id=168,iw=5,
                   title=""){
  # Get data of years y-1, y and y+1
  if (is.null(s)) s=get.scan.1d(var,"calend",loc,set)
  var=s$var; meth=s$method; loc=s$location; set=s$set
  iy=which(s$period==y)
  dum=c(); for (i in c(iy-1,iy,iy+1)) dum=c(dum,s$data[,1,i])
  # Get normals of corresp years
#  dumtr1=c(); for (i in c(iy-1,iy,iy+1)) dumtr1=c(dumtr1,s$trend[,1,i])
  mu=rep(s$scan[[paste0("Y",y)]]$fit1[,1,1],3)
  sd=rep(s$scan[[paste0("Y",y)]]$fit1[,1,2],3)
  # Get normals of reference
  mu0=rep(s$trend[,1,which(s$period==s$yref)],3)
  # Indices of days (id1,id2) corresponding to (id,iw)
  id1=id-floor((s$windows[iw]-1)/2)+365
  id2=id+floor(s$windows[iw]/2)+365
  # Days to plot
  idays=min(id1-30,366):max(id2+30,730)
  # Data to plot
  dum=dum[idays]; mu=mu[idays]; sd=sd[idays]; mu0=mu0[idays]
  # Plot
  par(mar=c(3,3,1.5,0.5),mgp=c(1.5,0.1,0),las=1,tcl=0.5)
  ylim=range(c(dum,mu+2*sd,mu-2*sd),na.rm=T)
  ylim=ylim+c(-0.1,0.1)*(ylim[2]-ylim[1])
  plot0(xlim=range(idays),ylim=ylim,frame.plot=T,xlab="Calendar days",ylab=var.name(var))
  # Shading  
  polygon(c(idays,rev(idays)),c(mu+sd,rev(mu+2*sd)),col="gray90",border=NA)
  polygon(c(idays,rev(idays)),c(mu-sd,rev(mu+sd)),col="gray80",border=NA)
  polygon(c(idays,rev(idays)),c(mu-sd,rev(mu-2*sd)),col="gray90",border=NA)
  polygon(c(id1-0.5,id2+0.5,id2+0.5,id1-0.5),rep(c(-999,999),each=2),col=rgb(0.5,0.75,1,0.5),border=NA)
  # Axis here because of vertical segments
  for (i in c(1,3)) my.calendar.axis(i,lab=(i==1))
  # Lines
  matplot(idays,cbind(mu,mu0,dum),add=T,type="l",lty=c(1,2,1),lwd=c(3,1.5,1),col=c(4,4,1))
  # Axes, legend and title
  for (i in c(2,4)) axis(i,at=seq(-100,100,5),lab=(i==2))
  for (i in c(2,4)) axis(i,at=seq(-100,100,1),lab=F,tcl=0.25)
  legend("topleft",bty="n",inset=0.01,leg=c(paste(y,"values"),paste(y,"normals"),"(+/- 1 and 2 sd)",paste(s$yref,"normals")),
         lty=c(1,1,0,2),lwd=c(1,3,0,1.5),col=c(1,4,0,4))
  title(paste(title,var,loc,y),adj=0)
  title(idays2leg(mymod(id1,365),mymod(id2,365)),adj=1,col.main=rgb(0.5,0.75,1))
}

#---------------------------------------------------------------------------------
# Figure local p1 (1d)
#---------------------------------------------------------------------------------
myfig.1d=function(s=NULL,var="TM",meth="calend",loc="France",set="mf",
                  y=2003,stat="p1",reverse.prob=F,
                  lev=BRKS1,col=COLSM,title="",image=T,colorbar=T,
                  show.events=1,min.value=NULL,min.win=NULL,max.win=NULL){
  # Get data
  if (is.null(s)) s=get.scan.1d(var,meth,loc,set)
  var=s$var; meth=s$method; loc=s$location; set=s$set
  sy=s$scan[[paste0("Y",y)]]
  mat=as.matrix(sy[[paste(stat)]])
  if (stat %in% c("p1","p0")){
    if (reverse.prob) mat=1-mat
    mat=100*mat
  }
  # Panel plot (image + colorbar)
  # Warning: colorbar first to allow for click on image in shiny app
  if (image & colorbar) layout(matrix(2:1,1,2),width=c(1,0.2))
  # Plot 1: colorbar
  if (colorbar){
    par(mar=c(3.5,0,1.5,1))
    levleg=lev[2:(length(lev)-1)]; if (stat %in% c("p1","p0")) levleg=perc2prob(levleg)
    mycolorbar(col,levleg,horiz=F,width=0.2,cex=1)
  }
  # Plot 2: image
  if (image){
    par(mar=c(3,3,1.5,0.5),mgp=c(1.5,0.1,0),las=1,tcl=0.5)
    image(1:nrow(mat),1:ncol(mat),mat,axes=F,frame.plot=T,main="",breaks=lev,col=col,
          xlab="Central day",ylab="Duration in days")
    # Axes
    for (i in c(1,3)) my.calendar.axis(i,lab=(i==1))
    for (i in c(2,4)) axis(i,c(-999,999),lab=F)
    abline(h=which(s$windows %in% c(5,10,30,90)),lty=3,col="gray")
    axis(2,at=1:length(s$windows),lab=s$windows)
    axis(4,at=1:length(s$windows),lab=F)
    # Titles
    title(paste(title,stat,var,loc,y),adj=0)
    title(method.name(meth),adj=1,col.main=method.color(meth))
    # Adding the global max
    if (show.events==1){
      ijmax=which(mat==max(mat,na.rm=T),arr.ind=T)
      points(ijmax,pch="X",font=2,cex=1.5)
    }
    # Adding max for each duration
    imax=as.numeric(apply(mat,2,which.max))
    points(imax,1:length(s$windows),pch=16,cex=0.6)
    # Adding the first 10 local max
    if (show.events>1){
      eve=get.max.events(s,stat.max=stat,years=y,only.max=(show.events==1),
                         min.value=min.value,min.win=min.win,max.win=max.win,reverse.prob=reverse.prob)$events
      neve=min(nrow(eve),show.events)
#      if (neve>0) points(eve$id[1:neve],eve$iw[1:neve],pch="X",font=2,cex=seq(1.5,0.5,le=10))
      if (neve>0) text(eve$id[1:neve],eve$iw[1:neve],1:neve,font=2)
      
    }
    if (!is.null(min.win)) abline(h=which(s$windows==min.win)-0.5,lty=2)
    if (!is.null(max.win)) abline(h=which(s$windows==max.win)+0.5,lty=2)
  }
}

#---------------------------------------------------------------------------------
# Figure local goodness of fit (1d)
#---------------------------------------------------------------------------------
myfig.gof=function(s=NULL,var="TM",meth="calend",loc="France",set="mf",
                   y=2003,id=1,iw=1){
  # Get data
  if (is.null(s)) s=get.scan.1d(var,meth,loc,set)
  var=s$var; meth=s$method; loc=s$location; set=s$set
  dat=s$data-s$trend+rep.abind(s$trend[,,which(s$period==y)],length(s$period))
  fit=s$scan[[paste0("Y",y)]]$fit1
  nx=dim(fit)[1]
  ny=dim(fit)[2]
  # Compute ks test
  ksfun=function(x) mykstest(x[-(1:3)],s$params$distrib,x[1:3])
  ks=apply(abind(fit,dat,along=3),1:2,ksfun)
#  if (s$params$annual.maxima){ for (i in 1:ncol(ks)) ks[,i]=ks[1,i]}
  # Replace NAs in fit for plot
  fit[which(is.na(fit))]=-999
  names=c("mu","sigma","ksi")
  # Panel plot
  layout(matrix(1:4,2,2))
  par(mar=c(2,2,2,1),mgp=c(1,0.1,0),las=1,tcl=0.3)
  for (i in 1:3) {
    if (s$params$annual.maxima){
      plot(fit[1,,i],type="l",lwd=2,main=names[i],xlab="",ylab="")
      abline(v=iw)
    } ; if (!s$params$annual.maxima){
      image.plot(1:nx,1:ny,fit[,,i],axes=F,main=names[i],xlab="",ylab="")
      points(id,iw,pch="X",font=2,cex=1.5)
    }
  }
  if (s$params$annual.maxima){
    plot(ks[1,],type="l",lwd=2,main="ks-test p-value",xlab="",ylab="",ylim=c(0,1))
    abline(v=iw); abline(h=c(0,0.05,0.1),lty=1:3)
  } ; if (!s$params$annual.maxima){
    image.plot(1:nx,1:ny,ks,axes=F,main="ks-test p-value",breaks=c(0,0.05,seq(0.1,1,0.1)),
               col=c("red","orange",brewer.pal(9,"GnBu")))
    points(id,iw,pch="X",font=2,cex=1.5)
  }
}

#---------------------------------------------------------------------------------
# Plots the history of most extreme events
# JJA: xlim=c(152,243)
# DJF: xlim=c(335,59) --> D year y-1 and JF year y
# July to June: xlim=c(182,181)
#---------------------------------------------------------------------------------
myfig.chronology=function(e=NULL,var='TM',meth='calend',loc='France',set='mf',
                          stat.max="p1",years=NULL,only.max=FALSE,
                          min.value=0.95,min.win=1,max.win=90,reverse.prob=F,
                          stat="p1",lev=BRKS0,col=COLSM,title="",yaxis=T,
                          xlim=c(1,365),show.events=0,image=T,colorbar=T,output=F){
  # Get events
  if (is.null(e))  e=get.max.events(var=var,meth=meth,loc=loc,set=set,stat.max=stat.max,
                                    years=years,only.max=only.max,min.value=min.value,
                                    min.win=min.win,max.win=max.win,reverse.prob=reverse.prob)
  var=e$var; meth=e$method; loc=e$location; set=e$set
  # Re-order wrt stat
  ee=e$events
  ee=ee[order(ee[,paste(stat)],decreasing=T),]
  # Limits of the plot
  xx=xlim; if (xlim[2]<xlim[1]) xx[2]=xx[2]+365
  yy=range(e$scan.years)+c(-2,2) ; if (xlim[2]<xlim[1]) yy[1]=yy[1]+1
  # Function that places a day in the plot
  id2dy=function(id){
    d=mymod(id,365); y=e$years[1]+(id-d)/365
    if (xlim[2]<xlim[1]){
      if (d>xlim[2]) y=y+1
      if (d<=xlim[2]) d=d+365
    }; c(d,y)
  }
  # Value to plot
  ee$val=ee[,paste(stat)]; if (stat %in% c("p1","p0")) ee$val=100*ee$val
  # Panel plot (image + colorbar)
  # Warning: colorbar first to allow for click on image in shiny app
  if (image & colorbar) layout(matrix(2:1,1,2),width=c(1,0.2))
  # Plot 1: colorbar
  if (colorbar){
    par(mar=c(3,0,1.5,1))
    levleg=lev[2:(length(lev)-1)]; if (stat %in% c("p1","p0")) levleg=perc2prob(levleg)
    mycolorbar(col,levleg,horiz=F,width=0.2,cex=1)
  }
  # Plot 2: image
  if (image){
    par(mar=c(3,3-2.5*(!yaxis),1.5,0.5),mgp=c(1.5,0.1,0),las=1,tcl=0.5)
    plot(0,type="n",xaxs="i",yaxs="i",xlim=xx,ylim=yy,
          axes=F,frame.plot=T,main="",xlab="Calendar days",ylab="")
    title(ylab="Years",mgp=c(2,0,0))
    # Part of the axes here
    for (i in c(1,3)) my.calendar.axis(i,lab=(i==1))
    ylab=seq(1950,2020,5); abline(h=ylab,lty=3,col='gray')
    # Prepare output if requested (for shiny click)
    if (output){
      xxx=xx[1]:xx[2]; yyy=yy[1]:yy[2]
      out=matrix(NA,nrow=length(xxx),ncol=length(yyy))
    }
    # Plot rectangles: loop on events in reversed order for text of top events
    for (i in nrow(ee):1){ ei=ee[i,]
      # coordinates
      dy0=id2dy(ei$id); d0=dy0[1]; y0=dy0[2]
      d1=id2dy(ei$day1)[1]; d2=id2dy(ei$day2)[1]
      # color
      ci=col[as.numeric(cut(ei$val,lev))]
      # rectangle
      if (d1<=d2) rect(d1-0.5,y0-0.5,d2+0.5,y0+0.5,col=ci)
      if (d1>d2){
        rect(d1-0.5,y0-0.5-(d1>d0),xx[2]+0.5,y0+0.5-(d1>d0),col=ci)
        rect(xx[1]-0.5,y0+0.5-(d1>d0),d2+0.5,y0+1.5-(d1>d0),col=ci)
      }
      #segments(ei$d0,ei$y0-0.5,ei$d0,y0+0.5,lwd=2)
      # top events
      if (i<=show.events) text(d1,y0-0.5,i,adj=c(1.3,0),font=2,cex=1.5)
      # fill output table
      if (output){
        if (d1<=d2) out[which(xxx %in% d1:d2),which(yyy==y0)]=i
        if (d1>d2){
          out[which(xxx %in% d1:xx[2]),which(yyy==(y0-(d1>d0)))]=i
          out[which(xxx %in% xx[1]:d2),which(yyy==(y0+1-(d1>d0)))]=i
        }
      }
    }
    # Rest of the axes
    axis(2,at=ylab,lab=yaxis & ylab); axis(4,at=ylab,lab=F)
    # Titles
    if (!is.null(title)) title(paste(title,stat,var,loc),adj=0)
    title(method.name(meth),adj=1,col.main=method.color(meth))
    # Output: return info for clickable graph
    if (output) return(list(x=xxx,y=yyy,index=out,events=ee))
  }
}

#---------------------------------------------------------------------------------
# A base for a Robinson map with grid, axes and titles
#---------------------------------------------------------------------------------
my.map.base=function(crs=crs.robinson$crs,
                     map=ne_coastline("medium","sf"),col.map="gray",
                     lon.lines=seq(-180,180,30),lat.lines=seq(-90,90,15),col.lines="gray",
                     lon.labels=seq(-120,120,60),lat.labels=seq(-60,60,30),
                     xline=-180,yline=-90,xadj=c(0.5,1.3),yadj=c(1.3,0.5),cex.label=0.75,
                     title=c("","",""),cex.title=rep(1.2,3),font.title=rep(2,3),
                     col.title=rep(1,3),adj.title=c(0,0.5,1)){
  # Projection
  if (is.null(crs) | crs=="") crs=st_crs(map)
  # Plot coast lines
  plot(st_geometry(st_transform(map,crs)),col=col.map)
  # Add lon-lat lines
  lines=graticule(lons=lon.lines,lats=lat.lines,proj=crs)
  edges=graticule(lons=range(lon.lines),lats=range(lat.lines),proj=crs)
  plot(lines,col=col.lines,add=TRUE); plot(edges,add=TRUE)
  # Add lon_lat labels
  labs=graticule_labels(lons=lon.labels,lats=lat.labels,xline=xline,yline=yline,proj=crs)
  text(subset(labs,labs$islon),lab=parse(text=labs$lab[labs$islon]),adj=xadj,cex=cex.label,xpd=NA)
  text(subset(labs,!labs$islon),lab=parse(text=labs$lab[!labs$islon]),adj=yadj,cex=cex.label,xpd=NA)
  # Add titles
  tits=graticule_labels(lons=c(min(lon.lines),mean(range(lon.lines)),max(lon.lines)),yline=max(lat.lines),proj=crs)
  for (i in 1:3) text(subset(tits,tits$islon)[i,],title[i],cex=cex.title[i],font=font.title[i],
                      col=col.title[i],adj=c(adj.title[i],-0.5),xpd=NA)
  # Returns arguments
  list(crs=crs,map=map,lon.lines=lon.lines,lat.lines=lat.lines,lon.labels=lon.labels,lat.labels=lat.labels,
       lines=lines,edges=edges,labs=labs,tits=tits)
}

#---------------------------------------------------------------------------------
# Figure p1max all domains (2d)
#---------------------------------------------------------------------------------
myfig.2d=function(e,stat="p1",lev=BRKSM,col=COLSM,title="",ltitle=T,rtitle=T,
                  image=T,colorbar=T,reverse.prob=F,crs=crs.robinson$crs){
  # Generate image
  if (image){
    # Projection
    if (is.null(crs) | crs=="") crs=st_crs(map)
    # Get map shapefile
    wraf=st_read(paste0(WRAFDIR,"shapefiles/region_fx-WRAF",e$wraf.set,"-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.shp"))
    wraf=suppressWarnings(st_set_crs(wraf,4326)) # WGS84
    wraf=st_transform(wraf,4326)
    # Get values to plot
    vals=as.numeric(e$events[,paste(stat)])
    if (stat %in% c("p1","p0")) vals=100*vals
    # Vector of colors corresponding to stat values
    cols=col[as.numeric(cut(vals,lev))]
  }
  # Plot image + colorbar
  # Warning: colorbar first to allow for click on image in shiny app
  if (image & colorbar) layout(matrix(2:1,1,2),width=c(1,0.15))
  # Plot 1: colorbar
  if (colorbar){
    par(mar=c(0.5,0,1.5,1))
    levleg=lev[2:(length(lev)-1)]; if (stat %in% c("p1","p0")) levleg=perc2prob(levleg)
    mycolorbar(col,levleg,horiz=F,width=0.2,cex=1)
  }
  # Plot 2: image
  if (image){
    par(mar=c(0.5,1.5,1.5,0.5))
    ltit=""; if (!is.null(title)){ ltit=title ; if (ltitle){
      yrs=e$selparams$years ; if (is.null(yrs)) yrs=e$scan.years
      dts=paste0(unique(range(yrs)),collapse="-")
      if (!is.null(e$selparams$only.months)){ mm=e$selparams$only.months
        if (length(e$selparams$only.months)==1) dts=paste(dts,MTHS[e$selparams$only.months])
        if (length(e$selparams$only.months)>1)  dts=paste(dts,paste0(MTHS1[e$selparams$only.months],collapse=""))
      }
      if (!is.null(e$selparams$id) & !is.null(e$selparams$iw))
        dts=paste(dts,idays2leg(e$events$day1[1],e$events$day2[2],short=T))
      ltit=paste(ltit,stat,e$var,dts)
    }}
    rtit=""; if (rtitle) rtit=method.name(e$method)
    adj.title=c(0,0.5,1); if (crs==crs.robinson$crs) adj.title=rep(0.5,3)
    map=my.map.base(crs=crs,title=c(ltit,"",rtit),col.title=c(1,1,method.color(e$method)),
                    adj.title=adj.title)
    plot(st_transform(wraf["ID"],map$crs),col=cols,add=T)
  }
}

#---------------------------------------------------------------------------------
# Figure p1max all domains (2d)
#---------------------------------------------------------------------------------
myfig.2d.ggplot=function(e,stat="p1",lev=BRKSM,col=COLSM,title="",ltitle="default",rtitle="default",
                         image=T,colorbar=T,reverse.prob=F,crs=crs.robinson){
  # Generate image
  if (image){
    # Get map shapefiles
    coast=ne_coastline("medium","sf")
    wraf=st_read(paste0(WRAFDIR,"shapefiles/region_fx-WRAF",e$wraf.set,"-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.shp"))
    wraf=suppressWarnings(st_set_crs(wraf,4326)) # WGS84
    wraf=st_transform(wraf,4326)
    # Get values to plot
    ee=e$events
    ee[c("p1","p0")]=100*ee[c("p1","p0")]
    ee=ee[,c(paste(stat),"loc")]
    colnames(ee)=c("value","Identifier")
    # Merge map and values
    colnames(ee)[which(colnames(ee)=="loc")]="Identifier"
    data=merge(wraf,ee,by="Identifier")
    data=data["value"]
    # Vector of colors corresponding to stat values
    cols=col[as.numeric(cut(c(data$value),lev))]
  }
  # Plot image + colorbar
  # Warning: no layout for basic plot and ggplot ... I use cowplot
  # Plot 1: colorbar
  if (colorbar){
    levleg=lev[2:(length(lev)-1)]; if (stat %in% c("p1","p0")) levleg=perc2prob(levleg)
    if (!image) mycolorbar(col,levleg,horiz=F,width=0.2,cex=1)
    if (image){
      plot1=function(){
        par(mar=c(0.5,0,1.5,1))
        mycolorbar(col,levleg,horiz=F,width=0.2,cex=1)
      }
    }
  }
  # Plot 2: map with ggplot2
  if (image){
    if (ltitle=="default"){
      yrs=e$selparams$years ; if (is.null(yrs)) yrs=e$scan.years
      dts=paste0(unique(range(yrs)),collapse="-")
      if (!is.null(e$selparams$only.months)){ mm=e$selparams$only.months
      if (length(e$selparams$only.months)==1) dts=paste(dts,MTHS[e$selparams$only.months])
      if (length(e$selparams$only.months)>1)  dts=paste(dts,paste0(MTHS1[e$selparams$only.months],collapse=""))
      }
      if (!is.null(e$selparams$id) & !is.null(e$selparams$iw))
        dts=paste(dts,idays2leg(e$events$day1[1],e$events$day2[2],short=T))
      ltitle=paste(stat,e$var,dts)
    }
    rtitlecol=1
    if (rtitle=="default") {rtitle=method.name(e$method) ; rtitlecol=method.color(e$method)}
    theme_set(theme_bw())
    plot2=ggplot() +
      geom_sf(data=st_geometry(coast),col="gray") +
      geom_sf(data=data,fill=cols,col=1) +
      #geom_sf(data,aes(fill=value),col=1) +
      coord_sf(xlim=crs$xlim,ylim=crs$ylim,crs=crs$crs) +
      labs(title=paste(title,ltitle),tag=rtitle,x="",y="",fill=stat) +
      theme(plot.title=element_text(hjust=0,size=15,face="bold"),
            plot.tag=element_text(color=rtitlecol,size=15,face="bold",hjust=1,vjust=1),
            plot.tag.position=c(1,1))
    if (!colorbar) return(plot2)
  }
  # Combine image + colorbar with cowplot
  if (image & colorbar) plot_grid(plot2,plot1,rel_widths=c(1,0.15))
}


##############################################
# Additional stuff.. work in progress..
##############################################

myfig.all.fits=function(var='TX',meth='annmax',wset=2){
  mu=sigma=ksi=ks=c();  for (loc in REGIONS[[paste0("WRAF",wset)]][1,]){
    s=get.scan.1d(var,meth,loc,dat='era5')
    y=s$scan.years[1]
    dat=s$data-s$trend+rep.abind(s$trend[,,which(s$period==y)],length(s$period))
    fit=s$scan[[paste0("Y",y)]]$fit1
    ksfun=function(x) mykstest(x[-(1:3)],s$params$distrib,x[1:3])
    ksl=apply(abind(fit,dat,along=3),1:2,ksfun)
    fit=apply(fit,2:3,mean,na.rm=T)
    mu=cbind(mu,fit[,1]); sigma=cbind(sigma,fit[,2]); ksi=cbind(ksi,fit[,3])
    ks=cbind(ks,apply(ksl,2,mean,na.rm=T))
  }
  list(locs=REGIONS[[paste0("WRAF",wset)]][1,],mu=mu,sigma=sigma,ksi=ksi,ks=ks)
}

check.ksi=function(wset=2,thres=0.5){
  out=c()
  for (loc in REGIONS[[paste0("WRAF",wset)]][1,]){
    SCAN=get.scan.1d('TX','annmax',loc,dat='era5')
    ksi=SCAN$scan$Y2001$fit1[1,,3]
    ii=which(abs(ksi)>thres)
    if (length(ii)>0){ for (i in ii) out=rbind(out,c(loc,i,ksi[i]))}
  }
  out
}

hist.all.stat.1d=function(var="TM",meth="calend",loc="France",set="mf",stat="p1",
                          years=NULL,id=NULL,iw=1,reverse.prob=FALSE){
  # Get all stats from scan 1d
  s=get.scan.1d(var,meth,loc,set)
  yrs=s$scan.years; if (!is.null(years)) yrs=years
  ss=c(); for (y in yrs) ss=abind(ss,as.matrix(s$scan[[paste0("Y",y)]][[paste(stat)]]),along=3)
  if (!is.null(id)) ss=array(ss[id,,],dim=c(length(id),dim(ss)[2],dim(ss)[3]))
  if (!is.null(iw)) ss=array(ss[,iw,],dim=c(dim(ss)[1],length(iw),dim(ss)[3]))
  if (meth %in% c("annmax","annmin")) ss=apply(ss,c(2,3),max,na.rm=T)
  # Returns a vector
  out=c(ss)
  if (reverse.prob) out=1-out
  out
}

hist.all.stat.2d=function(var="TM",meth="annmax",set="era5",wset=2,stat="p1",
                          years=NULL,id=NULL,iw=1:14,reverse.prob=FALSE){
  # Get all stats from scans 1d
  out=c()
  for (loc in REGIONS[[paste0("WRAF",wset)]][1,]){
    s=get.scan.1d(var,meth,loc,set)
    yrs=years; if (is.null(years)) yrs=s$scan.years
    for (y in yrs){
      sy=as.matrix(s$scan[[paste0("Y",y)]][[paste(stat)]])
      if (!is.null(id)) sy=as.matrix(sy[id,],nrow=length(id),ncol=ncol(sy))
      if (!is.null(iw)) sy=as.matrix(sy[,iw],nrow=nrow(sy),ncol=length(iw))
      if (meth %in% c("annmax","annmin")) sy=apply(sy,2,max,na.rm=T)
      out=c(out,sy)
    }
  }
  # Returns a vector
  if (reverse.prob) out=1-out
  out
}

hist.all.year.2d=function(var="TM",meth="annmax",set="era5",wset=2,years=NULL)
  get.scan.2d(var=var,meth=meth,set=set,wset=wset,years=years)$events$y

get.yearly.stat.2d=function(var="TM",meth="annmax",set="era5",wset=2,stat="p1",reverse.prob=F){
  yrs=get.scan.1d(var,meth,REGIONS[[paste0("WRAF",wset)]][1,1],set)$scan.years
  out=c()
  for (y in yrs)
    out=cbind(out,get.scan.2d(var,meth,set,wset,years=y)$events[,paste(stat)])
  if (reverse.prob) out=1-out
  out
}

get.all.trends=function(var="TM",meth="calend",set="era5",wset=2,id=1:365,iw=1){
  # Get all stats from scans 1d
  locs=REGIONS[[paste0("WRAF",wset)]][1,]
  out=c()
  for (loc in locs){
    s=get.scan.1d(var,meth,loc,set)
    fr=s$forced.response
    tr=apply(as.matrix(s$trend[id,iw,],nrow=length(id)),2,mean,na.rm=T)
    out=c(out,lm(tr~fr)$coef[2])
  }
  # Returns a df
  data.frame(loc=locs,val=out)
}

get.all.fits=function(var="TM",meth="calend",set="era5",wset=2,year=2022,id=1:365,iw=1){
  # Get all stats from scans 1d
  locs=REGIONS[[paste0("WRAF",wset)]][1,]
  out=c()
  for (loc in locs){
    s=get.scan.1d(var,meth,loc,set)$scan[[paste0("Y",year)]]
    fit=apply(as.matrix(s$fit1[id,iw,],nrow=length(id)),2,mean,na.rm=T)
    out=rbind(out,fit)
  }
  # Returns a df
  data.frame(loc=locs,mu=out[,1],sig=out[,2],ksi=out[,3])
}

extra.map=function(df,name,lev=NULL,col=NULL,title="",ltitle=T,rtitle=T,
                   image=T,colorbar=T,crs=crs.robinson$crs){
  wset=wraf.set(df$loc[1])
  s2d=get.scan.2d(wset=wset,years=2022)
  s2d$events[,paste(name)]=df[,paste(name)]
  if (is.null(lev)) lev=pretty(df[,paste(name)],10)
  if (is.null(col)){
    nlev=length(lev)
    nblue=floor(nlev/2)
    nreds=ceiling(nlev/2)-1
    col=c(rev(brewer.pal(nblue,"GnBu")),brewer.pal(nreds,"YlOrRd"))
  }
  myfig.2d(e=s2d,stat=name,lev=lev,col=col,title=title,ltitle=ltitle,rtitle=rtitle,
                  image=image,colorbar=colorbar,crs=crs)
}
  
  
