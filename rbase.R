###########################################
# Basic functions used for other scripts  #
# + loading of usual packages             #
###########################################

#-------------------------------------------------------------------------------------
# Packages
#-------------------------------------------------------------------------------------
library(ncdf4,quietly=TRUE,verbose=FALSE)
#library(clim.pact,quietly=TRUE,verbose=FALSE)
library(fields,quietly=TRUE,verbose=FALSE)
library(abind,quietly=TRUE,verbose=FALSE)
library(maps,quietly=TRUE,verbose=FALSE)
library(mapproj,quietly=TRUE,verbose=FALSE)
#library(sn,quietly=TRUE,verbose=FALSE)
#library(class,quietly=TRUE,verbose=FALSE)
#library(mclust,quietly=TRUE,verbose=FALSE)
#library(grImport,quietly=TRUE,verbose=FALSE)
library(chron,quietly=TRUE,verbose=FALSE)

#-------------------------------------------------------------------------------------
# Vectors/time-series related functions (x)
#-------------------------------------------------------------------------------------
# Basic functions to get/remove first/last/penul(n-1) element of a vector
first=function(x)   x[1]
nofirst=function(x) x[-1]
last=function(x)    x[length(x)]
nolast=function(x)  x[-length(x)]
penul=function(x)   x[length(x)-1]

# Difference between consecutive elements of a vector (eg t+1/t for time series)
autodiff=function(x) {
  out=numeric(0); if (length(x)>1) out=x[2:length(x)]-x[1:(length(x)-1)]
  out
}

# Function replacing vector values persisting less than n times by NAs (default)
mypersist=function(x,n,replace=NA){
  chg=c(0,which(is.na(autodiff(x)) | autodiff(x)!=0),length(x))
  na=which(autodiff(chg) < n)
  for (i in 1:length(na)) x[(chg[na[i]]+1):chg[na[i]+1]]=replace
  x
}

# Normalization of a time series
mynorm=function(x,mean=NULL,sd=NULL){
  if (is.null(mean)) mean=mean(x,na.rm=TRUE); if (is.null(sd)) sd=sd(x,na.rm=TRUE)
  (x-mean)/sd
}

# Biased formula for variance
mybiasedvar=function(x,na.rm=TRUE){
  if (na.rm) out=mybiasedvar(na.omit(x),na.rm=FALSE) else out=mean((x-mean(x))^2)
  out
}

# Computing percentages
mypercent=function(x) x*100/sum(x)

# Slope of the linear regression of a vector
mytrend=function(x){
  out=NA; if (any(!is.na(x))) out=lm(x~c(1:length(x)))$coef[2]
  out
}

# Linear fit of a vector
mylmfit=function(x,ref=NULL,predict.na=FALSE){
  if (is.null(ref)) ref=1:length(x)
  if (any(!is.na(x))){ fit=lm(x~ref,na.action=na.exclude)
    if (predict.na)  out=as.vector(predict(fit,as.data.frame(ref)))
    if (!predict.na) out=as.vector(predict(fit))
  }
  out
}

# Detrending a time series (wrt ref, only if significant wrt pval)
mydetrend=function(x,ref=NULL,pval=0.05){
  if (is.null(ref)) ref=1:length(x)
  if (mypvalue(x,ref)<pval) out=x-mylmfit(x~ref)
  else out=x-mean(x)
  out
}

# Splines accounting for NAs + possibility of deducing df from length(x) (sdf=length(x)/df)
# (Use smooth.na=F for not predicting at NA values of x)
mysmoothy=function(x,df=NULL,sdf=30,smooth.na=TRUE,periodic=FALSE){
  if (is.null(df) & is.null(sdf)) sdf=length(na.omit(x))
  if (is.null(df)) df=length(na.omit(x))/sdf
  if (length(unique(na.omit(x))) < max(4,df)) out=x#rep(NA,length(x))
  else {
    if (df>1){
      if (periodic) { x=rep(x,3) ; df=df*3 }
      sms=smooth.spline(which(!is.na(x)),na.omit(x),df=df)
      out=predict(sms$fit,1:length(x))$y
      if (!smooth.na) out[which(is.na(x))]=NA
      if (periodic) { out=matrix(out,ncol=3)[,2] ; df=df/3}
    }
    if (df<=1 & !periodic) out=mylmfit(x,predict.na=smooth.na)
    if (df<=1 & periodic){ out=rep(mean(x,na.rm=T),length(x)); if (!smooth.na) out[which(is.na(x))]=NA}
  }
  out
}

# Running averages
myrunmean=function(x,n,na.rm=T,fill.with.NAs=F,align="center"){
  out=c(); for (i in 1:(length(x)-n+1)) out=c(out,mean(x[i:(i+n-1)],na.rm=na.rm))
  if (fill.with.NAs){
    if (align=="center") out=c(rep(NA,floor((n-1)/2)),out,rep(NA,ceiling((n-1)/2)))
    if (align=="left")   out=c(out,rep(NA,n-1))
    if (align=="right")  out=c(rep(NA,n-1),out)
  }
  out
}

# Average by blocks (useful for seasonal averages, for instance)
myblockmean=function(x,n,na.rm=TRUE){
  if (length(x)/n != trunc(length(x)/n)) print("Error: length(x) must be a multiple of n.")
  else return(apply(matrix(x,nrow=n),2,mean,na.rm=na.rm))
}

myblockmax=function(x,n,na.rm=TRUE){
  if (length(x)/n != trunc(length(x)/n)) print("Error: length(x) must be a multiple of n.")
  else return(apply(matrix(x,nrow=n),2,max,na.rm=na.rm))
}

myblockmin=function(x,n,na.rm=TRUE){
  if (length(x)/n != trunc(length(x)/n)) print("Error: length(x) must be a multiple of n.")
  else return(apply(matrix(x,nrow=n),2,min,na.rm=na.rm))
}

#-------------------------------------------------------------------------------------
# Functions for comparing two vectors/time-series (x,y)
#-------------------------------------------------------------------------------------
# Euclidean distance between two vectors
myeuclid=function(x,y,na.rm=TRUE) sqrt(sum((x-y)^2,na.rm=na.rm))

# RMS difference between two vectors
myrmsd=function(x,y,na.rm=TRUE) sqrt(mean((x-y)^2,na.rm=na.rm))

# Mahalanobis distance between two populations (samples in lines)
mymaha=function(x,y,sym=TRUE){
  mux=apply(x,2,mean,na.rm=TRUE)
  muy=apply(y,2,mean,na.rm=TRUE)
  if (!sym)  out=(t(mux-muy)) %*% (solve(cov(x))) %*% (mux-muy)
  if (sym) out=0.5*(mymaha(x,y,sym=FALSE)+mymaha(y,x,sym=FALSE))
  out
}

# Correlation function accounting for NAs
mycor=function(x,y) {
  ix=which(is.na(x)); iy=which(is.na(y)); ixy=unique(c(ix,iy))
  if (length(ixy)>0) {x=x[-ixy];y=y[-ixy]}
  cor(x,y)
}

# P.value of a linear trend (if y=NULL) or a correlation test with y
mypvalue=function(x,y=NULL){
  if (is.null(y)) y=1:length(x)
  ix=which(is.na(x)); iy=which(is.na(y)); ixy=unique(c(ix,iy))
  if (length(ixy)>0) {x=x[-ixy];y=y[-ixy]}
  if (length(x)>2) return(cor.test(x,y)$p.value) else return(NA)
}

#-------------------------------------------------------------------------------------
# Weighted functions (eg for cos(latitude) weighting)
# Edit: IC 95% with bootstrap
#-------------------------------------------------------------------------------------
# Weighted mean sd var and skewness
mywmean=function(x,w,na.rm=TRUE,bootstrap=FALSE){
  nas=which(is.na(x)); if (length(nas)>0 & na.rm) {x=x[-nas]; w=w[-nas]}
  if (length(x)==0) return(NA)
  if (length(w)==1) w=rep(w,length(x))
  w=w/sum(w); out=sum(x*w)
  if (bootstrap) {tmp=c(); for (i in 1:1000){
    ii=sample(1:length(x),length(x),replace=TRUE); ww=w[ii]; xx=x[ii]
    ww=ww/sum(ww); tmp=c(tmp,sum(xx*ww))
  }; out=c(out,quantile(tmp,0.025),quantile(tmp,0.975))}
  out
}
mywsd=function(x,w,na.rm=TRUE,unbiased=TRUE,bootstrap=FALSE){
  nas=which(is.na(x)); if (length(nas)>0 & na.rm) {x=x[-nas]; w=w[-nas]}
  if (length(x)==0) return(NA)
  n=length(x); unbias=1; if (unbiased) unbias=n/(n-1)
  mx=mywmean(x,w); out=sqrt(unbias*mywmean((x-mx)^2,w))
  if (bootstrap) {tmp=c(); for (i in 1:1000){
    ii=sample(1:length(x),length(x),replace=TRUE); ww=w[ii]; xx=x[ii]
    mxx=mywmean(xx,ww); tmp=c(tmp,sqrt(unbias*mywmean((xx-mxx)^2,ww)))
  }; out=c(out,quantile(tmp,0.025),quantile(tmp,0.975))}
  out
}
mywvar=function(x,w,na.rm=TRUE,unbiased=TRUE,bootstrap=FALSE){
  nas=which(is.na(x)); if (length(nas)>0 & na.rm) {x=x[-nas]; w=w[-nas]}
  if (length(x)==0) return(NA)
  n=length(x); unbias=1; if (unbiased) unbias=n/(n-1)
  mx=mywmean(x,w); out=unbias*mywmean((x-mx)^2,w)
  if (bootstrap) {tmp=c(); for (i in 1:1000){
    ii=sample(1:length(x),length(x),replace=TRUE); ww=w[ii]; xx=x[ii]
    mxx=mywmean(xx,ww); tmp=c(tmp,unbias*mywmean((xx-mxx)^2,ww))
  }; out=c(out,quantile(tmp,0.025),quantile(tmp,0.975))}
  out
}
mywskew=function(x,w,na.rm=TRUE,unbiased=TRUE,bootstrap=FALSE){
  nas=which(is.na(x)); if (length(nas)>0 & na.rm) {x=x[-nas]; w=w[-nas]}
  if (length(x)==0) return(NA)
  n=length(x); unbias=1; if (unbiased) unbias=n^2/((n-1)*(n-2))
  mx=mywmean(x,w); sx=mywsd(x,w,unbiased=unbiased); out=unbias*mywmean(((x-mx)/sx)^3,w)
  if (bootstrap) {tmp=c(); for (i in 1:1000){
    ii=sample(1:length(x),length(x),replace=TRUE); ww=w[ii]; xx=x[ii]
    mxx=mywmean(xx,ww); sxx=mywsd(xx,ww,unbiased=unbiased);
    tmp=c(tmp,unbias*mywmean(((xx-mxx)/sxx)^3,ww))
  }; out=c(out,quantile(tmp,0.025),quantile(tmp,0.975))}
  out
}

# Weighted quantiles
mywqs=function(x,w,p=myquantiles,na.rm=TRUE){
  nas=which(is.na(x)); if (length(nas)>0 & na.rm) {x=x[-nas]; w=w[-nas]}
  if (length(x)==0) return(rep(NA,length(p)))
  x=x[order(x)]; w=w[order(x)]; out=c()
  for (pp in p) out=c(out,x[which.max(cumsum(w)/sum(w)>=pp)])
  out
}

# Weighted Euclidean distance
myweuclid=function(x,y,w,na.rm=TRUE) sqrt(sum(w*(x-y)^2,na.rm=na.rm))

# Weighted RMS difference
mywrmsd=function(x,y,w,na.rm=TRUE) sqrt(mywmean((x-y)^2,w,na.rm=na.rm))

# Weighted correlation
mywcor=function(x,y,w,na.rm=T){
  nas=which(is.na(x) | is.na(y))
  if (length(nas)>0 & na.rm) {x=x[-nas]; y=y[-nas]; w=w[-nas]}
  w=w/sum(w) ; x=x-sum(w*x) ; y=y-sum(w*y)
  sum(w*x*y)/sqrt(sum(w*x*x)*sum(w*y*y))
}

# Weighted coordinates for Taylor diagrams (x0 = reference)
mywtaylorcoords=function(x0,x,w) {
  sig0=mywsd(x0,w); sig=mywsd(x,w); r=mywcor(x0,x,w)
  sig/sig0*c(r,sqrt(1-r^2))
}

# Weighted nearest neighbours of dat among cls - (eg for weather regimes attribution step)
# Input: dat matrix ntime*nspace - cls matrix ncluster*nspace - weigths vector nspace
mywnn=function(dat,cls,weights,method=c("rmsd","cor","euclid")){
  meth=method[1]; nt=nrow(dat); ns=ncol(dat); nc=nrow(cls)
  if (ncol(cls)!=ns | length(weights)!=ns) return("Error in input dimensions")
  if (meth=="rmsd")   { fun=function(x,y) apply(x,1,mywrmsd,y,weights)   ; best=min }
  if (meth=="cor")    { fun=function(x,y) apply(x,1,mywcor,y,weights)    ; best=max }
  if (meth=="euclid") { fun=function(x,y) apply(x,1,myweuclid,y,weights) ; best=min }
  out=c(); for (i in 1:nt){
    vals=fun(cls,dat[i,]); icls=first(which(vals==best(vals)))
    out=c(out,icls,vals)
  }
  out=as.data.frame(matrix(out,nrow=nt,byrow=TRUE)); names(out)=c("class",paste(meth,1:nc,sep=""))
  out
}

#-------------------------------------------------------------------------------------
# Bootstrap functions
#-------------------------------------------------------------------------------------
# Mean of a time series with associated IC
mybootstrap=function(x,n=1000,p=0.05,nmin=length(x),na.rm=TRUE,replace=TRUE){
  if (nmin==length(x)) size=rep(length(x),n) else size=sample(nmin:length(x),n,replace=TRUE)
  boot=c(); for (i in 1:n){
    xi=x[sample(1:length(x),size[i],replace=replace)]
    boot=c(boot,mean(xi,na.rm=na.rm))
  }
  c(mean(boot),quantile(boot,p/2),quantile(boot,1-p/2))
}

# Mean of a difference between 2 time series with associated IC
mybootstrapdiff=function(x,y,n=1000,p=0.05,nmin=length(x),na.rm=TRUE,replace=TRUE){
  if (nmin==length(x)) size=rep(length(x),n) else size=sample(nmin:length(x),n,replace=TRUE)
  boot=c();  for (i in 1:n){
    xi=x[sample(1:length(x),size[i],replace=replace)]
    yi=y[sample(1:length(y),size[i],replace=replace)]
    boot=c(boot,mean(xi-yi,na.rm=na.rm))
  }
  c(mean(boot),quantile(boot,p/2),quantile(boot,1-p/2))
}

#-------------------------------------------------------------------------------------
# Array-related functions
#-------------------------------------------------------------------------------------
# Repetition of an array (x,y,t) along t, until length.out (if 2D-array, force 3D)
rep.abind=function(tab,length.out){
  if (length(dim(tab))==2){
    tmp=array(dim=c(dim(tab),1)); tmp[,,1]=tab; tab=tmp; rm(tmp)}
  d=dim(tab)[3]
  while (d<length.out){
    tab=abind(tab,tab,along=3)
    d=dim(tab)[3]
  }
  tab[,,1:length.out]
}

#-------------------------------------------------------------------------------------
# Time/Date-related functions - Dates are data.frames $m $d $y or yyyymmdd numbers
# Warning to CDO or NCO conventions for dates of monthly means (CDO: last day - NCO: median day)
#-------------------------------------------------------------------------------------

# Possible calendar types (Warning: modify ndays() below for weird calendars (eg. 360_day))
leapCaltypes=c("leap","standard","gregorian","proleptic_gregorian")
noleapCaltypes=c("noleap","365_day","360_day")

# Is year y leap (bissextile)?
bissex=function(y,caltype="leap"){
  if (caltype %in% noleapCaltypes) biss=FALSE
  else  biss=(trunc(y/400)*400==y | (trunc(y/4)*4==y & trunc(y/100)*100!=y))
  biss
}

# Nb of days of month (m,y) depending on calendar type
ndays=function(m,y,caltype="leap"){
  if (caltype=="360_day") n=rep(30,12)[m]
  else n=c(31,28+bissex(y,caltype),31,30,31,30,31,31,30,31,30,31)[m]
  n
}

# Nb of days of season (s,y)
ndays.season=function(s,y,caltype="leap"){
  n=0; for (m in season2months(s)) n=n+ndays(m,y,caltype)
  n
}

# My modulo (mymod(n,n)=n) - default: modulo 12 (for months)
mymod=function(x,n=12){
  mm=c(); for (i in 1:length(x)){
    if (x[i]>=1) mm=c(mm,x[i]-n*trunc((x[i]-1)/n))
    if (x[i]<1)  mm=c(mm,x[i]+n*(trunc(-x[i]/n)+1)) }
  mm
}

# Season (eg "JJA") to index of first month ($i=6) and nb of months ($n=3), to months (6:8)
# or to cdo-months ("6,7,8")
season2i=function(season){
  ii=c(); nn=c();
  mm=c("JF","FM","MA","AM","MJ","JJ","JA","AS","SO","ON","ND","DJ")
  for (s in season){
    split=strsplit(s,"")[[1]]  
    ii=c(ii,which(mm==paste(split[1:2],collapse="")))
    nn=c(nn,length(split))
  }
  data.frame(i=ii,n=nn)
}
season2months=function(season)    { ii=season2i(season); mymod(seq(ii$i,by=1,le=ii$n)) }
season2cdomonths=function(season) { mm=season2months(season); paste(mm,collapse=",")}

# Function similar to julday() but dealing with all calendar types
# (Odd reference to 1/1/-4713 in julday...)
# EDIT 2014: switch from julday to julian BUT reference left to -4713 for no-leap cases
myjulday=function(m,d,y,caltype="leap"){
  if (! caltype %in% c(leapCaltypes,noleapCaltypes)) print("ERROR: Unknown calendar type!")
  refy=-4713
  if (caltype %in% leapCaltypes)
    out=julian(m,d,y) #julday(m,d,y)
  if (caltype %in% noleapCaltypes){
    Ndyr=sum(ndays(1:12,1,caltype))
    dts=matrix(c(m,d,y),ncol=3)
    out=apply(dts,1,function(x) Ndyr*(x[3]-refy-(x[3]>0))+(x[1]!=1)*sum(ndays(1:(x[1]-1),x[3],caltype))+x[2])
  }
  out
}

# Function similar to caldat() but dealing with all caltypes (inverse function of myjulday)
# EDIT 2014: switch from caldat to month.day.year BUT reference left to -4713 for no-leap cases
mycaldat=function(jday,caltype="leap"){
  if (! caltype %in% c(leapCaltypes,noleapCaltypes)) print("ERROR: Unknown calendar type!")
  refy=-4713
  if (caltype %in% leapCaltypes) {
    out=month.day.year(jday) #out=caldat(jday)
    names(out)=c("m","d","y")
  }
  if (caltype %in% noleapCaltypes){
    Ndyr=sum(ndays(1:12,1,caltype)); CNdyr=cumsum(ndays(1:12,1,caltype))
    y=floor(refy+jday/Ndyr)+(floor(refy+jday/Ndyr)>=0)
    diny=jday-Ndyr*(y-refy-(y>0))
    y=y-(diny==0); diny=diny+(diny==0)*Ndyr
    dts=matrix(c(diny,y),ncol=2)
    m=apply(dts,1,function(x) min(which(sort(c(CNdyr,x[1]))==x[1])))
    dts=cbind(m,dts)
    d=apply(dts,1,function(x) x[2]-(x[1]!=1)*sum(ndays(1:(x[1]-1),x[3],caltype)))
    out=list(m=m,d=d,y=y)
  }
  out
}

# Function creating dates (data.frame) from m,d,y
makedates=function(m=NULL,d=NULL,y=NULL,avg="day",conv="cdo",caltype="leap"){
  out=c()
  if (length(y)!=0){ for (yy in y){ M=m; if (length(m)==0) M=1:12; for (mm in M){
      D=d; if (length(d)==0) D=1:ndays(mm,yy,caltype)
      if (avg=="mon" & conv=="cdo") D=last(D)
      if (avg=="mon" & conv=="nco") D=ceiling(median(D))
      if (avg=="mon" & conv=="1")   D=1
      for (dd in D) out=c(out,mm,dd,yy)
  }}}
  if (length(out)!=0) {
    out=data.frame(matrix(out,ncol=3,byrow=TRUE)); names(out)=c("m","d","y") }
  out
}

# Function creating dates corresponding to season over period
# (if inside=TRUE and season=DJF, JF of first year and D of last year will be included)
makedates.season=function(season,period,avg="day",inside=FALSE,caltype="leap",conv="cdo"){
  # Getting indices corresp. to season
  is=season2i(season)
  msta=min(is$i); mend=max(is$i)+is$n[which(is$i==max(is$i))]-1
  # Computing of dates
  dts=c(); for (y in period){ for (s in 1:length(season)){
    issta=is$i[s]; isend=is$i[s]+is$n[s]-1
    if (isend<=12)
      dts=c(dts,c(t(as.matrix(makedates(m=issta:isend,y=y,avg=avg,caltype=caltype,conv=conv)))))
    else {
      if (inside) dts=c(dts,c(t(as.matrix(makedates(m=1:mymod(isend),y=y,avg=avg,caltype=caltype,conv=conv)))),c(t(as.matrix(makedates(m=issta:12,y=y,avg=avg,caltype=caltype,conv=conv)))))
      else dts=c(dts,c(t(as.matrix(makedates(m=issta:12,y=y,avg=avg,caltype=caltype,conv=conv)))),c(t(as.matrix(makedates(m=1:mymod(isend),y=y+1,avg=avg,caltype=caltype,conv=conv)))))
    }
  }}
  dts=data.frame(matrix(dts,ncol=3,byrow=TRUE)); names(dts)=c("m","d","y")
  # Removing duplicated dates
  if (any(duplicated(dts))) dts=dts[-which(duplicated(dts)),]
  dts
}

# Convert date $m$d$y to time index wrt date0
mdy2i=function(date,date0,avg="day",caltype="leap"){
  if (avg=="day") {
    fun=function(x,y) myjulday(x[1],x[2],x[3],caltype)-myjulday(y[1],y[2],y[3],caltype)+1
    out=apply(as.matrix(date),1,fun,as.matrix(date0))
  }
  if (avg=="mon") out=apply(as.matrix(date),1,function(x,y) 12*(x[3]-y[3])+x[1],as.matrix(date0))
  out
}

# Convert date $m$d$y to time values wrt units
mdy2time=function(date,units,caltype="leap"){
  what=strsplit(units," ")[[1]][1]
  date0=as.numeric(strsplit(strsplit(units," ")[[1]][3],"-")[[1]])[c(2,3,1)]
  fun=function(x,y)
    myjulday(x[1],x[2],x[3],caltype)-myjulday(y[1],y[2],y[3],caltype)
  diff=apply(as.matrix(date),1,fun,date0)
  if (what=="days")    out=diff
  if (what=="hours")   out=24*diff
  if (what=="minutes") out=1440*diff
  if (what=="seconds") out=86400*diff
  out
}

# Reverse function
time2mdy=function(time,units,caltype="leap"){
  what=strsplit(units," ")[[1]][1]
  date0=as.numeric(strsplit(strsplit(units," ")[[1]][3],"-")[[1]])[c(2,3,1)]
  m0=date0[1]; d0=date0[2]; y0=date0[3]
  fun=function(x)
    as.data.frame(mycaldat(x+myjulday(m0,d0,y0,caltype),caltype))
  if (what=="days")    out=fun(time)
  if (what=="hours")   out=fun(trunc(time/24))
  if (what=="minutes") out=fun(trunc(time/1440))
  if (what=="seconds") out=fun(trunc(time/86400))
  names(out)=c("m","d","y")
  out
}

# Convert date $m$d$y to yyyymmdd (and inverse function)
mdy2yyyymmdd=function(date) 10000*date$y+100*date$m+date$d
yyyymmdd2mdy=function(x){
  yy=trunc(x/10000); mm=trunc((x-10000*yy)/100)
  data.frame(m=mm,d=x-100*trunc(x/100),y=yy)
}

# Convert date $m$d$y to cdo-date "yyyy-mm-dd"
mdy2cdo=function(date){
  out=c(); for (i in 1:nrow(date)){
    yyyy=paste(c(rep(0,4-nchar(date$y[i])),date$y[i]),collapse="")
    mm=paste(c(rep(0,2-nchar(date$m[i])),date$m[i]),collapse="")
    dd=paste(c(rep(0,2-nchar(date$d[i])),date$d[i]),collapse="")
    out=c(out,paste(yyyy,mm,dd,sep="-"))
  }
  out
}

#----------------------------------------------------------------------
# NETCDF functions (new version with package ncdf4)
#----------------------------------------------------------------------
# Create netcdf
mync=function(file,dat,var,lon,lon.units,lat,lat.units,time,time.units,caltype=NA){
  if (length(lon)==0)  lon=0  ; if (length(lon.units)==0)  lon.units=""
  if (length(lat)==0)  lat=0  ; if (length(lat.units)==0)  lat.units=""
  if (length(time)==0) time=0 ; if (length(time.units)==0) time.units=""
  LON=ncdim_def("lon",lon.units,lon,create_dimvar=TRUE)
  LAT=ncdim_def("lat",lat.units,lat,create_dimvar=TRUE)
  TIME=ncdim_def("time",time.units,time,unlim=TRUE,create_dimvar=TRUE,calendar=caltype)
  varnc=ncvar_def(var,"",list(LON,LAT,TIME),missval=1.e+30)
  ncnew=nc_create(file,varnc)
  if (length(dim(dat))==2) dat=array(dat,dim=c(dim(dat),1))
  ncvar_put(ncnew,var,dat,start=NA,count=dim(dat))
  nc_close(ncnew); rm(ncnew); gc()
#  system(paste("ncatted -a _FillValue,",var,",c,f,1.e+30 ",file,sep=""))
}

# Open netcdf > list (dealing with fill value when missing value not specified)
myno=function(file,varid=NA){
  nc=nc_open(file)
#  if (system(paste("ncdump -h",file,"| grep 'missing_value ='"))==1
#      & system(paste("ncdump -h",file,"| grep '_FillValue ='"))!=1 ){
#    system(paste("ncdump -h",file,"| grep '_FillValue ='  > tmp"))
#    tmp=scan("tmp",what="character"); system("rm -f tmp")
#    varname=strsplit(tmp[1],split=":")[[1]][1]
#    fv=strsplit(tmp[3],split="")[[1]]
#    fillval=paste(fv[-length(fv)],collapse="")
#    nc$var[[paste(varname)]]$missval=as.numeric(fillval)
#  }
  out=list(varname=last(names(nc$var)))
  out$lon=nc$dim$lon$vals; out$lonu=nc$dim$lon$units
  out$lat=nc$dim$lat$vals; out$latu=nc$dim$lat$units
  if (any(names(nc$dim)=="plev")) {out$lev=nc$dim$plev$vals; out$levu=nc$dim$plev$units}
  out$time=floor(as.numeric(nc$dim$time$vals)); out$timeu=nc$dim$time$units
  if (any(names(nc$dim$time)=="calendar")) out$caltype=nc$dim$time$calendar
  out$var=ncvar_get(nc,varid=varid); nc_close(nc); rm(nc); gc()
  out
}

# Only getting dimensions of a netcdf (not the var)
mynd=function(file){
  nc=nc_open(file)
  out=list(varname=names(nc$var))
  out$lon=nc$dim$lon$vals; out$lonu=nc$dim$lon$units
  out$lat=nc$dim$lat$vals; out$latu=nc$dim$lat$units
  if (any(names(nc$dim)=="plev")) {out$lev=nc$dim$plev$vals; out$levu=nc$dim$plev$units}
  out$time=floor(as.numeric(nc$dim$time$vals)); out$timeu=nc$dim$time$units
  if (any(names(nc$dim$time)=="calendar")) out$caltype=nc$dim$time$calendar
  nc_close(nc); rm(nc); gc()
  out
}

# View netcdf with ncview
mynv=function(file) system(paste("ncview",file,"&"))

# Applying a cdo command to one (or more) ifile(s) -> ofile
mycdo=function(cdocmd,ifiles,ofile)
  system(paste("cdo",cdocmd,ifiles,ofile))

#----------------------------------------------------------------------
# OLD FUNCTIONS ---> TO BE CHECKED
#----------------------------------------------------------------------
# Function for classifying elements of dat wrt cls (eg for attributing to one cluster)
# (dat matrix nobs*size, cls matrix ncls*size)
myclassif=function(dat,cls,method=c("euclidean","correl","corrank"),mindist=NULL,mincorr=NULL){
  meth=method[1]
  if (meth=="euclidean"){
   fun=function(x,y) apply(x,1,myeuclid,y)
   best=function(x) c(sort(x,index.return=TRUE)$ix[1],min(x)/ncol(dat))
  }
  if (meth=="correl"){
   fun=function(x,y) apply(x,1,cor,y)
   best=function(x) c(sort(x,decreasing=TRUE,index.return=TRUE)$ix[1],max(x))
  }
  if (meth=="corrank"){
   fun=function(x,y) apply(x,1,cor,y,method="spearman")
   best=function(x) c(sort(x,decreasing=TRUE,index.return=TRUE)$ix[1],max(x))
  }
  classif=c(); for (i in 1:nrow(dat)){
    mes=fun(cls,dat[i,]); j=best(mes)
    if (length(mindist)>0 && myeuclid(dat[i,],cls[j[1],])>mindist) j=c(NA,NA)
    if (length(mincorr)>0 && cor(dat[i,],cls[j[1],])<mincorr) j=c(NA,NA)
    classif=c(classif,j)
  }
  classif=data.frame(matrix(classif,ncol=2,byrow=TRUE))
  names(classif)=c("class","dist")
  classif
}

# NetCDF (1 seul pas de temps) en tableau LON LAT VAR
netcdf2dat=function(file.in,file.out="OUT.dat"){
  nc=open.ncdf(file.in)
  x=nc$dim$lon$vals
  y=nc$dim$lat$vals
  z=get.var.ncdf(nc)
  repx=rep(x,each=length(y))
  repy=rep(y,length(x))
  out=data.frame(LON=repx,LAT=repy,VAR=c(t(z)))
  write.table(out,file.out,quote=FALSE,row.names=FALSE)
}

# Fonction pour calculer moyennes saisonnieres d'indices
indices.seas=function(index){
  a=as.matrix(read.table(paste("/home/e2c2/jcattiau/",index,".dat",sep="")))
  a[a==-99.99]=NA
  nyears=dim(a)[1]-1
  b=matrix(NA,nyears,12)
  b[,1:10]=a[1:nyears,4:13]
  b[,11:12]=a[2:(nyears+1),2:3]
  ave=data.frame(years=a[1:nyears,1])
  seasons=c("MAM","JJA","SON","DJF")
  for (i in 1:4) ave[,paste(seasons[i])]=apply(b[,(3*(i-1)+1):(3*i)],1,mean)
  write.table(ave,paste("/home/e2c2/jcattiau/",index,"_seas.dat",sep=""),row.names=FALSE,quote=FALSE)
}

# Fonction qui calcule la similarite entre 2 cartes (matrices lon*lat), a ponderer par la latitude (racine du cos)
maps.compare=function(map1,map2,lat,method=c("cor","pval","esv","eulidean")){
  method=method[1]
  coslat=rep(sqrt(cos(lat*pi/180)),each=dim(map1)[1])
  cm1=c(map1)*coslat ; cm2=c(map2)*coslat
  if (method=="cor") out=cor(cm1,cm2)
  if (method=="corrank") out=cor(c(map1),c(map2),method="spearman")
  if (method=="pval") out=cor.test(cm1,cm2)$p.value
  if (method=="esv") out=1-sd(cm1-cm2,na.rm=TRUE)^2/sd(cm1,na.rm=TRUE)^2
  if (method=="euclidean") out=sqrt(sum((cm1-cm2)^2))
  out
}

# Fonction qui interpole une carte map lon lat sur lonref latref
myinterpol=function(map,lon,lat,lonref,latref){
  print("Warning: be sure lon/lonref are -180..180 and not 0..360")
  out=matrix(NA,length(lonref),length(latref))
  for (i in 1:length(lonref)){
    # Lons les plus proches
    LON=lonref[i]
    lon1=min(lon); if (any(lon<=LON)) lon1=max(lon[which(lon<=LON)])
    i1=which(lon==lon1)
    lon2=max(lon); if (any(lon>=LON)) lon2=min(lon[which(lon>=LON)])
    i2=which(lon==lon2)
    for (j in 1:length(latref)){
      # Lats les plus proches
      LAT=latref[j]
      lat1=min(lat); if (any(lat<=LAT)) lat1=max(lat[which(lat<=LAT)])
      j1=which(lat==lat1)
      lat2=max(lat); if (any(lat>=LAT)) lat2=min(lat[which(lat>=LAT)])
      j2=which(lat==lat2)
      # Poids
      wn=0.5; if (lon1!=lon2) wn=(lon2-LON)/(lon2-lon1)
      wt=0.5; if (lat1!=lat2) wt=(lat2-LAT)/(lat2-lat1)
      # Moyenne ponderee
      out[i,j]=wn*wt*map[i1,j1]+(1-wn)*wt*map[i2,j1]
               +wn*(1-wt)*map[i1,j2]+(1-wn)*(1-wt)*map[i2,j2]
    }
  }
  out
}

#----------------------------------------------------------------------
# Graphical functions
#----------------------------------------------------------------------
# Reset des parametres graphiques
resetPar=function() {
  dev.new(); op=par(no.readonly = TRUE); dev.off(); par(op) }

# Plot vide
plot0=function(xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",axes=FALSE,frame.plot=FALSE,main="")
  plot(1,type="n",frame.plot=frame.plot,axes=axes,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)

# mtext ameliore pour noms axes
mymtext=function(text,side,line=0,at=NA,srt = if(side == 2) 90  else if(side == 4) 270 else 0,...) {
  # dimensions of plotting region in user units
  usr=par('usr'); pin=par('pin')
  # user units per inch
  upi=c(usr[2]-usr[1],usr[4]-usr[3])/pin
  # margins values (raw and inches) to derive line in inches
  mar=par('mar'); mai=par('mai'); lin=mai/mar
  # default x and y positions
  xpos=at; ypos=at ; if (is.na(at)){ xpos=(usr[1]+usr[2])/2; ypos=(usr[3]+usr[4])/2 }
  if(1==side) ypos=usr[3]-line*lin[1]*upi[2]
  if(2==side) xpos=usr[1]-line*lin[2]*upi[1]
  if(3==side) ypos=usr[4]+line*lin[3]*upi[2]
  if(4==side) xpos=usr[2]+line*lin[4]*upi[1]
  text(x=xpos,y=ypos,text,xpd=NA,srt=srt,adj=c(0.5,-0.3+1.6*(side==1)),...)
}

# Boxplot personnalises
mybp=function(x,at=1,col="transparent",range=0.00,quantiles=c(1/6,5/6),boxwex=1,pch=1,lty=1,add=F,horiz=F){
  x=na.omit(x); y=boxplot(x,range=range,plot=F,horizontal=horiz)
  qx=quantile(x,quantiles); y$stats[c(2,4),1]=matrix(qx); noout=which(y$out>=qx[1] & y$out<=qx[2])
  if (length(noout)>0) { y$out=y$out[-noout]; y$group=y$group[-noout] }
  bxp(y,at=at,boxfill=col,add=add,boxwex=boxwex,pch=pch,
      frame.plot=F,axes=F,xlab="",ylab="",boxlty=0,whisklty=lty,staplelty=1,horizontal=horiz)
}

# Plot pour ajouter un confidence interval
myci=function(x,at,width=1,col=1,lty=1,lwd=1,horiz=FALSE){
  if (!horiz){ segments(at-width/2,max(x),at+width/2,max(x),col=col,lty=lty,lwd=lwd)
               segments(at-width/2,min(x),at+width/2,min(x),col=col,lty=lty,lwd=lwd)
               segments(at,min(x),at,max(x),col=col,lty=lty,lwd=lwd)}
  if (horiz) { segments(max(x),at-width/2,max(x),at+width/2,col=col,lty=lty,lwd=lwd)
               segments(min(x),at-width/2,min(x),at+width/2,col=col,lty=lty,lwd=lwd)
               segments(min(x),at,max(x),at,col=col,lty=lty,lwd=lwd)}
}

# Plot type histos (polygones) autour d'un moyenne
mybarplot=function(y,at=NULL,y0=0,width=1,pch=NULL,pbg=7,col=c("red","blue"),border=NULL,fill=NULL,lty=1,lwd=1,plot.ci=FALSE,ci.u=NULL,ci.l=NULL,ci.lty=1,ci.lwd=1,xlim=NULL,yaxis=FALSE,ylab="",ylim=NULL,y0line=TRUE,col.levs=NULL,col.vals=y-y0,add=FALSE,frame.plot=FALSE){
  # Gestion arguments entree
  if (plot.ci & (length(ci.u)!=length(y) | length(ci.l)!=length(y)))
    {print("Error specifying CI"); stop()}
  if (length(at)==0) at=1:length(y)
  if (length(width)==1) width=rep(width,length(y))
  if (length(pch)==1) pch=rep(pch,length(y))
  cls=col; if (length(col)==1) cls=rep(col,length(y)); if (length(col)==2){
    cls=col[1.5-0.5*sign(y-y0)]
    if (col==c("red","blue") && length(col.levs)>0){
      lcl=length(col.levs); csq=seq(1,0,-1/lcl); myred=rgb(1,csq,csq); myblue=rgb(csq,csq,1)
      cls=c(); for (i in 1:length(y)){
        ic=which(sort(c(abs(col.vals[i]),col.levs))==abs(col.vals[i]))
        cls=c(cls,cbind(myred,myblue)[ic,1.5-0.5*sign(y[i]-y0)])
      }
    }
  }
#  if (length(col)==length(y)) cls=col
  # Plot
  if (is.null(xlim)) xlim=range(at)+c(-1,1)*mean(width)
  if (is.null(ylim)) ylim=range(c(y0,y,ci.u,ci.l))
  if (!add) plot0(xlim=xlim,ylim=ylim,ylab=ylab,frame.plot=frame.plot)
  for (i in 1:length(y)){
    xi=at[i]+c(-0.5,0.5,0.5,-0.5)*width[i]; yi=c(y0,y0,y[i],y[i])
    polygon(xi,yi,col=cls[i],border=border,density=fill[i],lty=lty,lwd=lwd)
    if (plot.ci){
      xi=at[i]+c(-0.25,0.25)*width[i]
      segments(at[i],ci.l[i],at[i],ci.u[i],lty=ci.lty,lwd=ci.lwd)
      segments(xi[1],ci.l[i],xi[2],ci.l[i],lty=ci.lty,lwd=ci.lwd)
      segments(xi[1],ci.u[i],xi[2],ci.u[i],lty=ci.lty,lwd=ci.lwd)
    }
    if (length(pch)>0) points(at[i],y[i],pch=pch[i],bg=pbg)
  }
  if (yaxis) axis(2)
  if (y0line) lines(xlim,rep(y0,2),lty=1)
}

# Plot type histos (polygones) similaires a mybarplot, 
# mais avec plusieurs contributions par barres
# (y matrice, nrow=nb de contrib et ncol=nb de barres)
mymultibarplot=function(y,at=NULL,width=1,col=rainbow(nrow(y)),border=NULL,lty=1,lwd=1,yaxis=FALSE,ylab="",ylim=NULL,y0line=TRUE,add=FALSE){
  # Gestion arguments entree
  if (length(at)==0) at=1:ncol(y)
  if (length(width)==1) width=rep(width,ncol(y))
  if (length(col)==1) cls=rep(col,ncol(y))
  # Plot
  xlim=range(at)+c(-1,1)*mean(width)
  if (length(ylim)==0) ylim=max(apply(abs(y),2,sum))*c(-1/2,1/2)
  if (!add) plot0(xlim=xlim,ylim=ylim,ylab=ylab)
  for (i in 1:ncol(y)){
    xi=at[i]+c(-0.5,0.5,0.5,-0.5)*width[i]; yt=0; yb=0
    for (j in 1:nrow(y)){
      if (y[j,i]>=0) {yij=yt+c(0,0,y[j,i],y[j,i]); yt=yt+y[j,i]}
      if (y[j,i]<0)  {yij=yb+c(0,0,y[j,i],y[j,i]); yb=yb+y[j,i]}
      polygon(xi,yij,col=col[j],border=border,lty=lty,lwd=lwd)
    } 
  }
  if (yaxis) axis(2)
  if (y0line) lines(xlim,rep(0,2),lty=1)
}

# Plot avec couleur entre-courbes
myfillplot=function(x,y,col=rep("black",2),fillcol=c("blue","red"),lty=c(1,1),lwd=c(1,1),at=1:length(x),xlim=range(at),ylim=range(c(x,y)),xlab="",ylab="",frame.plot=FALSE,add=FALSE){
  if (!add) plot0(xlim,ylim,xlab,ylab,frame.plot)
  d=y-x; s=sign(d)
  fcs=function(s) fillcol[ceiling(s/2+1)]
  if (length(x)>1) {for (i in 1:(length(x)-1)){
    if (s[i]==s[i+1])
      polygon(at[c(i,i+1,i+1,i)],c(x[i:(i+1)],y[(i+1):i]),col=fcs(s[i]),border=NA)
    else {
      q=abs(d[i])/sum(abs(d[i:(i+1)]))
      a=at[i]+q*(at[i+1]-at[i]); b=x[i]+q*(x[i+1]-x[i])
      polygon(c(at[i],a,a,at[i]),c(x[i],b,b,y[i]),col=fcs(s[i]),border=NA)
      polygon(c(a,at[i+1],at[i+1],a),c(b,x[i+1],y[i+1],b),col=fcs(s[i+1]),border=NA)
    }
  }}
  lines(at,x,col=col[1],lty=lty[1],lwd=lwd[1])
  lines(at,y,col=col[2],lty=lty[2],lwd=lwd[2])
}

# Scatterplot avec ronds de couleurs indiquant une 3e variable
myscatts=function(x,y,c=NULL,t=NULL,xlim=NULL,ylim=NULL,levs=NULL,cols=NULL,colorbar=T,regr=T,regr.txt=T,pval=0.1,regr.ci=T,regr.ci.level=0.95,add=F,frame.plot=T,pch=21,cex=2,bg.def="white"){
  if (length(levs)==0) {
    levs=seq(-1,1,0.2); levs=levs[-which(levs==0)]}
  nl=length(levs)
  if (length(cols)==0){
    l1=floor(nl/4); l2=nl/2-l1
    bl1=rgb(rep(0,l1),rep(0,l1),seq(0.5,0.8,le=l1))
    bl2=rgb(seq(0.1,0.9,le=l2),seq(0.1,0.9,le=l2),rep(1,l2))
    re2=rgb(rep(1,l2),seq(0.9,0.1,le=l2),seq(0.9,0.1,le=l2))
    re1=rgb(seq(0.8,0.5,le=l1),rep(0,l1),rep(0,l1))
    cols=c(bl1,bl2,"white",re2,re1)
  }
  if (length(which(!is.na(x)))==0 | length(which(!is.na(y)))==0){
    xlim=c(0,1); ylim=c(0,1)
    if (!add) plot0(xlim=xlim,ylim=ylim,frame.plot=T)
    text(sum(xlim)/2,sum(ylim)/2,"no data")
  }
  if (length(which(!is.na(x)))>0 & length(which(!is.na(y)))>0){
    if (length(xlim)==0) xlim=max(abs(x),na.rm=T)*c(-1.2,1.2)
    if (length(ylim)==0) ylim=max(abs(y),na.rm=T)*c(-1.2,1.2)
    if (!add) { plot0(xlim,ylim,frame.plot=frame.plot)
                abline(h=0,lty=3); abline(v=0,lty=3) }
    if (regr && cor.test(x,y)$p.value<pval) {
      lm=lm(y~x)
      if (regr.ci) {
        xci=data.frame(x=seq(min(x,na.rm=T),max(x,na.rm=T),le=10))
        xci=data.frame(x=seq(xlim[1],xlim[2],le=100))
        lci=predict(lm,xci,level=regr.ci.level,interval="confidence")
        lci=as.data.frame(lci)
        lines(xci$x,lci$lwr,lty=2); lines(xci$x,lci$upr,lty=2)
        polygon(c(xci$x,rev(xci$x)),c(lci$lwr,rev(lci$upr)),border=NA,col="gray90")
      }
      if (regr.txt) {
        rt=paste("r =",round(cor(na.omit(x),na.omit(y)),2))
        text(0.1*xlim[1]+0.9*xlim[2],0.05*ylim[1]+0.95*ylim[2],rt,cex=0.9)
      }
      abline(lm,lty=2)
    }
    ci=bg.def; if (length(c)>0) {ci=c(); for (j in 1:length(c))
      ci=c(ci,cols[which(sort(c(c[j],levs))==c[j])])}
    points(x,y,pch=pch,col=1,bg=ci,lwd=1.5,cex=cex)
    if (length(t)>0) text(x,y,t,cex=0.8,pos=3)
  }
  if (colorbar){
    xcb=seq(0.45*xlim[1]+0.55*xlim[2],xlim[2],le=nl+2)
    ycb=c(ylim[1],0.95*ylim[1]+0.05*ylim[2])
    for (i in 1:(nl+1))
      polygon(xcb[c(i,i+1,i+1,i)],ycb[c(1,1,2,2)],col=cols[i])
    text(xcb[-c(1,length(xcb))],ycb[2],levs,adj=c(0,0.2),cex=0.6,srt=60)
  }
}

# Color bar
mycolorbar=function(cols,levs,horiz=TRUE,srt=0,cex=2,width=1,main="",frame.plot=F){
  nl=length(levs); xlim=c(0,nl+1); ylim=c(0,1)
  if (horiz){
    xcb=seq(xlim[1],xlim[2],le=nl+2);ycb=ylim[2]-c(width,0)
    plot0(xlim=xlim,ylim=ylim,frame.plot=frame.plot)
    for (i in 1:(nl+1)) polygon(xcb[c(i,i+1,i+1,i)],ycb[c(1,1,2,2)],col=cols[i])
    text(xcb[-c(1,length(xcb))],ycb[1],levs,pos=1,cex=cex,srt=srt)
    text(mean(xlim),ylim[1],main,pos=3,cex=cex)
  }
  if (!horiz){
    xcb=seq(xlim[1],xlim[2],le=nl+2);ycb=ylim[1]+c(0,width)
    plot0(xlim=ylim,ylim=xlim,frame.plot=frame.plot)
    for (i in 1:(nl+1)) polygon(ycb[c(1,1,2,2)],xcb[c(i,i+1,i+1,i)],col=cols[i])
    text(ycb[2],xcb[-c(1,length(xcb))],levs,pos=4,cex=cex,srt=srt)
    text(ylim[2],mean(xlim),main,pos=2,cex=cex,srt=90)
  }
}

# Carte a partir d'une var en 3 colonnes
# (argument regular pour jeu gridded -> indiquer resolution)
XYZmap=function(df,clevs=(-5:5)[-6],ccols=c(topo.colors(20)[c(1,5,7,8)],cm.colors(20)[c(1)],"white",heat.colors(20)[c(18,13,9,5,1)]),regular=NULL,xlim=NULL,ylim=NULL,border=FALSE,colorbar=TRUE,lwd.grid=1,inverse.cols=FALSE,fill=FALSE,fill.col="gray80",add=FALSE){
  map=maps::map
  if (inverse.cols) ccols=ccols[length(ccols):1]
  if (!add){
    if (length(xlim)==0) xlim=range(df[,1])
    if (length(ylim)==0) ylim=range(df[,2])
    map(xlim=xlim,ylim=ylim,lwd=1.5,fill=fill,col=fill.col,mar=par("mar"))
  }
  for (i in 1:dim(df)[1]){ if (!is.na(df[i,3])){
    icol=which(sort(c(clevs,df[i,3]))==df[i,3])
    if (length(regular)>0){
      x1=df[i,1]-regular/2 ; x2=x1+regular ; y1=df[i,2]-regular/2 ; y2=y1+regular
      bd=NA ; if (border) bd=NULL
      polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),border=bd,col=ccols[icol])
    }
    else { 
      points(df[i,1],df[i,2],col=ccols[icol],pch=16)
      if (border) points(df[i,1],df[i,2])
    }
  }}
  if (!add){
    map.axes()
    if (length(lwd.grid)>0) map.grid(labels=FALSE,col="gray40",lwd=lwd.grid)
    if (length(regular)>0) map(lwd=1.5,add=TRUE,fill=fill,col=fill.col)
    if (colorbar) image.plot(as.matrix(0:length(ccols)),add=TRUE,legend.only=TRUE,legend.mar=8,legend.shrink=1,legend.width=0.8,col=ccols,breaks=0:length(ccols),lab.breaks=c("",clevs,""))
  }
}


