#---------------------------------------------------------------------------------
# Load rbase.R (packages and base functions) + a few more packages
#---------------------------------------------------------------------------------
source("rbase.R")
library(extRemes)   # GEV
library(evd)        # GEV with fixed params
#library(ismev)      # other GEV package
#library(fGarch)     # Skew-normal
library(geosphere)  # geometry on Earth (e.g. distance)
library(data.table) # reads files faster than read.table()

#---------------------------------------------------------------------------------
# Directories of data common to all analyses (scan computation + visualization)
#---------------------------------------------------------------------------------
CMIPDIR="../data/cmip/"
WRAFDIR="../data/wraf/"
SCANDIR="./"

#---------------------------------------------------------------------------------
# Dealing with CMIP data to estimate the forced response
#---------------------------------------------------------------------------------
# HIST+SSP CMIP6 multi-model mean for estimating the forced response
cmipfile=paste0(CMIPDIR,"CMIP6_tas_hist585_yr_raw_ensmean.nc")
dum=mynd(cmipfile)
CMIP.PERIOD=time2mdy(dum$time,dum$timeu)$y

# Get tools and params to estimate the forced response in multi-model mean
source("forced_response.R")
RF=read.table(paste0(CMIPDIR,"CMIP6_Forcings.dat"),header=T)
EBMs=read.table(paste0(CMIPDIR,"CMIP5_ebm_params.dat"),header=T)
Enat=ebm_response(RF,EBMs,CMIP.PERIOD)$be

#---------------------------------------------------------------------------------
# Dealing with WRAF region-averaged data
#---------------------------------------------------------------------------------
source("wraf_source.R")

#---------------------------------------------------------------------------------
# A few useful functions
#---------------------------------------------------------------------------------
# Running mean accounting for NAs and centering
myrunavg=function(x,n) stats::filter(x,rep(1/n,n))

# Maximum accounting for NAs (if exceedance=below, computes minimum)
mymax=function(x,exceedance="above"){
  out=NA; if (any(!is.na(x))){
    if (exceedance=="above") out=max(x,na.rm=T)
    if (exceedance=="below") out=min(x,na.rm=T)
  }
  out
}

# Various fits of parametric laws, output is systematically 3 parameters (possibly NA)
myfitparams=function(dum,distrib,ksi=NULL){
  #dum[iY0]=NA ### to exclude the value of the year of interest (we don't want to do that)
  dum=na.omit(dum); out=rep(NA,3); if (length(dum)>0){
    if (distrib=="gauss")     out=c(mean(dum),sd(dum),NA)
    #  if (distrib=="gevlmo")    out=unname(fevd(dum,method="Lmoments")$results)
    #  if (distrib=="gevmle")    out=unname(fevd(dum,method="MLE")$results$par)
    #  if (distrib=="gevmle")    out=unname(gev.fit(dum,show=F)$mle)
    if (distrib=="gev"){
      out=unname(fgev(dum,std.err=F,method="CG")$param)
      if (out[2]<0) out=unname(fgev(dum,std.err=F,method="L-BFGS-B")$param)
    }
    if (distrib=="gevfix")    out=unname(fgev(dum,shape=ksi,std.err=F)$param)
    if (distrib=="gumbel")    out=c(unname(fevd(dum,type="Gumbel")$results$par),0)
    #  if (distrib=="gevminlmo") out=unname(fevd(-dum,method="Lmoments")$results)
    #  if (distrib=="gevminmle") out=unname(fevd(-dum,method="MLE")$results$par)
    #  if (distrib=="gevminmle") out=unname(gev.fit(-dum,show=F)$mle)
    if (distrib=="gevmin"){
      out=unname(fgev(-dum,std.err=F,method="CG")$param)
      if (out[2]<0) out=unname(fgev(-dum,std.err=F,method="L-BFGS-B")$param)
    }
    if (distrib=="gevminfix") out=unname(fgev(-dum,shape=ksi,std.err=F)$param)
    if (distrib=="gumbelmin") out=c(unname(fevd(-dum,type="Gumbel")$results$par),0)
    if (distrib=="snorm")     out=unname(snormFit(dum)$par)
  }
  if (!is.na(out[2]) & out[2]<=0){
    print("Warning: negative scale parameter found --> force parameters to NAs.")
    out=rep(NA,3)
  }
  out
}

# PDF function for various parametric laws
mypdffun=function(distrib,params){
  out=NULL; if (any(!is.na(params))){
    if (distrib=="gauss")
      out=function(x) dnorm(x,params[1],params[2])
    if (distrib %in% c("gev","gevmle","gevfix","gumbel"))
      out=function(x) devd(x,params[1],params[2],params[3])
    if (distrib %in% c("gevmin","gevminmle","gevminfix","gumbelmin"))
      out=function(x) devd(-x,params[1],params[2],params[3])
    if (distrib=="snorm")
      out=function(x) dsnorm(x,params[1],params[2],params[3])
  }
  out
}

# Percentile level for various parametric laws (exceedance above or below)
myperclevel=function(dum0,distrib,params,exceedance="above"){
  out=NA; if (!is.na(dum0) & any(!is.na(params))){
    if (distrib=="gauss")
      out=pnorm(dum0,params[1],params[2])
    if (distrib %in% c("gev","gevmle","gevfix","gumbel"))
      out=pevd(dum0,params[1],params[2],params[3])
    if (distrib %in% c("gevmin","gevminmle","gevminfix","gumbelmin"))
      out=1-pevd(-dum0,params[1],params[2],params[3])
    if (distrib=="snorm")
      out=psnorm(dum0,params[1],params[2],params[3])
  }
  if (exceedance=="above") return(out)
  if (exceedance=="below") return(1-out)
}

# KS-test for a given sample and a given distribution
mykstest=function(dum,distrib,params){
  out=NA; if (any(!is.na(dum)) & any(!is.na(params))){
    if (distrib=="gauss")
      out=ks.test(dum,pnorm,params[1],params[2])$p.value
    if (distrib %in% c("gev","gevmle","gevfix","gumbel"))
      out=ks.test(dum,pevd,params[1],params[2],params[3])$p.value
    if (distrib %in% c("gevmin","gevminmle","gevminfix","gumbelmin"))
      out=ks.test(-dum,pevd,params[1],params[2],params[3])$p.value
    if (distrib=="snorm")
      out=ks.test(dum,psnorm,params[1],params[2],params[3])$p.value
  }
  out  
}
