#---------------------------------------------------------------------------------
# Subroutine of compute.stat.1d() that computes the nD0 x nWIN matrix
# of event values from ts (vector)
#---------------------------------------------------------------------------------
compute.matrix.of.event.values=function(ts){
  # Initialize output
  dumY0=matrix(NA,nD0,nWIN)
  # Loop on time windows
  for (iw in 1:nWIN){ w=WINDOWS[iw]
    tsw=myrunavg(ts,w)
    # Loop on calendar days
    for (id in 1:nD0){ d=iDAYS0[id]
      dumY0[id,iw]=tsw[d]
    }
  }
  dumY0
}

#---------------------------------------------------------------------------------
# Subroutine of compute.stat.1d() that computes the nD0 x nWIN x nYRS array
# of data samples from ts (vector)
#---------------------------------------------------------------------------------
compute.array.of.data.samples=function(ts,annual.maxima=FALSE,calendar.flex=0,
                                       exceedance="above",first.year=TRUE,last.year=TRUE){
  # Initialize output
  dum=array(NA,dim=c(nD0,nWIN,nYRS))
  # Loop on time windows
  for (iw in 1:nWIN){ w=WINDOWS[iw]
    # Running mean of ts over w
    tsw=myrunavg(ts,w)
    # If annual.maxima, same sample for all calendar days (only id=1 is filled)
    if (annual.maxima){
      for (y in PERIOD) dum[1,iw,which(PERIOD==y)]=mymax(tsw[which(DATES$y==y)],exceedance)
    }
    # Else, loop on calendar days
    if (!annual.maxima){ for (id in 1:nD0){ d=iDAYS0[id]
      idays=which(DATES$m==DATES$m[d] & DATES$d==DATES$d[d])
      if (calendar.flex==0) dum[id,iw,]=tsw[idays]
      else {
        flex=-calendar.flex:calendar.flex
        idays=rep(idays,each=length(flex))+rep(flex,nYRS)
        idays[idays<1]=1 ; idays[idays>nrow(DATES)]=nrow(DATES)
        dum[id,iw,]=apply(matrix(tsw[idays],nrow=length(flex)),2,mymax,exceedance)
      }
    }}
    # Remove first and/or last years if requested (e.g. for annual.maxima when the year is incomplete)
    if (!first.year) dum[,,2]=NA
    if (!last.year)  dum[,,nYRS-1]=NA
  }
  dum
}

#---------------------------------------------------------------------------------
# Subroutine of compute.stat.1d() that computes the nD0 x nWIN x nYRS array
# of trend values from tr (vector of length ts or nYRS or matrix nWIN x nYRS)
#---------------------------------------------------------------------------------
compute.array.of.trend.values=function(ts,tr,annual.maxima=FALSE){
  # Initialize output
  dumtr=array(NA,dim=c(nD0,nWIN,nYRS))
  # If annual.maxima and daily trend, compute annual means of tr
  if (annual.maxima & length(ts)==length(tr)) tr=apply(matrix(tr,nrow=nD0),2,mean,na.rm=T)
  # If annual trend (no seasonal cycle), repeat annual value to compute dumtr
  if (length(tr)==nYRS){ for (i in 1:nYRS) dumtr[,,i]=tr[i] }
  # If annual trend with window dependence, repeat annual values for each window
  if (length(tr)==nWIN*nYRS) {for (i in 1:nWIN){ for (j in 1:nYRS) dumtr[,i,j]=tr[i,j]}}
  # If daily trend (seasonal cycle), compute running averages of tr to compute dumtr
  if (length(tr)==length(ts)){
    for (iw in 1:nWIN){ w=WINDOWS[iw]
      trw=myrunavg(tr,w)
      for (id in 1:nD0){ d=iDAYS0[id]
        idays=which(DATES$m==DATES$m[d] & DATES$d==DATES$d[d])
        dumtr[id,iw,]=trw[idays]
      }
    }
  }
  dumtr
}

#---------------------------------------------------------------------------------
# Function that returns a list of stats (p1, p0, pr) for a given location (station,
# lon-lat domain, country) and for all time windows requested in main.R (WINDOWS).
#   > ts and tr (vectors), or input with pre-computed dum, dumY0 and dumtr (arrays),
#     annual.maxima to fit annual max (default = calendar days), calendar.flex for
#     calendar local window, distrib (gauss, gev, etc.), fixed.ksi (EVT laws only),
#     smooth.fit.sdf to smooth fit parameters with mysmoothy(), and exceedance
#     to compute percentile level with myperclevel()
#   < stats in a data.frame object with colnames = WINDOWS
# Warning: the version provided below only includes a very simple detrending (the 
#          correction for climate change does not depend on the time window), which
#          can only be used with temperature variables (see C&R2018 for details)
#---------------------------------------------------------------------------------
compute.stat.1d=function(ts=NULL,tr=NULL,input=NULL,annual.maxima=FALSE,calendar.flex=0,
			 distrib="gauss",fixed.ksi=NULL,smooth.fit.sdf=NULL,
			 exceedance="above",first.year=TRUE,last.year=TRUE){
  #---- tr must be of size ts or nYRS
  if (! length(tr) %in% c(length(ts),nYRS)) return("Error: tr size must be ts, nYRS or nWINxnYRS.")
  #---- fixed.ksi must be of size 1 or nWIN
  if (! length(fixed.ksi) %in% c(0,1,nWIN)) return("Error: fixed.ksi must be size 1 or nWIN.")
  #---- Building arrays dum (samples), dumY0 (event values) and dumtr (trends)
  if (!is.null(input)) {dum=input$dum; dumY0=input$dumY0; dumtr=input$dumtr}
  else {
    dum=compute.array.of.data.samples(ts,annual.maxima,calendar.flex,exceedance,first.year,last.year)
    dumY0=compute.matrix.of.event.values(ts)
    dumtr=compute.array.of.trend.values(ts,tr,annual.maxima,first.year,last.year)
  }
  #---- Computing arrays for detrended data wrt. t=Y0 (dum1, for p1) and t=YREF (dum0, for p0)
  dum1=dum-dumtr+rep.abind(dumtr[,,iY0],nYRS)
  dum0=dum-dumtr+rep.abind(dumtr[,,iYREF],nYRS)
  #---- Fitting parameters of the requested distribution
  if (distrib %in% c("gev","gevmin") & !is.null(fixed.ksi)){
    fit1=dum1[,,1:3]; fit0=dum0[,,1:3]
    ksi=rep(fixed.ksi,length(WINDOWS))[1:length(WINDOWS)]
    for (iw in 1:nWIN){ for (id in 1:nD0){ if (any(!is.na(dum[id,iw,]))){
      fit1[id,iw,]=myfitparams(dum1[id,iw,],paste0(distrib,"fix"),ksi=ksi[iw])
      fit0[id,iw,]=myfitparams(dum0[id,iw,],paste0(distrib,"fix"),ksi=ksi[iw])
    }}}
  }
  else {
    fit1=aperm(apply(dum1,1:2,myfitparams,distrib),c(2,3,1))
    fit0=aperm(apply(dum0,1:2,myfitparams,distrib),c(2,3,1))
  }
  #---- Smoothing fits if requested
  if (!is.null(smooth.fit.sdf)){
    # Calendar: smoothing over days
    if (!annual.maxima){
      fit1=apply(fit1,2:3,mysmoothy,sdf=smooth.fit.sdf,periodic=(nD0==365))
      fit0=apply(fit0,2:3,mysmoothy,sdf=smooth.fit.sdf,periodic=(nD0==365))
    }
    # Annual max: smoothing over windows
    # EDIT for GEV: take ksi as the mean over windows and re-fit mu and sigma with ksi constant
    if (annual.maxima){
      if (! distrib %in% c("gev","gevmin")){
        fit1[1,,]=apply(fit1[1,,],2,mysmoothy,sdf=smooth.fit.sdf)
        fit0[1,,]=apply(fit0[1,,],2,mysmoothy,sdf=smooth.fit.sdf)
      }
      if (distrib %in% c("gev","gevmin")){
        ksi=mean(fit1[1,,3],na.rm=T)
        for (iw in 1:nWIN){ for (id in 1:nD0){ if (any(!is.na(dum[id,iw,]))){
          fit1[id,iw,]=myfitparams(dum1[id,iw,],paste0(distrib,"fix"),ksi)
          fit0[id,iw,]=myfitparams(dum0[id,iw,],paste0(distrib,"fix"),ksi)
        }}}
      }
    }
  }
  #---- If annual.maxima, fits are repeated over all days
  if (annual.maxima){ for (id in 2:nD0){ fit1[id,,]=fit1[1,,] ; fit0[id,,]=fit0[1,,] }}  
  #---- Computing p1 (percentile level of dumY0 within fit1)
  out.p1=abind(dumY0,fit1,along=3)
  out.p1=apply(out.p1,1:2,function(x) myperclevel(x[1],distrib,x[-1],exceedance))
  out.p1=as.data.frame(out.p1); names(out.p1)=WINDOWS
  #---- Computing p0 (percentile level of dumY0 within fit0)
  out.p0=abind(dumY0,fit0,along=3)
  out.p0=apply(out.p0,1:2,function(x) myperclevel(x[1],distrib,x[-1],exceedance))
  out.p0=as.data.frame(out.p0); names(out.p0)=WINDOWS
  #---- Computing probability ratio (warning: p0 and p1 are percentile levels here!)
  out.pr=(1-out.p1)/(1-out.p0)
  #---- Output
  list(p1=out.p1,p0=out.p0,pr=out.pr,
       data=dum,value=dumY0,trend=dumtr,#data1=dum1,data0=dum0,
       fit1=fit1,fit0=fit0,distrib=distrib)
}
