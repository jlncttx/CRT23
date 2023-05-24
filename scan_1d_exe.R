#-------------------------------------------------------------------------------------
# Computes scan for LOC / SCAN.YRS / METHOD
# Returns list 'SCAN' with results of compute.stat.1d() for all years + metadata
#-------------------------------------------------------------------------------------
VERBOSE=paste(VAR,LOC,paste(range(SCAN.YRS),collapse=" "),paste(METHOD,collapse=" "))

#-------------------------------------------------------------------------------------
# Arguments of compute.stat.1d() corresponding to METHOD
#-------------------------------------------------------------------------------------
if (METHOD=="calend"){
  annual.maxima=F ; calendar.flex=0     ; distrib="gauss"
  fixed.ksi=NULL  ; smooth.fit.sdf=30   ; exceedance="above"
  first.year=T    ; last.year=T
}
if (METHOD=="calend2"){
  annual.maxima=F ; calendar.flex=0     ; distrib="snorm"
  fixed.ksi=NULL  ; smooth.fit.sdf=30   ; exceedance="above"
  first.year=T    ; last.year=T
}
if (METHOD=="locmax"){
  annual.maxima=F ; calendar.flex=7     ; distrib="gauss"
  fixed.ksi=NULL  ; smooth.fit.sdf=30   ; exceedance="above"
  first.year=T    ; last.year=T
}
if (METHOD=="locmin"){
  annual.maxima=F ; calendar.flex=7     ; distrib="gauss"
  fixed.ksi=NULL  ; smooth.fit.sdf=30   ; exceedance="below"
  first.year=T    ; last.year=T
}
if (METHOD=="annmax"){ # GEV
  annual.maxima=T ; calendar.flex=0     ; distrib="gev"
  fixed.ksi=NULL  ; smooth.fit.sdf=NULL ; exceedance="above"
  first.year=T    ; last.year=LAST.YEAR.MAX
}
if (METHOD=="annmax1"){ # GEV avec ksi fixe a la moyenne des ksi
  annual.maxima=T ; calendar.flex=0     ; distrib="gev"
  fixed.ksi=NULL  ; smooth.fit.sdf=1    ; exceedance="above"
  first.year=T    ; last.year=LAST.YEAR.MAX
}
if (METHOD=="annmax2"){ # Gumbel
  annual.maxima=T ; calendar.flex=0     ; distrib="gev"
  fixed.ksi=0     ; smooth.fit.sdf=NULL ; exceedance="above"
  first.year=T    ; last.year=LAST.YEAR.MAX
}
if (METHOD=="annmax3"){ # Gaussien
  annual.maxima=T ; calendar.flex=0     ; distrib="gauss"
  fixed.ksi=NULL  ; smooth.fit.sdf=NULL ; exceedance="above"
  first.year=T    ; last.year=LAST.YEAR.MAX
}
if (METHOD=="annmin"){
  annual.maxima=T ; calendar.flex=0     ; distrib="gevmin"
  fixed.ksi=NULL  ; smooth.fit.sdf=NULL ; exceedance="below"
  first.year=T    ; last.year=LAST.YEAR.MIN
}
if (METHOD=="annmin1"){
  annual.maxima=T ; calendar.flex=0     ; distrib="gevmin"
  fixed.ksi=NULL  ; smooth.fit.sdf=1    ; exceedance="below"
  first.year=T    ; last.year=LAST.YEAR.MIN
}
if (METHOD=="annmin2"){
  annual.maxima=T ; calendar.flex=0     ; distrib="gevmin"
  fixed.ksi=0     ; smooth.fit.sdf=NULL ; exceedance="below"
  first.year=T    ; last.year=LAST.YEAR.MIN
}
if (METHOD=="annmin3"){
  annual.maxima=T ; calendar.flex=0     ; distrib="gauss"
  fixed.ksi=NULL  ; smooth.fit.sdf=NULL ; exceedance="below"
  first.year=T    ; last.year=LAST.YEAR.MIN
}

params=list(annual.maxima=annual.maxima,calendar.flex=calendar.flex,
            distrib=distrib,fixed.ksi=fixed.ksi,smooth.fit.sdf=smooth.fit.sdf,
            exceedance=exceedance,first.year=first.year,last.year=last.year)

#-------------------------------------------------------------------------------------
# Data time series ts (length 365 x nYRS)
#-------------------------------------------------------------------------------------
print(paste(VERBOSE,": preparing data"))
print(".. time series")
ts=XLOC$var

#-------------------------------------------------------------------------------------
# Forced response fr (length nYRS)
#-------------------------------------------------------------------------------------
print(".. forced response")
fr=TREND

#-------------------------------------------------------------------------------------
# Trend tr (size depends on METHOD)
#-------------------------------------------------------------------------------------
print(".. trend")

# Trend time series tr with seasonal cycle for calendar method (length 365 x nYRS)
if (METHOD %in% c("calend","calend2","locmax","locmin")) tr=fit_onto_forced_response(ts,fr,12,6)

# Trend time series tr for annmax/annmin methods, with a dependance to the time
# window, i.e. TX1d vs. TX10d have different trends (matrix nWIN x nYRS)
if (METHOD %in% c("annmax","annmax1","annmax2","annmax3")){
  # - compute annmax/annmin samples
  tsx=compute.array.of.data.samples(ts,ann=T,exc="above",first.year=first.year,last.year=last.year)[1,,]
  #  - fit annmax/annmin values onto forced response
  cfr=fr-mean(fr)
  f_tsx=t(apply(tsx,1,function(x) lm(x~cfr)$coef))
  #  - smooth scaling factors over time windows
  sf_tsx=mysmoothy(f_tsx[,2],sdf=2)
  #  - computes tr matrices
  tr=matrix(rep(f_tsx[,1],nYRS)+rep(sf_tsx,nYRS)*rep(cfr,each=nWIN),nWIN,nYRS)
}
if (METHOD %in% c("annmin","annmin1","annmin2","annmin3")){
  tsn=compute.array.of.data.samples(ts,ann=T,exc="below",first.year=first.year,last.year=last.year)[1,,]
  cfr=fr-mean(fr)
  f_tsn=t(apply(tsn,1,function(x) lm(x~cfr)$coef))
  sf_tsn=mysmoothy(f_tsn[,2],sdf=2)
  tr=matrix(rep(f_tsn[,1],nYRS)+rep(sf_tsn,nYRS)*rep(cfr,each=nWIN),nWIN,nYRS)
}

#-------------------------------------------------------------------------------------
# Initialize scan
#-------------------------------------------------------------------------------------
print(paste(VERBOSE,": initializing scan"))
Y0=PERIOD[1] ;  iY0=which(PERIOD==Y0) ; iDAYS0=which(DATES$y==Y0)

#-------------------------------------------------------------------------------------
# Data values (dum) samples (dumsp) and trends (dumtr) 
#-------------------------------------------------------------------------------------
print(".. computes compute.stat.1d() matrix inputs")

# Values (always the same)
dum=compute.array.of.data.samples(ts)

# Sample (depends on METHOD through parameters)
{
if (!annual.maxima & calendar.flex==0) dumsp=dum
else dumsp=compute.array.of.data.samples(ts,annual.maxima,calendar.flex,exceedance,first.year,last.year)
}

# Trend (depends on METHOD through tr computation above)
dumtr=compute.array.of.trend.values(ts,tr)

#-------------------------------------------------------------------------------------
# Scan
#-------------------------------------------------------------------------------------
print(paste(VERBOSE,": starting scan"))

# Initialize SCAN
SCAN=list(var=VAR,location=LOC,years=YEARS,period=PERIOD,scan.years=SCAN.YRS,yref=YREF,windows=WINDOWS,
          method=METHOD,params=params,data=dumsp,trend=dumtr,cmip.period=CMIP.PERIOD,cmip.mmm=XMMM,
          forced.response=TREND,scan=list())

# Loop on requested years
for (Y0 in SCAN.YRS){
  print(paste(LOC,": scan",METHOD,Y0))
  # Indices of year/days of the considered event (hard coded..)
  iY0=which(PERIOD==Y0) ; iDAYS0=which(DATES$y==Y0)
  # Matrix of event values
  dumY0=dum[,,iY0]
  # Compute stat 1d
  out=compute.stat.1d(input=list(dum=dumsp,dumY0=dumY0,dumtr=dumtr),
		      ann=annual.maxima,cal=calendar.flex,dis=distrib,
		      fix=fixed.ksi,smo=smooth.fit.sdf,exc=exceedance,
		      fir=first.year,las=last.year)
  # Save in SCAN
  SCAN$scan[[paste0("Y",Y0)]]=list(p1=out$p1,p0=out$p0,value=out$value,fit1=out$fit1,fit0=out$fit0)
}

