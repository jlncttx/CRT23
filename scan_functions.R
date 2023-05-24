#---------------------------------------------------------------------------------
# Special functions for multi-year 'scan' analysis
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Get a 'scan_1d_exe' object
#---------------------------------------------------------------------------------
get.scan.1d=function(var="TM",meth="calend",loc="France",set="mf"){
  out=NULL
  rfile=paste0(SCANDIR,"scans_1d/",set,"/scan_",loc,"_",var,"_",meth,".Rdata")
  if (file.exists(rfile)){
    load(rfile)
    SCAN$set=set
    out=SCAN[c(length(SCAN),1:(length(SCAN)-1))]
  }
  out
}

#---------------------------------------------------------------------------------
# Get a 'scan_2d_exe' object
#---------------------------------------------------------------------------------
get.scan.2d=function(var="TM",meth="annmax",set="era5",wset=2,years=NULL,only.max=TRUE){
  out=NULL
  rfile=paste0(SCANDIR,"scans_2d/",set,"/events_WRAF",wset,"_",var,"_",meth,".Rdata")
  if (file.exists(rfile)){
    load(rfile)
    if (only.max){
      out=YEARLY.EVENTS
      out$selparams$years=years
      if (!is.null(years)) out$events=out$events[,,which(out$scan.years %in% years)]
      if (length(dim(out$events))==3){
        dum=out$events[,,1]; for (i in 1:nrow(dum)){
          imax=which.max(out$events[i,which(colnames(out$events)==out$selparams$stat.max),])
          dum[i,]=out$events[i,,imax]
        }
        out$events=dum
      }
      dum=matrix(c(out$events),nrow=nrow(out$events))
      mud=type.convert(dum[,-ncol(dum)],dec=".")
      mud=as.data.frame(mud); names(mud)=colnames(out$events)[-ncol(dum)]
      out$events=mud
      out$events$loc=noquote(dum[,ncol(dum)])
    }
    if (!only.max){
      out=MAX.EVENTS
      if (!is.null(years)) out$events=out$events[which(out$events$y %in% years),]
    }
  }
  out
}

#---------------------------------------------------------------------------------
# Identifies the most extreme events within a 'scan_1d_exe' object
# Criteria: local max + non-overlap, + thresholds in arguments
#---------------------------------------------------------------------------------
get.max.events=function(s=NULL,var="TM",meth="annmax",loc="France",set="mf",
                        stat.max="p1",years=NULL,only.max=FALSE,only.months=NULL,
                        min.value=NULL,min.win=NULL,max.win=NULL,reverse.prob=F){
  # Data
  if (is.null(s)) s=get.scan.1d(var,meth,loc,set)
  var=s$var; meth=s$method; loc=s$location; set=s$set
  # Full maps of p1 and p0
  if (!is.null(years)) yrs=years else yrs=s$scan.years
  dum=list(); for (n in c("p1","p0")){
    dums=c(); for (y in yrs) dums=rbind(dums,s$scan[[paste0("Y",y)]][[paste(n)]])
    dum[[paste(n)]]=dums
  }
  # Events based on p or 1-p?
  if (reverse.prob){ dum$p1=1-dum$p1 ; dum$p0=1-dum$p0 }
  # Stat to maximize 
  dum0=dum[[paste(stat.max)]]
  # Restriction to values / windows / months
  if (!is.null(min.value) & !only.max) dum0[dum0<min.value]=NA
  if (!is.null(min.win)) dum0[,which(s$windows<min.win)]=NA
  if (!is.null(max.win)) dum0[,which(s$windows>max.win)]=NA
  if (!is.null(only.months)) {
    dts=makedates(y=yrs,avg='day',caltype='noleap')
    imths=which(! dts$m %in% only.months)
    dum0[imths,]=NA
  }
  # If only max, selecting global maximum (patch for ex-aequo, eg. p1=1)
  if (only.max){
    imax=which(dum0==max(dum0,na.rm=T),arr.ind=T)
    if (nrow(imax)>1){
      print(paste("Warning: ex-aequo with",stat.max,"=",max(dum0,na.rm=T),"!"))
      imax=t(as.matrix(imax[nrow(imax),]))
    }
  }
  # Else
  if (!only.max){
    # Detecting local maxima
    dum1=dum0; for (i in 1:nrow(dum1)){ for (j in 1:ncol(dum1)){ if (!is.na(dum0[i,j])){
      is=max(1,i-1):min(nrow(dum1),i+1)
      js=max(1,j-1):min(ncol(dum1),j+1)
      if (max(dum0[is,js],na.rm=T)!=dum0[i,j]) dum1[i,j]=NA
    }}}
    # Selecting local maxima in descending order and masking overlaps
    imax=c(); dumi=dum1; while (sum(!is.na(dumi))>0){
      imaxi=which(dumi==max(dumi,na.rm=T),arr.ind=T)
      # Patch for ex aequo (eg. p1=1)
      if (nrow(imaxi)>1){
        print(paste("Warning: ex-aequo with",stat.max,"=",max(dumi,na.rm=T),"!"))
        imaxi=t(as.matrix(imaxi[nrow(imaxi),]))
      }
      # Increment imax
      imax=rbind(imax,imaxi)
      # Mask ovelapping events
      win=s$windows[imaxi[2]]
      day1=imaxi[1]-trunc((win-1)/2); day2=day1+win-1
      for (i in 1:length(s$windows))
        dumi[max(c(1,day1-s$windows[i]+1)):min(c(nrow(dumi),day2+s$windows[i]-1)),i]=NA
    }
  }
  # Generating metadata from imax
  out=c(); if (length(imax)>0){
    for (i in 1:nrow(imax)){
      win=s$windows[imax[i,2]]
      day1=imax[i,1]-trunc((win-1)/2); day2=day1+win-1
      date1=mdy2yyyymmdd(time2mdy(day1,paste0("days since ",yrs[1]-1,"-12-31"),"noleap"))
      date2=mdy2yyyymmdd(time2mdy(day2,paste0("days since ",yrs[1]-1,"-12-31"),"noleap"))
      ymd=time2mdy(imax[i,1],paste0("days since ",yrs[1]-1,"-12-31"),"noleap")
      outi=c(dum$p1[imax[i,1],imax[i,2]],dum$p0[imax[i,1],imax[i,2]],day1,day2,date1,date2,win,ymd$y,ymd$m,ymd$d,imax[i,])
      out=rbind(out,outi)
    }
    out=as.data.frame(out,row.names="")
    names(out)=c("p1","p0","day1","day2","date1","date2","win","y","m","d","id","iw")
  }
  # Returns a list with metadata
  output=list()
  for (n in c("set","var","location","years","period","scan.years","yref","windows","method","params"))
    output[[paste(n)]]=s[[paste(n)]]
  output$selparams=list(stat.max=stat.max,years=years,only.max=only.max,min.value=min.value,
                        min.win=min.win,max.win=max.win,reverse.prob=reverse.prob)
  output$events=out
  output
}

#---------------------------------------------------------------------------------
# Identifies all most extreme events over a set of wraf regions
#---------------------------------------------------------------------------------
get.all.max.events=function(var="TX",meth="annmax",set="era5",wset=2,stat.max="p1",years=NULL,
                            only.max=TRUE,only.months=NULL,min.value=NULL,min.win=NULL,max.win=NULL,reverse.prob=F){
  # Get all scans and max events
  out=c()
  for (loc in REGIONS[[paste0("WRAF",wset)]][1,]){
    s=get.scan.1d(var,meth,loc,set)
    dum=get.max.events(s,stat.max=stat.max,years=years,only.max=only.max,only.months=only.months,
                       min.value=min.value,max.win=max.win,min.win=min.win,reverse.prob=reverse.prob)$events
    dum$loc=loc
    out=rbind(out,dum)
  }
  # Returns a list with metadata
  output=list()
  for (n in c("set","var","years","period","scan.years","yref","windows","method","params"))
    output[[paste(n)]]=s[[paste(n)]]
  output$wraf.set=wset
  output$selparams=list(stat.max=stat.max,years=years,only.max=only.max,only.months=only.months,
                        min.value=min.value,min.win=min.win,max.win=max.win,reverse.prob=reverse.prob)
  output$events=out
  output
}

#---------------------------------------------------------------------------------
# Identifies the most extreme events within multiple 'scan_1d_exe' objects (a family)
# See get.max.events for details
#---------------------------------------------------------------------------------
get.family.max.events=function(s=NULL,var="TM",meth="annmax",family=wraf.family("X.2"),set="era5",
                        stat.max="p1",years=NULL,only.max=FALSE,only.months=NULL,
                        min.value=NULL,min.win=NULL,max.win=NULL,reverse.prob=F){
  # Data
  if (is.null(s)){
    s=list(); for (loc in family) s[[paste(loc)]]=get.scan.1d(var,meth,loc,set)
  }
  var=s[[1]]$var; meth=s[[1]]$method; family=names(s); set=s[[1]]$set
  # Full maps of p1 and p0 (arrays with dim 3 for family)
  if (!is.null(years)) yrs=years else yrs=s[[1]]$scan.years
  dum=list(); for (n in c("p1","p0")){
    dumf=c(); for (loc in family){
      dums=c(); for (y in yrs) dums=rbind(dums,s[[paste(loc)]]$scan[[paste0("Y",y)]][[paste(n)]])
      dumf=abind(dumf,dums,along=3)
    }
    dum[[paste(n)]]=dumf
  }
  # Events based on p or 1-p?
  if (reverse.prob){ dum$p1=1-dum$p1 ; dum$p0=1-dum$p0 }
  # Stat to maximize 
  dum0=dum[[paste(stat.max)]]
  # Restriction to values / windows / months
  if (!is.null(min.value) & !only.max) dum0[dum0<min.value]=NA
  if (!is.null(min.win)) dum0[,which(s[[1]]$windows<min.win),]=NA
  if (!is.null(max.win)) dum0[,which(s[[1]]$windows>max.win),]=NA
  if (!is.null(only.months)) {
    dts=makedates(y=yrs,avg='day',caltype='noleap')
    imths=which(! dts$m %in% only.months)
    dum0[imths,,]=NA
  }
  # If only max, selecting global maximum (patch for ex-aequo, eg. p1=1)
  if (only.max){
    imax=which(dum0==max(dum0,na.rm=T),arr.ind=T)
    if (nrow(imax)>1){
      print(paste("Warning: ex-aequo with",stat.max,"=",max(dum0,na.rm=T),"!"))
      imax=t(as.matrix(imax[nrow(imax),]))
    }
  }
  # Else
  if (!only.max){
    # Detecting local maxima (for each loc)
    dum1=dum0; for (k in 1:dim(dum1)[3]){
      for (i in 1:dim(dum1)[1]){ for (j in 1:dim(dum1)[2]){ if (!is.na(dum0[i,j,k])){
        is=max(1,i-1):min(dim(dum1)[1],i+1)
        js=max(1,j-1):min(dim(dum1)[2],j+1)
        if (max(dum0[is,js,k],na.rm=T)!=dum0[i,j,k]) dum1[i,j,k]=NA
      }}}
    }
    # Selecting local maxima in descending order and masking overlaps
    imax=c(); dumi=dum1; while (sum(!is.na(dumi))>0){
      imaxi=which(dumi==max(dumi,na.rm=T),arr.ind=T)
      # Patch for ex aequo (eg. p1=1)
      if (nrow(imaxi)>1){
        print(paste("Warning: ex-aequo with",stat.max,"=",max(dumi,na.rm=T),"!"))
        imaxi=t(as.matrix(imaxi[nrow(imaxi),]))
      }
      # Increment imax
      imax=rbind(imax,imaxi)
      # Mask ovelapping events
      win=s[[1]]$windows[imaxi[2]]
      day1=imaxi[1]-trunc((win-1)/2); day2=day1+win-1
      for (i in 1:length(s[[1]]$windows))
        dumi[max(c(1,day1-s[[1]]$windows[i]+1)):min(c(nrow(dumi),day2+s[[1]]$windows[i]-1)),i,]=NA
    }
  }
  # Generating metadata from imax
  out=c(); if (length(imax)>0){
    for (i in 1:nrow(imax)){
      win=s[[1]]$windows[imax[i,2]]
      day1=imax[i,1]-trunc((win-1)/2); day2=day1+win-1
      date1=mdy2yyyymmdd(time2mdy(day1,paste0("days since ",yrs[1]-1,"-12-31"),"noleap"))
      date2=mdy2yyyymmdd(time2mdy(day2,paste0("days since ",yrs[1]-1,"-12-31"),"noleap"))
      ymd=time2mdy(imax[i,1],paste0("days since ",yrs[1]-1,"-12-31"),"noleap")
      loc=names(s)[imax[i,3]]
      outi=c(dum$p1[imax[i,1],imax[i,2],imax[i,3]],dum$p0[imax[i,1],imax[i,2],imax[i,3]],
             day1,day2,date1,date2,win,ymd$y,ymd$m,ymd$d,imax[i,1:2],loc)
      out=rbind(out,outi)
    }
    out=as.data.frame(out,row.names="")
    names(out)=c("p1","p0","day1","day2","date1","date2","win","y","m","d","id","iw","loc")
  }
  # Returns a list with metadata
  output=list()
  for (n in c("set","var","years","period","scan.years","yref","windows","method","params"))
    output[[paste(n)]]=s[[paste(n)]]
  output$family=family
  output$selparams=list(stat.max=stat.max,years=years,only.max=only.max,min.value=min.value,
                        min.win=min.win,max.win=max.win,reverse.prob=reverse.prob)
  output$events=out
  output
}

#---------------------------------------------------------------------------------
# Gets calendar event for a given date: same format as max.eve but only 1 event
# Warning: s must be a calendar scan!
#---------------------------------------------------------------------------------
get.cal.event=function(s=NULL,var="TM",loc="France",set="mf",y=2003,id=220,iw=9,reverse.prob=F){
  # Data
  if (is.null(s)) s=get.scan.1d(var,"calend",loc,set)
  var=s$var; meth=s$method; loc=s$location; set=s$set
  # Generating metadata from calendar window
  out=c()
  sy=s$scan[[paste0("Y",y)]]
  if (reverse.prob) {sy$p1=1-sy$p1; sy$p0=1-sy$p0}
  win=s$windows[iw]
  day1=id-trunc((win-1)/2); day2=day1+win-1
  date1=mdy2yyyymmdd(time2mdy(day1,paste0("days since ",y-1,"-12-31"),"noleap"))
  date2=mdy2yyyymmdd(time2mdy(day2,paste0("days since ",y-1,"-12-31"),"noleap"))
  ymd=time2mdy(id,paste0("days since ",y-1,"-12-31"),"noleap")
  out=rbind(out,c(sy$p1[id,iw],sy$p0[id,iw],day1,day2,date1,date2,win,ymd$y,ymd$m,ymd$d,id,iw))
  out=as.data.frame(out,row.names="")
  names(out)=c("p1","p0","day1","day2","date1","date2","win","y","m","d","id","iw")
  # Returns a list with metadata
  output=list()
  for (n in c("set","var","location","years","period","scan.years","yref","windows","method","params"))
    output[[paste(n)]]=s[[paste(n)]]
  output$selparams=list(years=y,id=id,iw=iw,reverse.prob=reverse.prob)
  output$events=out
  output
}

#---------------------------------------------------------------------------------
# Identifies all calendar events for a given date over a set of wraf regions
#---------------------------------------------------------------------------------
get.all.cal.events=function(var="TM",set="era5",wset=2,y=2003,id=220,iw=9,reverse.prob=F){
  # Get all scans and max events
  out=c()
  for (loc in REGIONS[[paste0("WRAF",wset)]][1,]){
    s=get.scan.1d(var,"calend",loc,set)
    dum=get.cal.event(s,y=y,id=id,iw=iw,reverse.prob=reverse.prob)$events
    dum$loc=loc
    out=rbind(out,dum)
  }
  # Returns a list with metadata
  output=list()
  for (n in c("set","var","years","period","scan.years","yref","windows","method","params"))
    output[[paste(n)]]=s[[paste(n)]]
  output$wraf.set=wset
  output$selparams=list(years=y,id=id,iw=iw,reverse.prob=reverse.prob)
  output$events=out
  output
}

#---------------------------------------------------------------------------------
# Gets calendar anomaly for a given event
# Warning: s must be a calendar scan!
#---------------------------------------------------------------------------------
get.anom.eve=function(s=NULL,var='TM',loc='France',set='mf',
                      y=2003,id=220,iw=9){
  # Data
  if (is.null(s)) s=get.scan.1d(var,"calend",loc,set)
  var=s$var; meth=s$method; loc=s$location; set=s$set
  # Values
  iy=which(s$period==y)
  value=s$scan[[paste0("Y",y)]]$value[id,iw]
  trend=s$trend[id,iw,]
  iref=which(s$period==s$yref)
  delta=trend[iy]-trend[iref]
  # Anomalies
  fit1=s$scan[[paste0("Y",y)]]$fit1[id,iw,]
  mu=fit1[1]; sd=fit1[2]
  ano=value-mu
  anonorm=ano/sd
  # Output
  data.frame(value=value,delta=delta,mu=mu,sd=sd,ano=ano,anonorm=anonorm)
}

#---------------------------------------------------------------------------------
# Gets fit and trend for a given event and a given scan
#---------------------------------------------------------------------------------
get.fit.eve=function(s=NULL,var="TM",meth="calend",loc="France",set="mf",
                      y=2003,id=220,iw=9){
  # Data
  if (is.null(s)) s=get.scan.1d(var,meth,loc,set)
  var=s$var; meth=s$method; loc=s$location; set=s$set
  # Values
  iy=which(s$period==y)
  value=s$scan[[paste0("Y",y)]]$value[id,iw]
  trend=s$trend[id,iw,]
  iref=which(s$period==s$yref)
  delta=trend[iy]-trend[iref]
  # Anomalies
  fit1=s$scan[[paste0("Y",y)]]$fit1[id,iw,]
  mu=fit1[1]; sigma=fit1[2]; ksi=fit1[3]
  # Output
  data.frame(value=value,delta=delta,mu=mu,sigma=sigma,ksi=ksi)
}

#---------------------------------------------------------------------------------
# Recomputes stats with confidence interval for a given event
#---------------------------------------------------------------------------------
get.stat.eve=function(s=NULL,var='TM',meth='annmax',loc='France',set='mf',
                      y=2003,id=220,iw=9,nboot=1000,level=0.9){
  # Data
  if (is.null(s)) s=get.scan.1d(var,meth,loc,set)
  var=s$var; meth=s$method; loc=s$location; set=s$set
  # Sample + value
  iy=which(s$period==y)
  value=s$scan[[paste0("Y",y)]]$value[id,iw]
  if (s$params$annual.maxima) data=s$data[1,iw,] else data=s$data[id,iw,]
  trend=s$trend[id,iw,]
  sample=na.omit(data-trend+trend[iy])
  iref=which(s$period==s$yref)
  delta=trend[iy]-trend[iref]
  # Arguments for fits / percentile levels
  dis=s$params$distrib
  ksi=s$params$fixed.ksi
  exc=s$params$exceedance
  # Best estimate
  be.fit1=myfitparams(sample,dis,ksi)
  be.fit0=myfitparams(sample-delta,dis,ksi)
  be.p1=myperclevel(value,dis,be.fit1,exc)
  be.p0=myperclevel(value,dis,be.fit0,exc)
  be.pr=(1-be.p1)/(1-be.p0)
  # Bootstrap
  bs.p1=c(); bs.p0=c(); bs.pr=c()
  for (i in 1:nboot){
    dum=sample(sample,length(sample),replace=T)
    f1=myfitparams(dum,dis,ksi)
    f0=myfitparams(dum-delta,dis,ksi)
    p1=myperclevel(value,dis,f1,exc)
    p0=myperclevel(value,dis,f0,exc)
    pr=(1-p1)/(1-p0)
    bs.p1=c(bs.p1,p1); bs.p0=c(bs.p0,p0); bs.pr=c(bs.pr,pr)
  }
  # Output
  out=apply(cbind(bs.p1,bs.p0,bs.pr),2,quantile,c(0.5,(1-level)/2,1-(1-level)/2),na.rm=T)
  out=rbind(c(be.p1,be.p0,be.pr),out)
  out=as.data.frame(out)
  names(out)=c("p1","p0","pr")
  out
}

#---------------------------------------------------------------------------------
# Converts list of events into a latex table
#---------------------------------------------------------------------------------
events2tab.anom=function(e,n=NULL,stat.max="p1",loc=F,locname=F){
  # Selection of n max events
  if (is.null(n)) sel=e$events
  if (! is.null(n)){
    n=min(n,nrow(e$events))
    sel=e$events[order(e$events[,paste(stat.max)],decreasing=T),][1:n,]
  }
  # Add info for table
  if (! "loc" %in% names(sel)) sel$loc=rep(e$location,nrow(sel))
  sel$n=paste0(1:nrow(sel),".")
  sel$dleg=sel$val=sel$ano=sel$n
  for (i in 1:nrow(sel)){ si=sel[i,]
    cday=mdy2yyyymmdd(si[c("m","d","y")])
    sel$dleg[i]=dates2leg(cday,cday,short=T)
    ano=get.anom.eve(var=e$var,loc=si$loc,set=e$set,y=si$y,id=mymod(si$id,365),iw=si$iw)
    sel$val[i]=round(ano$value,1)
    sel$ano[i]=signif(ano$anonorm,2)
  }
  # Generate table
  if (!loc) tab=matrix(paste(sel$n,"&",sel$y,"&",sel$dleg,"&",sel$win,"&",sel$val,"&",sel$ano,"\\\\"),nrow(sel),1)
  if (loc & !locname) tab=matrix(paste(sel$n,"&",sel$loc,"&",sel$y,"&",sel$dleg,"&",sel$win,"&",sel$val,"&",sel$ano,"\\\\"),nrow(sel),1)
  if (!loc & locname) tab=matrix(paste(sel$n,"&",wraf.name(sel$loc),"&",sel$y,"&",sel$dleg,"&",sel$win,"&",sel$val,"&",sel$ano,"\\\\"),nrow(sel),1)
  if (loc & locname)tab=matrix(paste(sel$n,"&",sel$loc,"&",wraf.name(sel$loc),"&",sel$y,"&",sel$dleg,"&",sel$win,"&",sel$val,"&",sel$ano,"\\\\"),nrow(sel),1)
  tab
}

#---------------------------------------------------------------------------------
# Converts list of events into a latex table
#---------------------------------------------------------------------------------
events2tab.stat=function(e,n=NULL,stat.max="p1",loc=F,locname=F,confint=F,nboot=1000,level=0.9){
  # Selection of n max events
  if (is.null(n)) sel=e$events
  if (! is.null(n)){
    n=min(n,nrow(e$events))
    sel=e$events[order(e$events[,paste(stat.max)],decreasing=T),][1:n,]
  }
  # Add info for table
  if (! "loc" %in% names(sel)) sel$loc=rep(e$location,nrow(sel))
  sel$n=paste0(1:nrow(sel),".")
  sel$dleg=sel$n
  for (i in 1:nrow(sel)) sel$dleg[i]=dates2leg(sel$date1[i],sel$date2[i],short=T)
  # Stats without confidence intervals
  sel$pp1=signif(1/(1-sel$p1),2)
  sel$pp0=signif(1/(1-sel$p0),2)
  sel$ppr=signif((1-sel$p1)/(1-sel$p0),2)
  # Confidence intervals if requested
  if (confint){
    for (i in 1:nrow(sel)){ si=sel[i,]
    ci=get.stat.eve(var=e$var,meth=e$method,loc=si$loc,set=e$set,
                    y=si$y,id=mymod(si$id,365),iw=si$iw,nboot=nboot,level=level)
    ci[,1:2]=signif(1/(1-ci[,1:2]),2); ci$pr=signif(ci$pr,2)
    sel$pp1[i]=paste0(sel$pp1[i]," [",ci$p1[3]," to ",ci$p1[4],"]")
    sel$pp0[i]=paste0(sel$pp0[i]," [",ci$p0[3]," to ",ci$p0[4],"]")
    sel$ppr[i]=paste0(sel$ppr[i]," [",ci$pr[3]," to ",ci$pr[4],"]")
    }
  }
  # Generate table
  if (!loc) tab=matrix(paste(sel$n,"&",sel$y,"&",sel$dleg,"&",sel$pp1,"&",sel$pp0,"&",sel$ppr,"\\\\"),nrow(sel),1)
  if (loc & !locname) tab=matrix(paste(sel$n,"&",sel$loc,"&",sel$y,"&",sel$dleg,"&",sel$pp1,"&",sel$pp0,"&",sel$ppr,"\\\\"),nrow(sel),1)
  if (!loc & locname) tab=matrix(paste(sel$n,"&",wraf.name(sel$loc),"&",sel$y,"&",sel$dleg,"&",sel$pp1,"&",sel$pp0,"&",sel$ppr,"\\\\"),nrow(sel),1)
  if (loc & locname) tab=matrix(paste(sel$n,"&",sel$loc,"&",wraf.name(sel$loc),"&",sel$y,"&",sel$dleg,"&",sel$pp1,"&",sel$pp0,"&",sel$ppr,"\\\\"),nrow(sel),1)
  tab
}