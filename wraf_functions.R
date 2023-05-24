
# WRAF sets (sizes) used
WRAF.SETS=c("10","5","2","05")

# WRAF regions information
REGIONS=list(); for (wset in WRAF.SETS){
  ifile=paste0(WRAFDIR,"region_fx-WRAF",wset,"-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc")
  dum=rbind(myno(ifile,"regionid")$var,myno(ifile,"nameshort")$var)
  REGIONS[[paste0("WRAF",wset)]]=dum
}

# Wraf regions identifiers
REGIDS=c(); for (wset in WRAF.SETS) REGIDS=c(REGIDS,REGIONS[[paste0("WRAF",wset)]][1,])

# COUNTRIES in WRAF format
COUNTRIES=myno(paste0(WRAFDIR,"countries_in_WRAF-v4-1_format.nc"),"nameshort")$var
COUNTRIES=gsub(" ",".",COUNTRIES)

##################################################

# Function that gives the set of a region
wraf.set=function(regids){
  out=c(); for (regid in regids){
    if (regid %in% REGIDS)    out=c(out,WRAF.SETS[length(strsplit(regid,".",fixed=T)[[1]])])
    if (regid %in% COUNTRIES) out=c(out,"cn")
  }
  out
}

# Functions that gives the name of a region id
wraf.name=function(regids){
  out=regids
  wset=wraf.set(regids)
  for (i in 1:length(out)){
    if (wset[i] %in% WRAF.SETS){
      dum=REGIONS[[paste0("WRAF",wset[i])]]
      out[i]=dum[2,which(dum[1,]==regids[i])]
    }
  }
  out
}

# Function that gives the family of a region
wraf.family=function(regid){
  out=regid
  if (regid %in% REGIDS){
    dum=strsplit(regid,".",fixed=T)[[1]]
    parents=c(); if (length(dum)>1){ for (i in 1:(length(dum)-1)) parents=c(parents,paste(dum[1:i],collapse="."))}
    children=REGIDS[startsWith(REGIDS,paste0(regid,"."))]
    out=c(parents,regid,children)
    out=out[which(out%in%REGIDS)]
  }
  out  
}

# Function that gives the extended family of a region
wraf.extended.family=function(regid,root.wset="5"){
  out=regid
  if (regid %in% REGIDS){
    dum=strsplit(regid,".",fixed=T)[[1]]
    root=min(length(dum),which(WRAF.SETS==root.wset))
    root.reg="?"; while(! root.reg %in% REGIDS){
      root.reg=paste(dum[1:root],collapse=".")
      root=root+1
    }
    out=wraf.family(root.reg)
  }
  out
}

# Function that gets data over a given region from pre-computed WRAF files associated with ifile
# (no filter on NAs, typically for ERA5 or CMIP data)
get.wraf.data=function(ifile,regid){
  if (wraf.set(regid)!="cn") dat.file=paste0(ifile,".WRAF",wraf.set(regid),".txt")
  else dat.file=paste0(ifile,".COUNTRIES.txt")
  as.data.frame(fread(dat.file,select=c("date",regid)))
}

# Same with a filter on the fraction of data vs. NAs (in %) + nb of data vs. NAs per year (in days)
# (typically for BEST data)
get.wraf.data.with.filters=function(ifile,regid,not.nas.frac=NULL,not.nas.days=NULL){
  nas.file=paste0(ifile,".NAs.WRAF",wraf.set(regid),".txt")
  if (!is.null(not.nas.frac) & !file.exists(nas.file)) {
    print("Error: NAs file is not present.")
    return(NULL)
  }
  if (wraf.set(regid)!="cn") dat.file=paste0(ifile,".WRAF",wraf.set(regid),".txt")
  else dat.file=paste0(ifile,".COUNTRIES.txt")
  out=as.data.frame(fread(dat.file,select=c("date",regid)))
  if (!is.null(not.nas.frac)){
    nas=as.data.frame(fread(nas.file,select=c("date",regid)))
    out[which(is.na(nas[,2]) | nas[,2]<not.nas.frac),2]=NA
  }
  if (!is.null(not.nas.days)){
    dts=yyyymmdd2mdy(out[,1])
    yrs=unique(dts)$y
    for (y in yrs){
      if (sum(!is.na(out[which(dts$y==y),2]))<not.nas.days) out[which(dts$y==y),2]=NA
    }
  }
  out
}

# Function that pre-processes NA filters on WRAF files
filter.wraf.data=function(ifile,wset,not.nas.frac=90,not.nas.days=350){
  dat.file=paste0(ifile,".WRAF",wset,".txt")
  nas.file=paste0(ifile,".NAs.WRAF",wset,".txt")
  out=as.data.frame(fread(dat.file))
  if (!is.null(not.nas.frac) & file.exists(nas.file)){
    nas=as.data.frame(fread(nas.file))
    out[is.na(nas) | nas<not.nas.frac]=NA
  }
  if (!is.null(not.nas.days)){
    dts=yyyymmdd2mdy(out[,1])
    yrs=unique(dts)$y
    for (i in 2:ncol(out)){ for (y in yrs){
      if (sum(!is.na(out[which(dts$y==y),i]))<not.nas.days) out[which(dts$y==y),i]=NA
    }}
  }
  outfile=paste0(ifile,".p",not.nas.frac,"d",not.nas.days,".WRAF",wset,".txt")
  write.table(out,outfile,quote=F,row.names=F,col.names=T)
}




