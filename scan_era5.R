#---------------------------------------------------------------------------------
# Source libraries / functions
#---------------------------------------------------------------------------------
source("scan_source.R")

#---------------------------------------------------------------------------------
# Namelist of the event
#---------------------------------------------------------------------------------
# Variable of interest
VAR="TM"
#VAR="TX"
#VAR="TN"

##################################################################################
for (VAR in c("TM","TX","TN")){ # LOOP ON VARIABLE
##################################################################################

# Years used for data samples
YEARS=1959:2022

# Years to scan (possibly a subset of YEARS)
SCAN.YRS=YEARS

# First/last years for annual-max/min analysis (are Txnday or Tnnday representative?)
LAST.YEAR.MAX=(max(YEARS)==2022)
LAST.YEAR.MIN=(max(YEARS)==2022)

# Time windows for searching the event (durations in days)
WINDOWS=c(1:8,10,12,15,20,25,30,40,50,60,75,90,120,180,270,365)
nWIN=length(WINDOWS)

# Location (example for initialization, see loop below)
LOC="X.2.2"

# Year of reference for counter-factual world
YREF=1959

# Output directory
ODIR=paste0(SCANDIR,"scans_1d/era5/")
system(paste("mkdir -p",ODIR))

#---------------------------------------------------------------------------------
# Get data: XLOC data.frame $date $var = daily values at LOC
#---------------------------------------------------------------------------------
# Local data for 1D-analysis (ERA5)
if (VAR=="TM") obsfile="../data/era5/tas_era5_day_raw_1959-2022.nc"
if (VAR=="TX") obsfile="../data/era5/tasmax_era5_day_raw_1959-2022.nc"
if (VAR=="TN") obsfile="../data/era5/tasmin_era5_day_raw_1959-2022.nc"
XLOC=get.wraf.data(obsfile,LOC)
names(XLOC)=c("date","var")
XLOC$var=XLOC$var-273.15

# Restrict to YEARS
XLOC=XLOC[which(yyyymmdd2mdy(XLOC$date)$y %in% YEARS),]

# Add a year of NAs before/after (edge effects) + delete feb29s
PERIOD=(min(YEARS)-1):(max(YEARS)+1)
dum=data.frame(date=mdy2yyyymmdd(makedates(y=PERIOD,avg="day",caltype="noleap")),var=NA)
dum$var[which(dum$date %in% XLOC$date)]=XLOC$var[which(XLOC$date %in% dum$date)]
XLOC=dum

#---------------------------------------------------------------------------------
# Get forced response: TREND vector of yearly values over PERIOD 
#---------------------------------------------------------------------------------
# Estimate forced response at LOC (from XMMM = multi-model mean)
XMMM=get.wraf.data(cmipfile,LOC)[,2]-273.15
TREND=forced_response_by_gam_fit(XMMM,Enat,10)$all[which(CMIP.PERIOD %in% PERIOD)]

#---------------------------------------------------------------------------------
# A few helpful assigns used in other R subroutines
#---------------------------------------------------------------------------------
nWIN=length(WINDOWS)                        # number of time windows
nYRS=length(PERIOD)                         # number of years
DATES=yyyymmdd2mdy(XLOC$date)               # dates as data.frame $m $d $y
nD0=365                                     # number of days of the event period (here, whole year)
iYREF=which(PERIOD==YREF)                   # index of the year of reference

#---------------------------------------------------------------------------------
# Scan 1d of a region
#---------------------------------------------------------------------------------
if (2==1){ ###

# Load functions
source("scan_1d_sub.R")

# Method to use
METHOD="calend"

# Source scan subroutine
source("scan_1d_exe.R")

# Save result
save(list="SCAN",file=paste0(ODIR,"scan_",LOC,"_",VAR,"_",METHOD,".Rdata"))

} ###

#---------------------------------------------------------------------------------
# All scans 1d (loop on sets/regions)
#---------------------------------------------------------------------------------
if (2==1){ ###

# Load functions
source("scan_1d_sub.R")

# Loop on regions
for (w in WRAF.SETS){
for (loc in REGIONS[[paste0("WRAF",w)]][1,]){

  # Location
  LOC=loc

  # Get data at LOC, same selection of dates as above
  XLOC=get.wraf.data(obsfile,LOC)
  names(XLOC)=c("date","var")
  XLOC$var=XLOC$var-273.15
  XLOC=XLOC[which(yyyymmdd2mdy(XLOC$date)$y %in% YEARS),]
  dum=data.frame(date=mdy2yyyymmdd(DATES),var=NA)
  dum$var[which(dum$date %in% XLOC$date)]=XLOC$var[which(XLOC$date %in% dum$date)]
  XLOC=dum

  # Get trend at LOC
  XMMM=get.wraf.data(cmipfile,LOC)[,2]-273.15
  TREND=forced_response_by_gam_fit(XMMM,Enat,10)$all[which(CMIP.PERIOD %in% PERIOD)]

  # Scan with loop on methods
  for (m in c("calend",#"locmax","locmin",
              "annmax",#"annmax1","annmax2",
              "annmin")){#,"annmin1")){,"annmin2")){
    # Method to use
    METHOD=m
    # Source scan subroutine
    source("scan_1d_exe.R")
    # Save result
    save(list="SCAN",file=paste0(ODIR,"scan_",LOC,"_",VAR,"_",METHOD,".Rdata"))
  }

}}

} ###

#---------------------------------------------------------------------------------
# Scans 2d (loop on sets)
#---------------------------------------------------------------------------------
#if (2==1){ ###

# Output directory
ODIR2=paste0(SCANDIR,"scans_2d/era5/")
system(paste("mkdir -p",ODIR2))  

# Load functions
source("scan_functions.R")
  
# Loop on regions
for (w in WRAF.SETS){
    
# Loop on methods
for (m in c("annmax","annmin")){

  # Select the max event per year
  for (y in SCAN.YRS){
    outy=get.all.max.events(var=VAR,meth=m,wset=w,years=y,only.max=T,max.win=30)
    if (y==YEARS[1]) {out=outy; out$selparams$years=NULL}
    else out$events=abind(out$events,outy$events,along=3)
  }
  YEARLY.EVENTS=out

  # Select all events with high percentile values
  out=get.all.max.events(var=VAR,meth=m,wset=w,only.max=F,max.win=30,min.value=0.95)
  MAX.EVENTS=out
      
  # Save result
  save(list=c("YEARLY.EVENTS","MAX.EVENTS"),file=paste0(ODIR2,"events_WRAF",w,"_",VAR,"_",m,".Rdata"))
    
}}
  
#} ###

##################################################################################
} ### END LOOP ON VARIABLE
##################################################################################

