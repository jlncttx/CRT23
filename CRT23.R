setwd("/home/cattiaux/Articles/SCAN2021/R")

source("scan_source.R")
source("scan_figures.R")
source("scan_functions.R")

figdir="../article/figures/"
tabdir="../article/tables/"

#-------------------------------------------------------------------------------------------
# Local analysis (France)
#-------------------------------------------------------------------------------------------

# Figure A1: detrending procedure
pdf(paste0(figdir,"fig_detrend.pdf"),width=8.5,height=6)
myfig.detrend(var='TM',loc='France',set='mf',y=2022,f.df=12,h.df=6,title=c("a)","b)","c)","d)"))
dev.off()

# Figure C1: example of an event
# France
pdf(paste0(figdir,"fig_event_France.pdf"),width=9,height=7.5)
layout(cbind(c(1,2),c(1,3),c(5,4)),widths=c(1,1,0.2))
myfig.eve(var='TM',loc='France',set='mf',y=2022,id=168,iw=5,title="a)")
myfig.1d(var='TM',loc='France',set='mf',meth = "calend",y=2022,min.win=3,title="b)",colorbar=F,show.events=5,min.value=0.99)
myfig.1d(var='TM',loc='France',set='mf',meth = "annmax",y=2022,max.win=30,title="c)",colorbar=F,show.events=5,min.value=0.75)
myfig.1d(image=F,colorbar=T)
dev.off()
# HW 1.2.1
pdf(paste0(figdir,"fig_event_hw.pdf"),width=9,height=7.5)
layout(cbind(c(1,2),c(1,3),c(5,4)),widths=c(1,1,0.2))
myfig.eve(var='TM',loc='1.2.1',set='era5',y=2021,id=181,iw=3,title="a)")
myfig.1d(var='TM',loc='1.2.1',set='era5',meth = "calend",y=2021,min.win=3,title="b)",colorbar=F,show.events=5,min.value=0.99)
myfig.1d(var='TM',loc='1.2.1',set='era5',meth = "annmax",y=2021,max.win=30,title="c)",colorbar=F,show.events=5,min.value=0.75)
myfig.1d(image=F,colorbar=T)
dev.off()
# CS 1.2.1
pdf(paste0(figdir,"fig_event_cs.pdf"),width=9,height=7.5)
layout(cbind(c(1,2),c(1,3),c(5,4)),widths=c(1,1,0.2))
myfig.eve(var='TM',loc='1.2.1',set='era5',y=2022,id=355,iw=5,title="a)")
myfig.1d(var='TM',loc='1.2.1',set='era5',meth = "calend",y=2022,min.win=3,title="b)",reverse.prob=T,colorbar=F,show.events=5,min.value=0.99)
myfig.1d(var='TM',loc='1.2.1',set='era5',meth = "annmin",y=2022,max.win=30,title="c)",colorbar=F,show.events=5,min.value=0.75)
myfig.1d(image=F,colorbar=T)
dev.off()

# Figure 1: methodology
# Calend
pdf(paste0(figdir,"fig_meth_calend.pdf"),width=9,height=3)
myfig.meth(var='TM',meth='calend',loc='France',set='mf',y=2022,id=168,iw=5,title=c("a)","b)"),reverse.prob=F)
dev.off()
# Annmax
pdf(paste0(figdir,"fig_meth_annmax.pdf"),width=9,height=3)
myfig.meth(var='TM',meth='annmax',loc='France',set='mf',y=2022,id=168,iw=5,title=c("c)","d)"),reverse.prob=F)
dev.off()

# Figure 2: local chronology
hw1=get.max.events(var='TM',meth='calend',loc='France',set='mf',min.value=0.99,min.win=3)
hw2=get.max.events(var='TM',meth='annmax',loc='France',set='mf',min.value=0.75,max.win=30)
cs1=get.max.events(var='TM',meth='calend',loc='France',set='mf',min.value=0.99,min.win=3,reverse.prob=T)
cs2=get.max.events(var='TM',meth='annmin',loc='France',set='mf',min.value=0.75,max.win=30)

pdf(paste0(figdir,"fig_1d.pdf"),width=6,height=7.5)
layout(cbind(c(1,3),c(2,4),c(6,5)),widths=c(1,0.5,0.2))
myfig.chronology(hw1,title="a) Heat waves |",colorbar=F,show=5)
myfig.chronology(hw2,title=NULL,xlim=c(160,250),colorbar=F,yaxis=F,show=5)
myfig.chronology(cs1,title="b) Cold spells |",xlim=c(182,181),colorbar=F,show=5)
myfig.chronology(cs2,title=NULL,xlim=c(325,80),colorbar=F,yaxis=F,show=5)
myfig.chronology(cs2,image=F)
dev.off()

# Table 1: calendar events
a=get.scan.1d(var='TM',meth='calend',loc='France',set='mf')
hot=get.max.events(a,only.months=1,only.max=T,min.win=3,max.win=30)
for (m in 2:12) hot$events=rbind(hot$events,get.max.events(a,only.months=m,only.max=T,min.win=3,max.win=30)$events)
cld=get.max.events(a,only.months=1,only.max=T,min.win=3,max.win=30,reverse.prob=T)
for (m in 2:12) cld$events=rbind(cld$events,get.max.events(a,only.months=m,only.max=T,min.win=3,max.win=30,reverse.prob=T)$events)
#
write.table(events2tab.anom(hot),paste0(tabdir,"tab_1d_hot.tex"),quote=F,row.names=F,col.names=F)
write.table(events2tab.anom(cld),paste0(tabdir,"tab_1d_cld.tex"),quote=F,row.names=F,col.names=F)

# Table 2: all-time events
write.table(events2tab.stat(hw2,5,confint=T,nboot=10000),paste0(tabdir,"tab_1d_hw.tex"),quote=F,row.names=F,col.names=F)
write.table(events2tab.stat(cs2,5,confint=T,nboot=10000),paste0(tabdir,"tab_1d_cs.tex"),quote=F,row.names=F,col.names=F)

# For the text: max events per months/windows
a=get.scan.1d(meth='calend')
all=c(); for (m in 1:12) all=rbind(all,get.max.events(a,only.months=m,only.max=T,min.win=3,max.win=30)$events)
a=get.scan.1d(meth='annmax')
all=c(); for (w in a$windows) all=rbind(all,get.max.events(a,only.max=T,max.win=w,min.win=w)$events)
a=get.scan.1d(meth='annmin')
all=c(); for (w in a$windows) all=rbind(all,get.max.events(a,only.max=T,max.win=w,min.win=w)$events)

#-------------------------------------------------------------------------------------------
# Global analysis
#-------------------------------------------------------------------------------------------

# Figure 3: HW at 2 Mkm2
ex1=get.all.cal.events(y=2021,id=181,iw=3)
yy=list(2021,2003:2022,1959:2022)
tt=c("a)","c)","d)")
pdf(paste0(figdir,"fig_2d_hw.pdf"),width=9,height=4.5)
layout(cbind(c(2,3),c(1,4),c(6,5)),widths=c(1,1,0.15))
myfig.2d(ex1,title="b)",colorbar=F)
for (i in 1:length(yy)){
  ee=get.scan.2d("TM","annmax","era5","2",yy[[i]])
  myfig.2d(ee,title=tt[i],colorbar=F)
}
myfig.2d(ee,image=F)
dev.off()

# Figure 4: CS at 2 Mkm2
ex2=get.all.cal.events(y=2022,id=355,iw=5,reverse.prob = T)
yy=list(2022,2003:2022,1959:2022)
tt=c("a)","c)","d)")
pdf(paste0(figdir,"fig_2d_cs.pdf"),width=9,height=4.5)
layout(cbind(c(2,3),c(1,4),c(6,5)),widths=c(1,1,0.15))
myfig.2d(ex2,title="b)",colorbar=F)
for (i in 1:length(yy)){
  ee=get.scan.2d("TM","annmin","era5","2",yy[[i]])
  myfig.2d(ee,title=tt[i],colorbar=F)
}
myfig.2d(ee,image=F)
dev.off()

# Table 3: events at 2 Mkm2
hw3=get.scan.2d("TM","annmax","era5","2",2003:2022)
cs3=get.scan.2d("TM","annmin","era5","2",2003:2022)
write.table(events2tab.stat(hw3,10,loc=T,locname=T,confint=T,nboot=1000),paste0(tabdir,"tab_2d_hw_wraf2.tex"),
            quote=F,row.names=F,col.names=F)
write.table(events2tab.stat(cs3,10,loc=T,locname=T,confint=T,nboot=1000),paste0(tabdir,"tab_2d_cs_wraf2.tex"),
            quote=F,row.names=F,col.names=F)

# Figure 5: zoom for 3 events
pdf(paste0(figdir,"fig_zoom_ggplot.pdf"),width=12,height=8)
#layout(cbind(c(1,3,5,7),c(2,4,6,8),c(10,10,10,9)),widths=c(1,1,0.15))
crs.Can=crs.lambert(-105,55,3000,3000)
crs.Eur=crs.lambert(10,53,2000,2000)
crs.Rus=crs.lambert(45,60,2500,2500)
pan=list(); i=1
for (w in WRAF.SETS){
  tits=c("","","")
  if (w==WRAF.SETS[1]) tits=c("a)","b) ","c)")
  if (w==WRAF.SETS[2]) tits=c("",wraf.name("7.1"),"")
  if (w==WRAF.SETS[3]) tits=c("","",wraf.name("X.2.2"))
  if (w==WRAF.SETS[4]) tits=c(wraf.name("1.2.1.3"),"","")
  ltit="" ; if (w==WRAF.SETS[1]) ltit="default"
  rtit=w
  pan[[i]] = myfig.2d.ggplot(get.scan.2d("TM","annmax","era5",w,2021),crs=crs.Can,colorbar=F,
                             title=tits[1],ltitle=ltit,rtitle=rtit)
  pan[[i+1]]=myfig.2d.ggplot(get.scan.2d("TM","annmax","era5",w,2010),crs=crs.Rus,colorbar=F,
                             title=tits[2],ltitle=ltit,rtitle=rtit)
  pan[[i+2]]=myfig.2d.ggplot(get.scan.2d("TM","annmax","era5",w,2003),crs=crs.Eur,colorbar=F,
                             title=tits[3],ltitle=ltit,rtitle=rtit)
#  pan[[i+1]]=myfig.2d.ggplot(get.scan.2d("TM","annmin","era5",w,2022),crs=crs.Can,colorbar=F)
  i=i+3
}
pan[[i]]=NULL
pan[[i+1]]=NULL
pan[[i+2]]=function(){
  par(mar=c(0,0.5,0,0.5))
  myfig.2d.ggplot(get.scan.2d("TM","annmin","era5",w,2022),crs=crs,image=F)
}
plot_grid(plotlist=pan,nrow=3,byrow=F,rel_widths=c(1,1,1,1,0.5),axis="t")
dev.off()

# Table 4: all events with filter on wraf families
family.filter=function(e,root.wset="5"){
  out=c(); dum=e; while(nrow(dum)>0){
    eve=dum[1,]
    fam=which(dum$loc %in% wraf.extended.family(eve$loc))
    fam=dum[fam,]; del=c(); for (i in 1:nrow(fam)){
      if (any(eve$date1:eve$date2 %in% fam$date1[i]:fam$date2[i])) del=c(del,which(dum$loc==fam$loc[i]))
    }
    out=rbind(out,eve)
    dum=dum[-del,]
  }
  out
}
# HW
per=2003:2022
hw4=get.scan.2d("TM","annmax","era5","10",per)
sel=hw4$events
for (wset in WRAF.SETS[-1]) sel=rbind(sel,get.scan.2d("TM","annmax","era5",wset,per)$events)
sel=sel[order(sel$p1,decreasing=T),]
hw4$events=family.filter(sel)
write.table(events2tab.stat(hw4,10,loc=T,locname=T,confint=T,nboot=1000),paste0(tabdir,"tab_2d_hw_all.tex"),
            quote=F,row.names=F,col.names=F)
# CS
per=2003:2022
cs4=get.scan.2d("TM","annmin","era5","10",per)
sel=cs4$events
for (wset in WRAF.SETS[-1]) sel=rbind(sel,get.scan.2d("TM","annmin","era5",wset,per)$events)
sel=sel[order(sel$p1,decreasing=T),]
cs4$events=family.filter(sel)
write.table(events2tab.stat(cs4,10,loc=T,locname=T,confint=T,nboot=1000),paste0(tabdir,"tab_2d_cs_all.tex"),
            quote=F,row.names=F,col.names=F)

# Figure B1: histograms of p1 over all events
pdf(paste0(figdir,"hist_p1.pdf"),width=9,height=6)
par(mar=c(2,2.5,1.5,0.5),mgp=c(1.5,0.5,0),tcl=0.2,las=0)
layout(matrix(1:9,3,3,byrow=T))
s=get.scan.1d(loc="1",set="era5")
for (m in c("calend","annmax","annmin")){
  if (m=="calend") {iw=c(3,11,19); tit=c("a)","b)","c)")}
  if (m=="annmax") {iw=c(1,5,14) ; tit=c("d)","e)","f)")}
  if (m=="annmin") {iw=c(1,5,14) ; tit=c("g)","h)","i)")}
  win=s$windows[iw]
  for (i in 1:3){
    vals=1-hist.all.stat.2d(var="TM",meth=m,set="era5",wset=2,stat="p1",iw=iw[i]) # percentile levels vs. prob
    hist(vals,breaks=seq(0,1,0.05),main="",xlab="",col=method.color2(m))
    title(paste0(tit[i]," p1 ",win[i],"-day"),adj=0)
    if (i==1) title("2 Mkm2 regions",adj=1)
    if (i==3) title(method.name(m),col.main=method.color(m),adj=1)
  }
}
dev.off()

# Figure B2a: histogram of selected years
pdf(paste0(figdir,"hist_recs.pdf"),width=6,height=6)
par(mar=c(2.2,2.2,1.5,0.5),mgp=c(1.2,0.2,0),tcl=0.5,las=0)
layout(matrix(1:2,2,1))
for (m in c("annmax","annmin")){
  h1=hist.all.year.2d(var="TM",meth=m,set="era5",wset="2")
  h2=hist.all.year.2d(var="TM",meth=m,set="era5",wset="05")
  hh1=hist(h1,breaks=1958:2022+0.5,plot=F)
  hh2=hist(h2,breaks=1958:2022+0.5,plot=F)
  ylim=c(0,1.1*max(c(hh1$density,hh2$density)))
  hist(h1,breaks=1958:2022+0.5,prob=T,col=method.color2(m),main="",xlab="Year",ylim=ylim)
  lines(hh2$mids,hh2$density,type="h",lwd=2)
  axis(1,at=seq(1960,2020,10))
  axis(1,at=seq(1958,2024,2),tcl=0.25,label=NA)
  axis(2)
  if (m=="annmax") title("a) Year of record heat waves",adj=0)
  if (m=="annmin") title("c) Year of record cold spells",adj=0)
  title(method.name(m),col.main=method.color(m),adj=1)
  legend("topright",bty="n",lty=c(NA,1),lwd=c(0,2),pch=c(22,NA),pt.cex=2,
         pt.bg=method.color2(m),leg=c("2 Mkm2","0.5 Mkm2"))
}
dev.off()

# Figure B2b: time series of p1
pdf(paste0(figdir,"ts_p1.pdf"),width=6,height=6)
par(mar=c(2.2,2.2,1.5,0.5),mgp=c(1.2,0.2,0),tcl=0.5,las=0)
layout(matrix(1:2,2,1))
for (m in c("annmax","annmin")){
  vals=1-get.yearly.stat.2d(var="TM",meth=m,set="era5",wset="2",stat="p1") # percentile levels vs. p1
  vals=t(apply(vals,2,quantile,c(0.5,0.05)))
  matplot(1959:2022,vals,type="l",lty=c(1,2),lwd=c(2,1),col=method.color(m),
          main="",xlab="Year",ylab="p1",ylim=c(0,1.1*max(vals)))
  if (m=="annmax") title("b) Yearly min p1 | 2 Mkm2",adj=0)
  if (m=="annmin") title("d) Yearly min p1 | 2 Mkm2",adj=0)
  title(method.name(m),col.main=method.color(m),adj=1)
  legend("topleft",bty="n",leg="Median",lty=1,lwd=2,col=method.color(m))
  legend("topright",bty="n",leg="Quantile 5%",lty=2,lwd=1,col=method.color(m))
}
dev.off()

# Figure B3: examples of data / fits issues
# Annmax
pdf(paste0(figdir,"bad_data_annmax.pdf"),width=9,height=3)
myfig.meth(var='TM',meth='annmax',loc='3.2',set='era5',y=1963,id=294,iw=12,title=c("a)","b)"),reverse.prob=F)
dev.off()
# Annmin
pdf(paste0(figdir,"bad_data_annmin.pdf"),width=9,height=3)
myfig.meth(var='TM',meth='annmin',loc='8.2',set='era5',y=1960,id=17,iw=3,title=c("c)","d)"),reverse.prob=F)
dev.off()
