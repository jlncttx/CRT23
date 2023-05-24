#
# (C) ebm_fit.R (Soulivanh) + ebm_response.R (Aurelien)
# 
# 15/04/2021


# Subroutine: EBM model function (from Soulivanh, notations Olivier)
hmodel <- function(FF, c, c0, lamb, gamm){
  N <- length(FF)
  dt <- 1; #-- timestep (year)
  #-1- Numerical solutions (explicit scheme)
  T <- numeric(N+1);
  To <- numeric(N+1);
  # H <- numeric(N+1);
  T[1] <- 0;
  To[1] <- 0;
  for(n in 2:(N+1)){
    T[n]  <-  (T[n-1] + dt/c*(FF[n-1] - lamb*T[n-1] - gamm*(T[n-1]-To[n-1])));
    To[n] <-  (To[n-1] + dt/c0*(gamm*(T[n-1]-To[n-1])));
    # H[n]  <-  gamm*(T[n] - To[n]);
  }
  # ans <- cbind(T, To, H)
  # colnames(ans) <- c("T", "To", "H")
  # ans
  T[-1]
}


# Main routine: EBM response as function of EBM params (adapted from Aurelien)
# Inputs:
#   - RF_time_series is a data.frame $year $RF_value
#   - ebm_params     is a data.frame with multiple sets of params F c c0 lamb gamm used in hmodel()
#   - period         is the period of interest (must be included in RF_time_series)
#   - nres           is a nb for randomly resampled responses (0 by default means no resampling)
# Output:
#   data.frame $year $be (best estimate) and then $res1..nres (if nres>0)

ebm_response = function(RF_time_series,ebm_params,period,nres=0){
  # Calculate response for each set of params
  all=c(); for (i in 1:nrow(ebm_params)){
    alli=hmodel(RF_time_series[,2],ebm_params$c[i],ebm_params$c0[i],
		  ebm_params$lamb[i],ebm_params$gamm[i])
    all=cbind(all,alli[RF_time_series$year %in% period])
  }
  # Best-estimate
  out=cbind(period,apply(all,1,mean))
  names=c("year","be")
  # Random resampling of responses
  if (nres>0){
    iparam=sample(1:nrow(ebm_params),nres,replace=T)
    out=cbind(out,all[,iparam])
    names=c(names,paste0("nres",1:nres))
  }
  # Output
  out=as.data.frame(out)
  names(out)=names
  out
}


#---
# GAM
library(gam)
forced_response_by_gam_fit=function(data,Enat,df){
  time = 1:length(data)
  fit = gam(data ~ s(time,df-1) + Enat)
#  fit = gam::gam(data ~ gam::s(time,df-1) + Enat)
  beta_nat = fit$coefficients[3]
  x_nat_ini = beta_nat*Enat
#  x_ant_ini = mysmoothy(data-x_nat_ini,df=df)#sdf=sdf)
  x_ant_ini = data - x_nat_ini - fit$residuals
  x_ant = x_ant_ini - x_ant_ini[1]
  x_all = mean(data) + x_nat_ini - mean(x_nat_ini) + x_ant - mean(x_ant)
  x_nat = x_all - x_ant
  data.frame(data=data,all=x_all,ant=x_ant,nat=x_nat,res=fit$residuals)
}

#---
# Fit a-la-Alix avec g connue (x(d,y) = f(d) + g(y).h(d) + eps)
fit_onto_forced_response=function(x,g,f.df,h.df){
  g=g-mean(g)
  fit=t(apply(matrix(x,nrow=365),1,function(x) lm(x~g)$coeff))
  sfit=cbind(mysmoothy(fit[,1],df=f.df,per=T),mysmoothy(fit[,2],df=h.df,per=T))
  rep(sfit[,1],length(g))+rep(sfit[,2],length(g))*rep(g,each=365)
}

