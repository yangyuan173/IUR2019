
IUR_core<-function(size,bootsdata,measureorg){
  #
  s2_star = apply(bootsdata,1,var) 
  s2_w = sum((size-1)*s2_star)/sum(size-1)
  nF=length(size)
  n_prime=(sum(size)-sum(size^2)/sum(size))/(nF-1)#nij=1
  T_mean=sum(size*measureorg)/sum(size)
  s2_t=sum(size*(measureorg-T_mean)^2)/(n_prime*(nF-1))
  s2_b=s2_t-s2_w
  IUR=s2_b/s2_t
  Fv=1/(1-IUR)
  p<-1-pf(Fv,nF-1,sum(size-1))
  results=c(s2_b,s2_w,IUR,p,nF,n_prime,T_mean,Fv)
  return(results)
}


IUR_boot <- function(data=data, q2=NULL,q3=NULL){
  time0=proc.time()[3]
  data=data[order(data$provfs),]
  m2 <- as.vector(sapply(split(data$ptyear,factor(data$provfs)),length)) 
  #m2 is the sample size for unique facility; length(m)=number of facility
  patient_year_pre <- as.vector(sapply(split(data$ptyear,factor(data$provfs)),sum)) 
  #m is the sample size for unique facility; length(m)=number of facility
  
  patient_year=rep(patient_year_pre,m2)
  facility=cbind(data$provfs,patient_year)
  
  facility=as.data.frame(facility)
  facilities<-facility[!duplicated(facility),]
  
  facilities=data.frame(facilities,size=m2)
  colnames(facilities)=c("facility","patient_year","size")
  
  obsadm<-as.vector(xtabs(numerator~factor(provfs),data))#obsadm[1:10]
  expadm<-as.vector(xtabs(denominator~factor(provfs),data))#expadm[1:10]
  data_summary<-cbind(facilities,obsadm,expadm)
  
  data_summary$measure<-data_summary$obsadm/data_summary$expadm
  
  data_pre_boot<-data
  
  m<-length(m2)
  mark<-cumsum(m2)
  total<-sum(m2)
  
  loop=0
  nloop=100
  repeat {
    loop=loop+1
    sample1<-NULL
    sample1<-c(sample1,sample.int(m2[1],size=m2[1],replace=TRUE))
    for(i in 2:m){ 
      sample1<-c(sample1,mark[i-1]+sample.int(m2[i],size=m2[i],replace=TRUE))  
    }
    data_pre_boot3<-data_pre_boot[sample1,]
    obsadm<-as.vector(xtabs(numerator~factor(provfs),  data_pre_boot3))
    expadm<-as.vector(xtabs(denominator~factor(provfs),  data_pre_boot3))
    
    databoot<-obsadm/expadm
    data_summary<-cbind(data_summary,databoot)
    if(loop==nloop) break 
  }
  
  #  modify criterion to choose working data
  data_restrict1<-data_summary[data_summary$patient_year>=5,]
  m3<-m2[data_summary$patient_year>=5]
  
  dataoriginal<-data_summary[,c(2:6)]
  nrow(dataoriginal)
  bootsdata<-data_restrict1[,-(1:6)]
  measureorg<-data_restrict1[,6]
  
  measureorg.summary<-c(length(measureorg),mean(measureorg),sd(measureorg),min(measureorg),max(measureorg))
  
  if(is.null(q2)){
    quantile_m<-quantile(data_restrict1$patient_year,prob=c(0.333,0.667))
    q2=quantile_m[1]
    q3=quantile_m[2]
  }
  
  label1<-cut(data_restrict1$size,breaks=c(-1,c(q2,q3),max(dataoriginal$size)))#STrR 1 year
  label2<-table(label1)
  label3<-names(label2)
  
  IUR_all<-matrix(0,nr=4,nc=2)
  IUR_all[1,]<-IUR_core(m3,bootsdata,measureorg)[c(3,5)]
  
  for(i in 1:3){
    IUR_result=IUR_core(m3[label1==label3[i]],bootsdata[label1==label3[i],],measureorg[label1==label3[i]])
    IUR_all[i+1,]=IUR_result[c(3,5)]
  }
  rownames(IUR_all)=c("Total","Group 1","Group 2","Group 3")
  colnames(IUR_all)=c("IUR","Group size")
  run.time=proc.time()[3]-time0
  return(list(IUR=IUR_all,cutoff=c(min(data_summary$patient_year),q2,q3,max(data_summary$patient_year)),run.time=run.time))
}

################## complete ################
stra_sampling<-function(fac){
  #sampling with replacement within each strata.
  #fac: a sorted factor vector, e.g. facility.
  n = length(fac)
  index = 1:n
  sample_index = rep(0,n)
  for(fac_temp in fac){
    sample_index[fac==fac_temp] = sample(index[fac==fac_temp],replace = TRUE)
  }
  return(sample_index)
}


IUR_bootdata <- function(size,measure.boot,measure.org){
  #         size: facility size vector. Each item corresponds to a facility size.
  # measure.boot: bootstraped measure vector.
  #  measure.org: original measure vector, e.g. SMR and SHR.
  s2_star = apply(measure.boot,1,var)      # within-facility variances
  s2_w = sum((size-1)*s2_star)/sum(size-1) # mean within-facility variance
  nF = length(size) # number of facilities.
  n_prime = (sum(size)-sum(size^2)/sum(size))/(nF-1)  
  T_mean = sum(size*measure.org)/sum(size) # mean orginal measure
  s2_t = sum(size*(measure.org-T_mean)^2)/(n_prime*(nF-1)) # total variance
  s2_b = s2_t-s2_w # between-facility variance
  IUR = s2_b/s2_t # total IUR
  IUR.fac = s2_b/(s2_b+s2_w/size) # facility level IUR
  return(list(IUR = IUR,nF = nF,IUR.fac = IUR.fac))
  #      IUR: total IUR
  #       nF: number of facilities
  # IUR.fac: faclity level IUR
}


cal_SMR <- function(Obs,Exp,fac){
  #Note: caculate SMR type measure. Data should be sorted by facitlity (fac).
  # Obs: observed outcomes;
  # Exp: expected outcomes;
  # fac: facility id vector;
  Obs.sum = as.vector(sapply(split(Obs,factor(fac)),sum))
  Exp.sum = as.vector(sapply(split(Exp,factor(fac)),sum))
  measure = Obs.sum/Exp.sum
  return(measure)
}


IUR_bootstrap<-function(Obs, Exp, fac, n.boot = 100, stratify.var = NULL,stratify.cut = NULL, measure.fun = cal_SMR, seed = 123){
  # Note: IUR using bootstrap method. The data should be sorted by facitlity (fac).
  #          Obs: observed outcomes;
  #          Exp: expected outcomes;
  #          fac: facility id vector;
  #       n.boot: the number of bootstraps;
  # stratify.var: stratification variable;
  # stratify.cut: stratification cutoff. It should be a vector with two cutoff points. If it is NULL, tertiles will be used.
  #  measure.fun: function to calculate the measure of interest. For now, we only implement SMR type of function.
  #         seed: an integer number to generate random numbers.
  set.seed(seed)
  # following 5 lines are to sort the data by fac. They could be removed if we sort the data before this function.
  Obs = Obs[order(fac)]
  Exp = Exp[order(fac)]
  if(!is.null(stratify.var)){
    stratify.var = stratify.var[order(fac)]
    fac = fac[order(fac)]
  }else{
    fac = fac[order(fac)]
    stratify.var = fac
  }
  
  size = as.vector(sapply(split(fac,factor(fac)),length))
  measure.origin = measure.fun(Obs,Exp,fac)
  measure.boot = matrix(NA,nrow = length(size),ncol = n.boot)
  
  # bootstraping
  loop = 0
  nloop = n.boot
  repeat {
    loop = loop+1
    sample.id = stra_sampling(fac)
    measure.boot[,loop] = measure.fun(Obs[sample.id],Exp[sample.id],fac[sample.id])
    if(loop == nloop) break
  }
  # calculating IUR
  fit.IUR = IUR_bootdata(size,measure.boot,measure.origin)
  IUR = c(fit.IUR$IUR,fit.IUR$nF)
  IUR.fac = fit.IUR$IUR.fac
  
  if(is.null(stratify.cut)){
    stra_sum = as.vector(sapply(split(stratify.var,factor(fac)),sum))
    stratify.cut = quantile(stra_sum,prob = c(0.33,0.67))
  }
  label1 = (stra_sum <= stratify.cut[1]) # right-closed interval
  label2 = (stratify.cut[1] < stra_sum & stra_sum <= stratify.cut[2])
  label3 = (stratify.cut[2] < stra_sum) 
  label = list(label1,label2,label3)
  for(i in 1:3){
    fit.IUR = IUR_bootdata(size[label[[i]]],measure.boot[label[[i]],],measure.origin[label[[i]]])
    IUR_group = c(fit.IUR$IUR,fit.IUR$nF)
    IUR  = rbind(IUR,IUR_group)
  }
  colnames(IUR) = c("IUR","Group size")
  return(list(IUR = IUR, IUR.fac = IUR.fac))
}


