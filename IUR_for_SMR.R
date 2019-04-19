####### Important Notes #####
# The data format should be .sas7bdat. Normmally the data are four-year patient-level data.
# The must-have variables inlude 'dial_drd', 'expectda', 'provfs', 'DIAL_yar', and 'year'.
# The optional variable is stratify.var, which is used to calculate the sub-group IURs. If not provided, facility size will be used.
# Explanation of those variables can be found below or on Gillab IUR document.
# If small facilities (facility expected death < 5) should be removed, this part should be done before calculating IUR.
# Also, short facilities should be removed before running the code.
# String is case sensitive in R. The variable names should be exactly the same as in the data.
rm(list=ls())
uniqname = 'yuanyang'
output_path = 'K:/Users/kecc-yuanyang/IUR0321'
data_path = 'K:/Projects/Dialysis_Reports_Shared/Data/SMR/special_request'
raw_data_name = 'SMRSHR_2014to2017.sas7bdat'

source()
obs_death = "dial_drd" #observed death
exp_death = "expectda" #expeceted death
facility = "provfs" #facility id
death_yar = "DIAL_yar" #death year
year = "year" #data year

############ Things you can change, but do not unless there are special needs. ###########
n.boot = 100
stratify.var = NULL
stratify.cut = NULL
seed = 123

##########################################################################################
##				Do NOT make any edits below! 						    ##
##########################################################################################

R_packages<-c("git2r","sas7bdat")

sapply(R_packages,function(x){
  if(!(x%in%rownames(installed.packages(lib.loc = "K:/Shared/xiaosonz/RPackages")))){
    install.packages(x,lib="K:/Shared/xiaosonz/RPackages")
  }
  library(x,character.only=T,lib.loc = "K:/Shared/xiaosonz/RPackages")
})

library(git2r)
library(sas7bdat)

setwd(output_path)

url.repo <- "git@gitlab.kecc.sph.umich.edu:R_PROJECTS_BIO/IUR.git" # URL to Gitlab repo
if (!dir.exists(paste0("C:/Users/kecc-",uniqname,"/.ssh"))) # check whether or not the folder containing SSH key pair exists
  stop("Folder containing SSH key pair NOT found!")
# create credentials by means of SSH public/private key pair
creds <- cred_ssh_key(publickey=paste0("C:/Users/kecc-",uniqname,"/.ssh/id_rsa.pub"),
                      privatekey=paste0("C:/Users/kecc-",uniqname,"/.ssh/id_rsa"))
if (!dir.exists("./IUR")) {
  clone(url=url.repo,local_path="./IUR",branch="master",credentials=creds)
} else {
  pull(repo=repository("./IUR"),credentials=creds)
}

source('./IUR/IUR_R_functions.R')


if(dir.exists("./IUR")){
  if( Sys.info()["sysname"] == "Windows"){
    shell(" rd /s /q IUR")
  }else{
    system("rm -fr IUR")
  }
}

measure.fun = cal_SMR
sink(file="IUR_SMR_log.txt",append=TRUE)
cat('\n')
Sys.time()
cat("Reading data ... \n")
t0 = proc.time()[3]
if(file.exists(paste(data_path,"/SMRSHR.Rdata",sep=""))) {
  load(paste(data_path,"/SMRSHR.Rdata",sep=""))
}else{
  dat <- read.sas7bdat(paste(data_path,"/",raw_data_name,sep="")) # takes about 8 minutes
  save(dat,file=paste(data_path,"/SMRSHR.Rdata",sep=""))
}
t1 = proc.time()[3]-t0
cat(paste("Reading data completed! Running time:",round(t1/60,digits = 2),"minutes.\n"))

cat("Checking absence of variables ...\n ")
var.list = c(obs_death,exp_death,death_yar, facility,year)
var.match = match(var.list,names(dat))
if(any(is.na(var.match))){
  stop(paste("Variable '", var.list[which(is.na(var.match))], "' NOT found!", sep=""),call.=F)
}
cat("Checking absence of variables completed!\n")

cat("Cleaning data ...\n ")
SMR_data0 = dat[,var.list]
n.all = nrow(SMR_data0)
n.0yar = sum(SMR_data0[,death_yar]==0)
SMR_data = SMR_data0[SMR_data0[,death_yar]>0,]
n.na = nrow(SMR_data)-sum(complete.cases(SMR_data))
SMR_data = SMR_data[complete.cases(SMR_data),]
short_fac = SMR_data[,facility] =='Short'
SMR_data = SMR_data[!short_fac,]
cat(paste("Cleaning data completed!\n There are",n.all,"observations in the original data;\n"
          ,n.0yar,"were removed due to 0 death YAR;\n",n.na,"were removed due to NAs;\n"
          ,sum(short_fac),"were removed because their facility ids are 'Short'.\n"))
year.list = unique(SMR_data[,year])

cat("Calculating IURs ...\n ")
t0 = proc.time()
for(yr in year.list){
  sub_index = (SMR_data[,year] == yr)
  SMR_sub = SMR_data[sub_index,]
  SMR_sub[,facility] = factor(SMR_sub[,facility])
  SMR_sub = SMR_sub[order(SMR_sub[,facility]),]
  SMR_sub$expdeath = rep(sapply(split(SMR_sub[,exp_death],SMR_sub[,facility]),sum),
                          sapply(split(SMR_sub[,exp_death],SMR_sub[,facility]),length)
  )
  small_fac = SMR_sub[,"expdeath"] < 3
  SMR_sub = SMR_sub[!small_fac,]
  cat(paste("In",yr, "data,",sum(small_fac),"were removed because facility expected death < 3.\n"))
  
  IUR_temp = IUR_bootstrap(Obs = SMR_sub[,obs_death], Exp = SMR_sub[,exp_death], 
                           fac = SMR_sub[,facility], n.boot = n.boot, stratify.var = stratify.var,
                           stratify.cut = stratify.cut, measure.fun = cal_SMR, seed = seed)
  write.csv(IUR_temp$IUR,file = paste(output_path,"/IUR_SMR_",yr,"_",format(Sys.Date(),"%b%d%Y"),".csv",sep=""))
  write.csv(IUR_temp$IUR.fac,file = paste(output_path,"/facility_IUR_SMR_",yr,"_",format(Sys.Date(),"%b%d%Y"),".csv",sep=""))
}

t1 = proc.time()
(t1-t0)[3]
cat(paste("Yearly IUR calculation is completed! Average running time:",round(((t1-t0)[3])/60/4,digits = 2),"minutes.\n"))

SMR_data[,facility] = factor(SMR_data[,facility])
SMR_data = SMR_data[order(SMR_data[,facility]),]
SMR_data$expdeath = rep(sapply(split(SMR_data[,exp_death],SMR_data[,facility]),sum),
                        sapply(split(SMR_data[,exp_death],SMR_data[,facility]),length)
)
small_fac = SMR_data[,"expdeath"] < 3
SMR_data = SMR_data[!small_fac,]
cat(paste("In 4-year data,",sum(small_fac),"were removed because facility expected death < 3.\n"))

IUR_temp = IUR_bootstrap(Obs = SMR_data[,obs_death], Exp = SMR_data[,exp_death], fac = SMR_data[,facility], 
                         n.boot = n.boot, stratify.var = stratify.var,stratify.cut = stratify.cut, 
                         measure.fun = cal_SMR, seed = seed)
write.csv(IUR_temp$IUR,file = paste(output_path,"/IUR_SMR_all","_",format(Sys.Date(),"%b%d%Y"),".csv",sep=""))
write.csv(IUR_temp$IUR.fac,file = paste(output_path,"/facility_IUR_SMR_all","_",format(Sys.Date(),"%b%d%Y"),".csv",sep=""))
t2 = proc.time()
(t2-t1)[3]
cat(paste("Four-year IUR calculation is completed! Running time:",round(((t2-t1)[3])/60,digits = 2),"minutes.\n"))
sink()
