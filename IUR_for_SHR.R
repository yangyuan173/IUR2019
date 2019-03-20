####### Important Notes #####
# The data format should be .sas7bdat. Normmally the data are four-year data.
# The needed variables inlude 'shrty', 'h_admits', 'expectda', 'h_dy_yar', 'provfs', and 'year'.
# The data should be patient-level data. The facility-level SHR should be matched with patient-level data.
# Explanation of those variables can be found on Gillab C statistic document.
# If small facilities (PYA > 5) should be removed, this part should be done before calculating C statistic.
# Also, short facilities should be removed before running the code.
# String is case sensitive in R. The variable names should be exactly the same as in the data.
rm(list=ls())
uniqname = 'yuanyang'
output_path = 'K:/Users/kecc-yuanyang/C_statistics0214'
data_path = 'K:/Projects/Dialysis_Reports_Shared/Data/SMR/special_request'
raw_data_name = 'SMRSHR_2014to2017.sas7bdat'

obs_admission = "h_admits" #observed admission
exp_admission = "expectda" #expeceted admission
#hosp_yar = "h_dy_yar" #hospital YAR
facility = "provfs" #facility id
year = "year" #data year


############ Things you can change, but do not unless there are special needs. ###########
n.boot = 100
stratify.var = NULL
stratify.cut = NULL
measure.fun = cal_SMR
seed = 123

##########################################################################################
##				Do NOT make any edits below! 						    ##
##########################################################################################

# R_packages<-c("Rcpp","RcppArmadillo","git2r","sas7bdat")
# 
# sapply(R_packages,function(x){
#   if(!(x%in%rownames(installed.packages(lib.loc = "K:/Shared/xiaosonz/RPackages")))){
#     install.packages(x,lib="K:/Shared/xiaosonz/RPackages")
#   }
#   library(x,character.only=T,lib.loc = "K:/Shared/xiaosonz/RPackages")
# })

setwd(output_path)
# url.repo <- "git@gitlab.kecc.sph.umich.edu:R_PROJECTS_BIO/C_Statistic.git" # URL to Gitlab repo
# if (!dir.exists(paste0("C:/Users/kecc-",uniqname,"/.ssh"))) # check whether or not the folder containing SSH key pair exists
#   stop("Folder containing SSH key pair NOT found!") 
# # create credentials by means of SSH public/private key pair
# creds <- cred_ssh_key(publickey=paste0("C:/Users/kecc-",uniqname,"/.ssh/id_rsa.pub"), 
#                       privatekey=paste0("C:/Users/kecc-",uniqname,"/.ssh/id_rsa"))
# if (!dir.exists("./C_Statistic")) {
#   clone(url=url.repo,local_path="./C_Statistic",branch="master",credentials=creds)
# } else {
#   pull(repo=repository("./C_Statistic"),credentials=creds)
# }
# 
# source('./C_Statistic/C_statistic_R_functions.R')
# Rcpp::sourceCpp('./C_Statistic/C_stat_cpp.cpp')
# 
# if(dir.exists("./C_Statistic")){
#   if( Sys.info()["sysname"] == "Windows"){ 
#     shell(" rd /s /q C_Statistic")
#   }else{ 
#     system("rm -fr C_Statistic")
#   }
# }

#sink(file="IUR_SHR_log.txt",append=TRUE)
Sys.time()
cat("Reading data ... \n")
t0 = proc.time()[3]
if(file.exists(paste(data_path,"/SMRSHR.Rdata",sep=""))) {
  load(paste(data_path,"/SMRSHR.Rdata",sep=""))
}else{
  dat = read.sas7bdat(paste(data_path,"/",raw_data_name,sep="")) # takes about 8 minutes
  save(dat,file=paste(data_path,"/SMRSHR.Rdata",sep=""))
}
t1 = proc.time()[3]-t0
cat(paste("Reading data completed! Running time:",round(t1/60,digits = 2),"minutes.\n"))

cat("Checking absence of variables ...\n ")
var.list = c(obs_admission,exp_admission,hosp_yar,facility,year)
var.match = match(var.list,names(dat))
if(any(is.na(var.match))){
  stop(paste("Variable '", var.list[which(is.na(var.match))], "' NOT found!", sep=""),call.=F)
}
cat("Checking absence of variables completed!\n")

cat("Cleaning data ...\n ")
SHR_data0 = dat[,var.list]
n.all = nrow(SHR_data0)
n.0yar = sum(SHR_data0[,hosp_yar]==0)
SHR_data = SHR_data0[SHR_data0[,hosp_yar]>0,]
n.na = nrow(SHR_data)-sum(complete.cases(SHR_data))
SHR_data = SHR_data[complete.cases(SHR_data),]
cat(paste("Cleaning data completed!\n There are",n.all,"observations in the original data;\n"
          ,n.0yar,"were removed due to 0 hospital YAR;\n",n.na,"were removed due to NAs.\n"))

year.list = unique(SHR_data[,year])

for(yr in year.list){
  sub_index = (SHR_data[,year] == sub_index)
  IUR_temp = IUR_bootstrap(Obs = SHR_data[sub_index,obs_admission], Exp = SHR_data[sub_index,exp_admission], fac = SHR_data[sub_index,facility], n.boot = n.boot, stratify.var = stratify.var,stratify.cut = stratify.cut, measure.fun = cal_SMR, seed = seed)
  write.csv(IUR_temp,file = paste(output_path,"/IUR_",yr,".csv",sep=""))
}

IUR_temp = IUR_bootstrap(Obs = SHR_data[,obs_admission], Exp = SHR_data[,exp_admission], fac = SHR_data[,facility], n.boot = n.boot, stratify.var = stratify.var,stratify.cut = stratify.cut, measure.fun = cal_SMR, seed = seed)
write.csv(IUR_temp,file = paste(output_path,"/IUR_all.csv",sep=""))


