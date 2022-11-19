#Scripts for fitting all 3 models to anonymized data
library(TMB)
library(INLA)
library(sf)

# Logit function
logit=function(p){log(p/(1-p))}
# Inverse logist function
invlogit=function(t){exp(t)/(1+exp(t))}

model_dir<-"C:/Users/mcdonaldra/Documents/GitHub/SEHBAM-paper-scripts/src/"
other_dir<-"C:/Users/mcdonaldra/Documents/GitHub/SEHBAM-paper-scripts/"

dyn.unload(dynlib(paste0(model_dir,"SEHBAM")))
TMB::compile(paste0(model_dir,"SEHBAM.cpp"))
dyn.load(dynlib(paste0(model_dir,"SEHBAM")))

load(paste0(other_dir,"anon_data_sehbam.RData"))

params<-list()
params$log_sigma_epsilon<-0
params$log_sigma_upsilon<-0
params$log_R0<-log(200)
params$log_B0<-log(10000)
params$log_m0<-log(0.1)
params$log_qR<-rep(log(0.3),data$n_h)
params$beta_I<-rep(0,data$n_h)
# params$beta_IR<-rep(0,data$n_h)
params$mu_I<- 150
params$mu_IR<-rep(20,data$n_h)
params$log_var_I<-10
params$log_var_IR<-rep(10,data$n_h)

params$omega_B<-matrix(rep(0,data$n_m*((data$n_t+1))),ncol=(data$n_t+1))
params$omega_R<-matrix(rep(0,data$n_m*(data$n_t)),ncol=data$n_t)

Range_B<-60
Range_R<-60
Range_m<-60
Sigma_B<-0.3
Sigma_R<-0.3
Sigma_m<-0.3
params$log_kappa_B<-log(sqrt(8)/Range_B)
params$log_tau_B<-log(1/((Sigma_B)*sqrt(4*pi)*exp(params$log_kappa_B)))
params$log_kappa_R<-log(sqrt(8)/Range_R)
params$log_tau_R<-log(1/((Sigma_R)*sqrt(4*pi)*exp(params$log_kappa_R)))


#For anisotropy, simply use the ones obtained from model fit for simplicity
params$log_H_input_B<-c(0,0)
params$log_H_input_R<-c(0,0)

params$omega_m<-matrix(rep(0,data$n_m*((data$n_t+1))),ncol=(data$n_t+1))

params$log_kappa_m<-log(sqrt(8)/Range_m)
params$log_tau_m<-log(1/((Sigma_m)*sqrt(4*pi)*exp(params$log_kappa_m)))
params$log_H_input_m<-c(0,0)

params$log_qI<-rep(NA,data$n_s)
for (i in 1:length(params$log_qI)){
  if (data$h_s[i]<2) {
    params$log_qI[i]<-log(0.3)
  } else if (data$h_s[i]==2){
    params$log_qI[i]<-log(0.45)
  }
}

params$log_qR<-rep(log(0.3),data$n_h)

params$log_S<-rep(log(0.3),data$n_h)

random<-c("omega_B","omega_R","omega_m","log_qI")

maps<-list(log_H_input_R=c(factor(NA),factor(NA)))

non_r<-names(params)[-which(names(params) %in% random)]

obj <- MakeADFun( data=data, parameters=params,random=random, map = maps )
Opt<-optimx::optimr(obj$par,obj$fn,obj$gr,
                    control=list(maxit=5000000,maxfeval=5000000),
                    method="nlminb")
while (Opt$message=="iteration limit reached without convergence (10)") {
      obj$par<-obj$env$last.par[which(names(obj$env$last.par) %in% non_r)]
      Opt <- optimx::optimr(obj$par,obj$fn,obj$gr,control=list(maxit=5000000,maxfeval=5000000),method="nlminb")
    }
rep<-sdreport(obj)
Report<-obj$report()

save(obj,Opt,rep,Report,file=paste0(other_dir,"manuscript_SEHBAM_fit2.RData"))

load(paste0(other_dir,"anon_data_SEAM.RData"))

params<-list()
params$log_sigma_epsilon<-0
params$log_sigma_upsilon<-0
params$log_R0<-log(500)
params$log_B0<-log(10000)
params$log_m0<-log(0.1)
params$log_qI<-rep(log(0.3),data$n_s)
params$log_qR<-log(0.3)
params$logit_p_I<-0
params$logit_p_IR<-0

params$omega_B<-matrix(rep(0,data$n_m*((data$n_t+1))),ncol=(data$n_t+1))
params$omega_R<-matrix(rep(0,data$n_m*(data$n_t)),ncol=data$n_t)

Range_B<-40
Range_R<-60
Range_m<-20
Sigma_B<-0.2
Sigma_R<-0.2
Sigma_m<-0.1
params$log_kappa_B<-log(sqrt(8)/Range_B)
params$log_tau_B<-log(1/((Sigma_B)*sqrt(4*pi)*exp(params$log_kappa_B)))
params$log_kappa_R<-log(sqrt(8)/Range_R)
params$log_tau_R<-log(1/((Sigma_R)*sqrt(4*pi)*exp(params$log_kappa_R)))


#For anisotropy, simply use the ones obtained from model fit for simplicity
params$log_H_input_B<-c(0,0)
params$log_H_input_R<-c(0,0)

params$omega_m<-matrix(rep(0,data$n_m*((data$n_t+1))),ncol=(data$n_t+1))

params$log_kappa_m<-log(sqrt(8)/Range_m)
params$log_tau_m<-log(1/((Sigma_m)*sqrt(4*pi)*exp(params$log_kappa_m)))
params$log_H_input_m<-c(0,0)

params$log_S<-log(0.4)

random<-c("omega_B","omega_R","omega_m")

maps<-list(log_H_input_R=c(factor(NA),factor(NA)))

non_r<-names(params)[-which(names(params) %in% random)]

obj <- MakeADFun( data=data, parameters=params,random=random, map = maps )
Opt<-optimx::optimr(obj$par,obj$fn,obj$gr,
                    control=list(maxit=5000000,maxfeval=5000000),
                    method="nlminb")
while (Opt$message=="iteration limit reached without convergence (10)") {
  obj$par<-obj$env$last.par[which(names(obj$env$last.par) %in% non_r)]
  Opt <- optimx::optimr(obj$par,obj$fn,obj$gr,control=list(maxit=5000000,maxfeval=5000000),method="nlminb")
}
rep<-sdreport(obj)
Report<-obj$report()

save(obj,Opt,rep,Report,file=paste0(other_dir,"SEBDAM_BF_fit_60_knots_habitat_prior_converged2.RData"))

load(paste0(other_dir,"anon_data_TLM.RData"))

params<-list()
params$log_sigma_tau<--1
params$log_sigma_phi<--1
params$log_sigma_m<--1
params$log_sigma_epsilon<--1
params$log_sigma_upsilon<--1
params$log_q_I<--1
params$log_q_R<--1
params$logit_p_I<-0
params$logit_p_IR<-0
params$log_S<-0


params$log_B<-rep(log(10000),length(unique(data$t_i))+1)
params$log_R<-rep(log(2000),length(unique(data$t_i))+1)
params$log_input_m<-rep(log(0.1),length(unique(data$t_i))+1)

random<-c("log_B","log_R","log_input_m")

maps<-list()

non_r<-names(params)[-which(names(params) %in% random)]

obj <- MakeADFun( data=data, parameters=params,random=random, map = maps )
Opt<-optimx::optimr(obj$par,obj$fn,obj$gr,
                    control=list(maxit=5000000,maxfeval=5000000),
                    method="nlminb")
while (Opt$message=="iteration limit reached without convergence (10)") {
  obj$par<-obj$env$last.par[which(names(obj$env$last.par) %in% non_r)]
  Opt <- optimx::optimr(obj$par,obj$fn,obj$gr,control=list(maxit=5000000,maxfeval=5000000),method="nlminb")
}
rep<-sdreport(obj)
Report<-obj$report()

save(obj,Opt,rep,Report,file=paste0(other_dir,"TLM_BF_fit_manuscript_diff_prior2.RData"))

