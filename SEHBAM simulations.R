library(INLA)
library(TMB)

model_dir<-"C:/Users/mcdonaldra/Documents/GitHub/SEHBAM-paper-scripts/src/"
fig_dir<-"C:/Users/mcdonaldra/Documents/GitHub/SEHBAM-paper-scripts/"


dyn.unload(dynlib(paste0(model_dir,"SEHBAM")))
TMB::compile(paste0(model_dir,"SEHBAM.cpp"))
dyn.load(dynlib(paste0(model_dir,"SEHBAM")))

library(sp)
library(dplyr)
library(raster)
library(ggplot2)

# Logit function
logit=function(p){log(p/(1-p))}
# Inverse logist function
invlogit=function(t){exp(t)/(1+exp(t))}

#Testing sehbam simulation?

nknot<-40

#Create a 100X100km square
x_coord<-c(150,250,250,150)
y_coord<-c(150,150,250,250)
xy<-cbind(x_coord,y_coord)
p = Polygon(xy)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))
# plot(sps)

#Randomly sample number of locations for tow locations
NY<-15
yearly_tows<-200
set.seed(30)
random_loc<-spsample(sps,NY*yearly_tows,"random")
#Setting up knots based on these locations
# nknot<-25
# set.seed(26)
set.seed(110) #for 40 years
knots<-kmeans(cbind(random_loc@coords[,1],random_loc@coords[,2]),nknot,nstart=25)
knots.loc<-as.data.frame(knots[[2]])
knotID<-knots$cluster

#Create the mesh
utm.prj4s <- CRS("+init=epsg:32619")
sp_knots<-SpatialPoints(knots.loc,proj4string = utm.prj4s)
mesh<- inla.mesh.2d(sp_knots,
                    max.edge=c(10,30),cutoff=2, #original is 30 and 1.5
                    boundary=inla.sp2segment(sps))
#Create anisotropic spde object for model
spde<-inla.spde2.matern(mesh)
# Triangle info
Dset = 1:2
TV = mesh$graph$tv           # Triangle to vertex indexing
V0 = mesh$loc[TV[,1],Dset]   # V = vertices for each triangle
V1 = mesh$loc[TV[,2],Dset]
V2 = mesh$loc[TV[,3],Dset]
E0 = V2 - V1                      # E = edge for each triangle
E1 = V0 - V2
E2 = V1 - V0
# Calculate Areas
TmpFn = function(Vec1, Vec2) abs(det( rbind(Vec1, Vec2) ))
Tri_Area = rep(NA, nrow(E0))
for(i in 1:length(Tri_Area)) Tri_Area[i] = TmpFn( E0[i,],E1[i,] )/2   # T = area of each triangle
spde_aniso <- list(
  "n_s"      = spde$n.spde,
  "n_tri"    = nrow(TV),
  "Tri_Area" = Tri_Area,
  "E0"       = E0,
  "E1"       = E1,
  "E2"       = E2,
  "TV"       = TV - 1,
  "G0"       = spde$param.inla$M0,
  "G0_inv"   = as(diag(1/diag(spde$param.inla$M0)), "dgTMatrix"))

#To create area, must create grid and attribute each cell to closest knot
grid<-raster(extent(sps),resolution=c(1,1))
grid<-raster::extend(grid,c(1,1))
gridpoly<-rasterToPolygons(grid)
clipgrid<-raster::intersect(gridpoly,sps)

centlon<-c()
centlat<-c()
centsize<-c()
for (i in 1:length(clipgrid@polygons)) {
  centlon<-c(centlon,clipgrid@polygons[[i]]@labpt[1])
  centlat<-c(centlat,clipgrid@polygons[[i]]@labpt[2])
  centsize<-c(centsize,clipgrid@polygons[[i]]@area)
}
centroids<-cbind(centlon,centlat)
distmat<-pointDistance(centroids,sp_knots@coords,lonlat=F)
polyknotID<-c()
for (i in 1:length(clipgrid@polygons)) {
  polyknotID<-c(polyknotID,which(distmat[i,]==min(distmat[i,])))
}
centroids<-cbind(centroids,polyknotID,centsize,c(as.character(1:length(clipgrid@polygons))))
centroids<-as.data.frame(centroids)
centroids$knotID<-as.vector(centroids$knotID)
colnames(centroids)<-c("lon","lat","knotID","area","id")
centroidsID<-centroids[,c(3,5)]
centroids$area<-as.numeric(as.character(centroids$area))
stratarea<-aggregate(area~knotID,sum,data=centroids)

#Setup for plotting purposes
clipfort<-fortify(clipgrid)
plotloc<-as.data.frame(sp_knots@coords)
plotloc$group<-rep(1,length(plotloc[,1]))
plotloc$knotID<-as.factor(1:nknot)
plotloc<-left_join(plotloc,stratarea,by="knotID")
stratplot<-left_join(clipfort,centroidsID,by="id")

library(sf)
sf_knots<-st_as_sf(sp_knots)
st_crs(sf_knots)<-NA
sf_loc<-st_as_sf(random_loc)
distmat<-st_distance(sf_loc,sf_knots)
polyknotID<-rep(NA,nrow(distmat))
for (i in 1:nrow(distmat)){
  polyknotID[i]<-which(distmat[i,]==min(distmat[i,]))
}
polyknotID<-polyknotID-1

n_hab<-3

data<-list()

data$model<-"SEHBAM"

data$options_vec<-c(2,1,0,1,1,1,0)

data$logI<-rep(1,NY*yearly_tows)
data$logIR<-rep(1,NY*yearly_tows)
data$area<-stratarea$area
data$C<-matrix(rep(1,(NY+1)*nknot),ncol=(NY+1),nrow=nknot)
data$n_I<-rep(1,NY*yearly_tows)
data$n_IR<-rep(1,NY*yearly_tows)

data$n_i<-length(data$logI)
data$n_t<-NY
data$n_s<-nknot
data$n_m<-mesh$n
data$n_h<-n_hab

#If even sampling
# data$s_i<-rep(rep((1:nknot)-1,yearly_tows/nknot),NY)
data$s_i<-polyknotID
#If random, so doesn't have to be the same effort everywhere
# data$s_i<-sample(1:nknot,length(data$I),replace=T)
data$t_i<-rep(0:(NY-1),each=yearly_tows)
data$v_i<-mesh$idx$loc-1
data$h_i<-rep(NA,length(data$s_i))
# for (i in 1:length(data$logI)) {
#   if (data$s_i[i] %in% c(0:4)) data$h_i[i]<-0
#   if (data$s_i[i] %in% c(5:9)) data$h_i[i]<-1
#   if (data$s_i[i] %in% c(10:14)) data$h_i[i]<-2
#   if (data$s_i[i] %in% c(15:19)) data$h_i[i]<-3
#   if (data$s_i[i] %in% c(20:24)) data$h_i[i]<-4
# }
for (i in 1:length(data$logI)) {
  if (data$s_i[i] %in% c(0:13)) data$h_i[i]<-0
  if (data$s_i[i] %in% c(14:26)) data$h_i[i]<-1
  if (data$s_i[i] %in% c(27:39)) data$h_i[i]<-2
}

data$X<-matrix(nrow=data$n_i,ncol=data$n_h)
# for (i in 1:nrow(data$X)){
#   if (data$h_i[i]==0) data$X[i,]<-c(1,0,0,0,0)
#   if (data$h_i[i]==1) data$X[i,]<-c(1,1,0,0,0)
#   if (data$h_i[i]==2) data$X[i,]<-c(1,0,1,0,0)
#   if (data$h_i[i]==3) data$X[i,]<-c(1,0,0,1,0)
#   if (data$h_i[i]==4) data$X[i,]<-c(1,0,0,0,1)
# }
for (i in 1:nrow(data$X)){
  if (data$h_i[i]==0) data$X[i,]<-c(1,0,0)
  if (data$h_i[i]==1) data$X[i,]<-c(1,1,0)
  if (data$h_i[i]==2) data$X[i,]<-c(1,0,1)
}

data$s_a<-rep(0,each=data$n_s)

#Growth rates
set.seed(122)
data$gI<-as.matrix(rbind(sample(seq(0.8,1.2,by=0.05),NY,replace=T)))
data$gR<-as.matrix(rbind(sample(seq(0.8,1.2,by=0.05),NY,replace=T)))

data$h_s<-c(rep(0:2,each=nknot/3),2)

data$plug_exploit<-0.1

data$mesh_obj<-spde_aniso

data$L<-rep(1,NY*yearly_tows)
data$n_bin<-rep(1,NY*yearly_tows)
data$prior_pars<-c(50,130,50,130,50,70)
data$prior_pars_1<-rep(NA,nknot)
data$prior_pars_2<-rep(NA,nknot)
for (i in 1:nknot){
  if (data$h_s[i]==0) {data$prior_pars_1[i]<-50; data$prior_pars_2[i]<-130}
  if (data$h_s[i]==1) {data$prior_pars_1[i]<-50; data$prior_pars_2[i]<-130}
  if (data$h_s[i]==2) {data$prior_pars_1[i]<-50; data$prior_pars_2[i]<-70}
}

params<-list()
params$log_sigma_epsilon<-log(1)
params$log_sigma_upsilon<-log(1)
params$log_R0<-6
params$log_B0<-7
params$log_m0<-log(0.1)
params$log_qR<-c(log(0.1),log(0.15),log(0.05))
params$beta_I<-c(-3.14, -0.21, 1.92)
params$mu_I<- 192
params$mu_IR<-c(40,15,7.5)
params$log_var_I<-10.8
params$log_var_IR<-c(9.33,6.94,6.29)

params$omega_B<-matrix(rep(0,mesh$n*((data$n_t+1))),ncol=(data$n_t+1))
params$omega_R<-matrix(rep(0,mesh$n*(data$n_t)),ncol=data$n_t)

Range_B<-75
Range_R<-59
Range_m<-15
Sigma_B<-0.3
Sigma_R<-0.5
Sigma_m<-0.1
params$log_kappa_B<-log(sqrt(8)/Range_B)
params$log_tau_B<-log(1/((Sigma_B)*sqrt(4*pi)*exp(params$log_kappa_B)))
params$log_kappa_R<-log(sqrt(8)/Range_R)
params$log_tau_R<-log(1/((Sigma_R)*sqrt(4*pi)*exp(params$log_kappa_R)))


#For anisotropy, simply use the ones obtained from model fit for simplicity
params$log_H_input_B<-c(0.63,0.96)
params$log_H_input_R<-c(0,0.02)

params$omega_m<-matrix(rep(0,mesh$n*((data$n_t+1))),ncol=(data$n_t+1))

params$log_kappa_m<-log(sqrt(8)/Range_m)
params$log_tau_m<-log(1/((Sigma_m)*sqrt(4*pi)*exp(params$log_kappa_m)))
params$log_H_input_m<-c(-0.5,0.5)

params$log_qI<-rep(NA,data$n_s)
# for (i in 1:length(params$log_qI)){
#   if (data$h_s[i]==0) params$log_qI[i]<-log(rbeta(1,data$prior_pars[1],data$prior_pars[2]))
#   if (data$h_s[i]==1) params$log_qI[i]<-log(rbeta(1,data$prior_pars[3],data$prior_pars[4]))
#   if (data$h_s[i]==2) params$log_qI[i]<-log(rbeta(1,data$prior_pars[5],data$prior_pars[6]))
#   if (data$h_s[i]==3) params$log_qI[i]<-log(rbeta(1,data$prior_pars[7],data$prior_pars[8]))
#   if (data$h_s[i]==4) params$log_qI[i]<-log(rbeta(1,data$prior_pars[9],data$prior_pars[10]))
# }
set.seed(122)
for (i in 1:length(params$log_qI)){
  if (data$h_s[i]==0) params$log_qI[i]<-log(rbeta(1,data$prior_pars[1],data$prior_pars[2]))
  if (data$h_s[i]==1) params$log_qI[i]<-log(rbeta(1,data$prior_pars[3],data$prior_pars[4]))
  if (data$h_s[i]==2) params$log_qI[i]<-log(rbeta(1,data$prior_pars[5],data$prior_pars[6]))
}

set.seed(122)
params$log_S<-c(log(0.75),log(0.6),log(0.5))

random<-c("omega_B","omega_R","omega_m","log_qI")

maps<-list(log_H_input_R=c(factor(NA),factor(NA)))

obj <- MakeADFun( data=data, parameters=params,random=random, map = maps , DLL="SEHBAM")
test_dat<-obj$simulate(complete=T)

fit_params<-list()
fit_params$log_sigma_epsilon<--1
fit_params$log_sigma_upsilon<--1
fit_params$log_R0<-log(300)
fit_params$log_B0<-log(5000)
fit_params$log_m0<-log(0.1)
fit_params$log_qR<-rep(-1,data$n_h)
fit_params$beta_I<-rep(1,data$n_h)
fit_params$mu_I<- 200
fit_params$mu_IR<-rep(20,data$n_h)
fit_params$log_var_I<-6
fit_params$log_var_IR<-rep(5,data$n_h)

fit_params$omega_B<-matrix(rep(0,mesh$n*((data$n_t+1))),ncol=(data$n_t+1))
fit_params$omega_R<-matrix(rep(0,mesh$n*(data$n_t)),ncol=data$n_t)

fit_params$log_kappa_B<--1
fit_params$log_tau_B<-1
fit_params$log_kappa_R<--1
fit_params$log_tau_R<-1


#For anisotropy, simply use the ones obtained from model fit for simplicity
fit_params$log_H_input_B<-c(0,0)
fit_params$log_H_input_R<-c(0,0)

fit_params$omega_m<-matrix(rep(0,mesh$n*((data$n_t+1))),ncol=(data$n_t+1))

fit_params$log_kappa_m<--1
fit_params$log_tau_m<-1
fit_params$log_H_input_m<-c(0,0)

fit_params$log_qI<-rep(log(0.3),data$n_s)
fit_params$log_S<-rep(-1,data$n_h)

non_r<-names(params)[-which(names(params) %in% random)]

nrep<-200
Report_list<-list()
sim_save<-list()
diff_list<-list()
diff_list2<-list()
rep_list<-list()
msg_list<-list()
for (i in 1:nrep) {
  tryCatch({
    simdata <- obj $ simulate(complete=T)
    sim_save[[i]]<-simdata
    fit_params$log_B0<-log(max(exp(simdata$logI),na.rm=T))+1
    obj2 <- MakeADFun ( simdata , fit_params , random=random, map=maps )
    Opt2 <- try(optimx::optimr(obj2$par,obj2$fn,obj2$gr,
                               control=list(maxit=100000,maxeval=100000),method="nlminb"),T)
    while (Opt2$message=="iteration limit reached without convergence (10)") {
      obj2$par<-obj2$env$last.par[which(names(obj2$env$last.par) %in% non_r)]
      Opt2 <- try(optimx::optimr(obj2$par,obj2$fn,obj2$gr,
                                 control=list(maxit=100000,maxeval=100000),method="nlminb"),T)
    }
    rep_list[[i]] <- sdreport(obj2,bias.correct=F)
    Report_list[[i]] = obj2$report()
    if (!is.null(rep_list[[i]])){
      msg_list[[i]] <- Opt2$message
    }
    if ((i %% 5)==0){save(rep_list,sim_save,msg_list,Report_list,file=paste0("interrupt_sehbam_",i,"_good.RData"))}
  }, error=function(e) {})
}

#For setting 1
params$log_qI<-rep(NA,data$n_s)
set.seed(122)
for (i in 1:length(params$log_qI)){
  if (data$h_s[i]==0) params$log_qI[i]<-log(rbeta(1,data$prior_pars[1],data$prior_pars[2]))
  if (data$h_s[i]==1) params$log_qI[i]<-log(rbeta(1,data$prior_pars[3],data$prior_pars[4]))
  if (data$h_s[i]==2) params$log_qI[i]<-log(rbeta(1,data$prior_pars[5],data$prior_pars[6]))
}

converge<-0
false<-0
x_conv<-0
singular<-0
null_vec<-nrep
no_conv<-0
chosen_ones<-c()
for (i in 1:nrep){
  if(is.null(msg_list[[i]])) no_conv<-no_conv+1
  if (!is.null(msg_list[[i]])){null_vec<-null_vec-1
  if (msg_list[[i]]=="relative convergence (4)") {converge<-converge+1
  chosen_ones<-c(chosen_ones,i)}
  if (msg_list[[i]]=="false convergence (8)") {false<-false+1}
  if(msg_list[[i]]=="singular convergence (7)"){singular<-singular+1}
  if(msg_list[[i]]=="both X-convergence and relative convergence (5)"){x_conv<-x_conv+1
  chosen_ones<-c(chosen_ones,i)}
  if(msg_list[[i]]=="nlminb failed"){no_conv=no_conv+1}
  }
}
conv_frame<-data.frame(converge=converge,false=false,singular=singular,x=x_conv,no_conv=no_conv)

qI_val_sim<-ggplot()+geom_point(aes(x=c(7,21,33),y=c(0.24491,0.28538,0.44075),
                                    col=as.factor(c(0,1,2))),shape=17,cex=3)+
  geom_point(aes(x=1:40,y=exp(params$log_qI),
                 col=as.factor(data$h_s)),shape=16)+
  theme_bw()+
  ylab("Value")+xlab("Knot")+
  scale_color_discrete(name="Habitat")+
  guides(color = guide_legend(override.aes = list(shape = 16,cex=1.5) ) )
ggsave(filename=paste0(fig_dir,"qI_sim_val.png"),plot=qI_val_sim,
       width=6,height=5)

g_val_sim<-ggplot()+
  geom_point(aes(x=1:NY,y=data$gI[1,],
                 col="Commercial Size Growth",
                 shape="Commercial Size Growth"))+
  geom_line(aes(x=1:NY,y=data$gI[1,],
                col="Commercial Size Growth"))+
  geom_point(aes(x=1:NY,y=data$gR[1,],
                 col="Recruit Growth",
                 shape="Recruit Growth"))+
  geom_line(aes(x=1:NY,y=data$gR[1,],
                col="Recruit Growth"))+
  theme_bw()+
  scale_color_discrete(name="")+
  scale_shape_discrete(name="")+
  ylab("Value")+xlab("Year")
ggsave(filename=paste0(fig_dir,"growth_sim_val.png"),plot=g_val_sim,
       width=6,height=5)

sigma_epsilon<-c(rep(NA,nrep))
sigma_upsilon<-c(rep(NA,nrep))
S<-matrix(ncol=n_hab,nrow=nrep)
qI<-matrix(rep(NA,nrep*nknot),ncol=nknot)
qR<-matrix(ncol=n_hab,nrow=nrep)
mu_I<-c(rep(NA,nrep))
mu_IR<-matrix(ncol=n_hab,nrow=nrep)
log_var_I<-c(rep(NA,nrep))
var_IR<-matrix(ncol=n_hab,nrow=nrep)
B0<-c(rep(NA,nrep))
R0<-c(rep(NA,nrep))
m0<-c(rep(NA,nrep))
kappa_B<-c(rep(NA,nrep))
kappa_R<-c(rep(NA,nrep))
kappa_m<-c(rep(NA,nrep))
tau_B<-c(rep(NA,nrep))
tau_m<-c(rep(NA,nrep))
tau_R<-c(rep(NA,nrep))
log_H_input_B1<-c(rep(NA,nrep))
log_H_input_B2<-c(rep(NA,nrep))
log_H_input_m1<-c(rep(NA,nrep))
log_H_input_m2<-c(rep(NA,nrep))
SigmaO_B<-rep(NA,nrep)
SigmaO_R<-rep(NA,nrep)
SigmaO_m<-rep(NA,nrep)
for (i in chosen_ones){
  sigma_epsilon[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="sigma_epsilon")]
  sigma_upsilon[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="sigma_upsilon")]
  S[i,]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="S")]
  qR[i,]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="qR")]
  mu_I[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="mu_I")]
  mu_IR[i,]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="mu_IR")]
  log_var_I[i]<-log(rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="var_I")])
  var_IR[i,]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="var_IR")]
  B0[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="B0")]
  R0[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="R0")]
  m0[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="m0")]
  kappa_B[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="kappa_B")]
  kappa_R[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="kappa_R")]
  kappa_m[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="kappa_m")]
  tau_B[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="tau_B")]
  tau_R[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="tau_R")]
  tau_m[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="tau_m")]
  log_H_input_B1[i]<-Report_list[[i]]$log_H_input_B[1]
  log_H_input_B2[i]<-Report_list[[i]]$log_H_input_B[2]
  log_H_input_m1[i]<-Report_list[[i]]$log_H_input_m[1]
  log_H_input_m2[i]<-Report_list[[i]]$log_H_input_m[2]
  SigmaO_B[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="SigmaO_B")]
  SigmaO_R[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="SigmaO_R")]
  SigmaO_m[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="SigmaO_m")]
  # }
}

for (i in chosen_ones){
  for (j in 1:nknot){
    qI[i,j]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="qI")][j]
  }
}

qI_diffs<-matrix(rep(NA,nrep*nknot),ncol=nknot)
qR_diffs<-matrix(ncol=n_hab,nrow=nrep)
mu_IR_diffs<-matrix(ncol=n_hab,nrow=nrep)
log_var_IR_diffs<-matrix(ncol=n_hab,nrow=nrep)
S_diffs<-matrix(ncol=n_hab,nrow=nrep)
for (i in chosen_ones){
  qI_diffs[i,]<-qI[i,]-exp(sim_save[[i]]$log_qI)
  qR_diffs[i,]<-qR[i,]-sim_save[[i]]$qR
  mu_IR_diffs[i,]<-mu_IR[i,]-params$mu_IR
  log_var_IR_diffs[i,]<-log(var_IR[i,])-params$log_var_IR
  S_diffs[i,]<-S[i,]-sim_save[[i]]$S
}

all_params<-data.frame(sigma_epsilon,sigma_upsilon,
                       mu_I,log_var_I,
                       B0,R0,m0,kappa_B,
                       kappa_R,kappa_m,tau_B,tau_m,tau_R,
                       log_H_input_B1,log_H_input_B2,
                       log_H_input_m1,log_H_input_m2,
                       SigmaO_B,SigmaO_R,SigmaO_m)

hist_plot<-data.frame(value=c(sigma_epsilon,sigma_upsilon,
                              mu_I,log_var_I,B0,R0,m0,
                              kappa_B,kappa_R,kappa_m,
                              tau_B,tau_m,tau_R,
                              log_H_input_B1,log_H_input_B2,
                              log_H_input_m1,log_H_input_m2,
                              SigmaO_B,SigmaO_R,SigmaO_m),
                      parameter=rep(c("sigma_epsilon","sigma_upsilon",
                                      "mu_I","log_var_I","B0",
                                      "R0","m0",
                                      "kappa_B",
                                      "kappa_R",
                                      "kappa_m","tau_B",
                                      "tau_m","tauR",
                                      "H_input_B1","H_input_B2",
                                      "H_input_m1","H_input_m2",
                                      "SigmaO_B","SigmaO_R","SigmaO_m")
                                    ,each=nrep),
                      true=rep(c(exp(params$log_sigma_epsilon),
                                 exp(params$log_sigma_upsilon),params$mu_I,
                                 params$log_var_I,
                                 exp(params$log_B0),
                                 exp(params$log_R0),exp(params$log_m0),
                                 exp(params$log_kappa_B),exp(params$log_kappa_R),
                                 exp(params$log_kappa_m),
                                 exp(params$log_tau_B),exp(params$log_tau_m),exp(params$log_tau_R),
                                 params$log_H_input_B[1],params$log_H_input_B[2],
                                 params$log_H_input_m[1],params$log_H_input_m[2],
                                 1/((exp(params$log_tau_B))*sqrt(4*pi)*exp(params$log_kappa_B)),
                                 1/((exp(params$log_tau_R))*sqrt(4*pi)*exp(params$log_kappa_R)),
                                 1/((exp(params$log_tau_m))*sqrt(4*pi)*exp(params$log_kappa_m))),each=nrep))

hist_plot$parameter<-factor(hist_plot$parameter,
                            labels=c(expression(B [0]),expression("H" ["input"]^"B1"),
                                     expression("H" ["input"]^"B2"),
                                     expression("H" ["input"]^"m1"),expression("H" ["input"]^"m2"),
                                     expression(kappa ["B"]),expression(kappa ["m"]),
                                     expression(kappa ["R"]),expression(log (sigma["I"]^2)),
                                     expression("m" [0]),
                                     expression(mu["I"]),
                                     expression("R" [0]),
                                     expression(sigma [epsilon]),
                                     expression(sigma[upsilon]),
                                     expression(sigma[B]^Omega),
                                     expression(sigma[m]^Omega),
                                     expression(sigma[R]^Omega),
                                     expression(tau ["B"]),
                                     expression(tau ["m"]),expression(tau["R"])))

par_hist_set_1<-ggplot(hist_plot)+geom_histogram(aes(x=value),col="black",fill="grey")+
  geom_vline(aes(xintercept = true,col="True Value"))+
  theme(strip.background = element_blank(),legend.position = "none")+
  facet_wrap(~parameter,scales="free",labeller=label_parsed)+
  theme_bw()+
  xlab("Estimated Value")+ylab("Frequency")
ggsave(filename=paste0(fig_dir,"par_hist_1.png"),plot=par_hist_set_1,
       width=9,height=8)

plot_qI<-data.frame(value=as.vector(qI_diffs),
                    param=rep("qI",nknot*nrep))
plot_qR<-data.frame(value=as.vector(qR_diffs),
                    param=rep("qR",n_hab*nrep))
plot_mu_IR<-data.frame(value=as.vector(mu_IR_diffs),
                       param=rep("mu_IR",n_hab*nrep))
plot_log_var_IR<-data.frame(value=as.vector(log_var_IR_diffs),
                            param=rep("log_var_IR",n_hab*nrep))
plot_S<-data.frame(value=as.vector(S_diffs),
                   param=rep("S",n_hab*nrep))

diff_hist_plot<-rbind(plot_qI,plot_qR,plot_mu_IR,plot_log_var_IR,plot_S)

diff_hist<-ggplot(diff_hist_plot)+geom_histogram(aes(x=value),col="black",fill="grey")+
  geom_vline(aes(xintercept = 0),lty="dashed")+
  theme(strip.background = element_blank(),legend.position = "none")+
  facet_wrap(~param,scales="free",labeller=label_parsed)+
  theme_bw()+
  xlab("Estimated Value")+ylab("Frequency")
ggsave(filename=paste0(fig_dir,"diff_hist_1.png"),plot=diff_hist,
       width=7,height=6)

qI_by_knot<-ggplot(plot_qI2)+geom_histogram(aes(x=value),col="black",fill="grey")+
  geom_vline(aes(xintercept = 0),lty="dashed")+
  theme(strip.background = element_blank(),legend.position = "none")+
  facet_wrap(~knot,scales="free")+
  theme_bw()+
  xlab("Estimated Value")+ylab("Frequency")
ggsave(filename=paste0(fig_dir,"qI_by_knot.png"),plot=qI_by_knot,
       width=10,height=9)
#Knots to look more closely at are 4, 7, 12, 13, 14, 16,19,22,25,28,33 to 36,40

qI_by_hab<-ggplot(plot_qI2)+geom_histogram(aes(x=value),col="black",fill="grey")+
  geom_vline(aes(xintercept = 0),lty="dashed")+
  theme(strip.background = element_blank(),legend.position = "none")+
  facet_wrap(~hab,scales="free")+
  theme_bw()+
  xlab("Estimated Value")+ylab("Frequency")
ggsave(filename=paste0(fig_dir,"qI_by_hab_1.png"),plot=qI_by_hab,
       width=6,height=4)

C_avg<-matrix(nrow=nrep,ncol=nknot)
C_max<-matrix(nrow=nrep,ncol=nknot)
for (i in chosen_ones){
  C_avg[i,]<-rowMeans(sim_save[[i]]$C[,-(NY+1)])
  C_max[i,]<-max(sim_save[[i]]$C[,-(NY+1)])
}

check_qIdiffs_C<-data.frame(C=as.vector(C_avg),max=as.vector(C_max),qI_diff=as.vector(qI_diffs))
check_qI_C<-data.frame(C=as.vector(C_avg),max=as.vector(C_max),qI=as.vector(qI))
#Not related to C itself, interesting
#Maybe proportion of zeroes?

C_prop<-matrix(nrow=nrep,ncol=nknot)
for (i in chosen_ones){
  for (j in 1:nknot){
    barpidou<-0
    for (b in 1:NY){
      if (sim_save[[i]]$C[j,b]==0) barpidou<-barpidou+1
    }
    C_prop[i,j]<-barpidou/NY
  }
}

prop_C_qI<-data.frame(C=as.vector(C_prop),qI=as.vector(qI))
prop_C_qIdiff<-data.frame(C=as.vector(C_prop),qI_diff=as.vector(qI_diffs))
#So it is not the catch itself that is causing problems

#Maybe it is amount of data?

I_quant<-matrix(nrow=nrep,ncol=nknot)
for (i in chosen_ones){
  for (j in 1:nknot){
    I_quant[i,j]<-length(which(!is.na(sim_save[[i]]$logI[which(sim_save[[i]]$s_i==j-1)])))
  }
}

I_qI<-data.frame(I=as.vector(I_quant),qI=as.vector(qI))
I_qIdiff<-data.frame(I=as.vector(I_quant),qI=as.vector(qI_diffs))
#Nope, maybe avg size of observation?

I_magn<-matrix(nrow=nrep,ncol=nknot)
I_max<-matrix(nrow=nrep,ncol=nknot)
for (i in chosen_ones){
  for (j in 1:nknot){
    I_magn[i,j]<-mean(sim_save[[i]]$logI[(which(sim_save[[i]]$s_i==j-1))],na.rm=T)
    I_max[i,j]<-max(sim_save[[i]]$logI[(which(sim_save[[i]]$s_i==j-1))],na.rm=T)
  }
}
I_magn_qI<-data.frame(I=as.vector(I_magn),qI=as.vector(qI))
I_max_qI<-data.frame(I=as.vector(I_max),qI=as.vector(qI))

diff_B<-matrix(ncol=nrep,nrow=(NY+1))
diff_R<-matrix(ncol=nrep,nrow=NY)
diff_m<-matrix(ncol=nrep,nrow=(NY+1))
perc_diff_m<-matrix(ncol=nrep,nrow=(NY+1))
for (i in chosen_ones) {
  diff_B[,i]<-(Report_list[[i]]$totB-sim_save[[i]]$totB)/sim_save[[i]]$totB
  diff_R[,i]<-(Report_list[[i]]$totR-sim_save[[i]]$totR)/sim_save[[i]]$totR
  diff_m[,i]<-(Report_list[[i]]$mean_m-sim_save[[i]]$mean_m)/sim_save[[i]]$mean_m
}
Bquant<-matrix(ncol=3,nrow=(NY+1))
Rquant<-matrix(ncol=3,nrow=NY)
mquant<-matrix(ncol=3,nrow=(NY+1))
for(i in 1:(NY)){
  Bquant[i,]<-quantile(diff_B[i,],probs=c(0.25,0.5,0.75),na.rm=T)
  if (i <= NY){Rquant[i,]<-quantile(diff_R[i,],probs=c(0.25,0.5,0.75),na.rm=T)}
  mquant[i,]<-quantile(diff_m[i,],probs=c(0.25,0.5,0.75),na.rm=T)
}
Bquant<-as.data.frame(Bquant)
colnames(Bquant)<-c("t","med","n")
Rquant<-as.data.frame(Rquant)
colnames(Rquant)<-c("t","med","n")
mquant<-as.data.frame(mquant)
colnames(mquant)<-c("t","med","n")
Bmean<-mean(diff_B,na.rm=T)
Rmean<-mean(diff_R,na.rm=T)
mmean<-mean(diff_m,na.rm=T)

tot_B_diff_1<-ggplot(Bquant)+geom_ribbon(aes(x=1:(NY+1),ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:(NY+1),y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  # xlab("Year")+ylab("Predicted - True Biomass")+
  # ylim(c(-0.02,0.12))+
  theme_bw()+
  xlab("Year")+ylab("(Predicted - True Total Biomass)/True Total Biomass")
ggsave(filename=paste0(fig_dir,"B_diff_year_1.png"),plot=tot_B_diff_1,
       width=6,height=4)

tot_R_diff_1<-ggplot(Rquant)+geom_ribbon(aes(x=1:NY,ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:NY,y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  # ylim(c(-0.035,0.09))+
  # xlab("Year")+ylab("Predicted - True Recruitment")
  theme_bw()+
  xlab("Year")+ylab("(Predicted - True Total Recruitment)/True Total Recruitment")
ggsave(filename=paste0(fig_dir,"R_diff_year_1.png"),plot=tot_R_diff_1,
       width=6,height=4)

mean_m_diff_1<-ggplot(mquant)+geom_ribbon(aes(x=1:(NY+1),ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:(NY+1),y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  theme_bw()+
  xlab("Year")+ylab("Predicted - True Mean Mortality")
ggsave(filename=paste0(fig_dir,"m_diff_year_1.png"),plot=mean_m_diff_1,
       width=6,height=4)

knot_diff_B<-array(dim=c(nrep,NY+1,nknot))
knot_diff_R<-array(dim=c(nrep,NY,nknot))
knot_diff_m<-array(dim=c(nrep,NY+1,nknot))
for (i in which(!is.na(SigmaO_B))) {
  for (j in 1:(NY+1)){
    for (k in 1:nknot){
      knot_diff_B[i,j,k]<-(Report_list[[i]]$B[k,j]-sim_save[[i]]$B[k,j])/sim_save[[i]]$B[k,j]
      knot_diff_m[i,j,k]<-(Report_list[[i]]$m[k,j]-sim_save[[i]]$m[k,j])/sim_save[[i]]$m[k,j]
      if (j <= NY) knot_diff_R[i,j,k]<-(Report_list[[i]]$R[k,j]-sim_save[[i]]$R[k,j])/sim_save[[i]]$R[k,j]
    }
  }
}
Bquant_big<-matrix(ncol=3,nrow=(NY+1))
Rquant_big<-matrix(ncol=3,nrow=NY)
mquant_big<-matrix(ncol=3,nrow=(NY+1))
for(i in 1:(NY+1)){
  Bquant_big[i,]<-quantile(knot_diff_B[,i,],probs=c(0.25,0.5,0.75),na.rm=T)
  if (i <= NY){Rquant_big[i,]<-quantile(knot_diff_R[,i,],probs=c(0.25,0.5,0.75),na.rm=T)}
  mquant_big[i,]<-quantile(knot_diff_m[,i,],probs=c(0.25,0.5,0.75),na.rm=T)
}
Bquant_big<-as.data.frame(Bquant_big)
colnames(Bquant_big)<-c("t","med","n")
Rquant_big<-as.data.frame(Rquant_big)
colnames(Rquant_big)<-c("t","med","n")
mquant_big<-as.data.frame(mquant_big)
colnames(mquant_big)<-c("t","med","n")
Bmean_big<-mean(knot_diff_B,na.rm=T)
Rmean_big<-mean(knot_diff_R,na.rm=T)
mmean_big<-mean(knot_diff_m,na.rm=T)

#Temporally over all knots
knot_B_diff_year_1<-ggplot(Bquant_big)+geom_ribbon(aes(x=1:(NY+1),ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:(NY+1),y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  # xlab("Year")+ylab("Predicted - True Biomass")+
  # ylim(c(-0.02,0.12))+
  theme_bw()+
  xlab("Year")+ylab("(Predicted - True Biomass Density)/True Biomass Density")
ggsave(filename=paste0(fig_dir,"B_knot_diff_year_1.png"),plot=knot_B_diff_year_1,
       width=6,height=4)

knot_R_diff_year_1<-ggplot(Rquant_big)+geom_ribbon(aes(x=1:NY,ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:NY,y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  # ylim(c(-0.035,0.09))+
  # xlab("Year")+ylab("Predicted - True Recruitment")
  theme_bw()+
  xlab("Year")+ylab("(Predicted - True Recruitment)/True Recruitment")
ggsave(filename=paste0(fig_dir,"R_knot_diff_year_1.png"),plot=knot_R_diff_year_1,
       width=6,height=4)

knot_m_diff_year_1<-ggplot(mquant_big)+geom_ribbon(aes(x=1:(NY+1),ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:(NY+1),y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  theme_bw()+
  xlab("Year")+ylab("Predicted - True Mean Mortality")
ggsave(filename=paste0(fig_dir,"m_knot_diff_year_1.png"),plot=knot_m_diff_year_1,
       width=6,height=4)

#Look for confounding
knot_R<-array(dim=c(nrep,NY,nknot))
knot_m<-array(dim=c(nrep,NY,nknot))
for (i in which(!is.na(SigmaO_B))) {
  for (j in 1:(NY)){
    for (k in 1:nknot){
      knot_m[i,j,k]<-Report_list[[i]]$m[k,j]
      knot_R[i,j,k]<-Report_list[[i]]$R[k,j]
    }
  }
}
vec_R<-as.vector(knot_R)
vec_m<-as.vector(knot_m)
sub_vec_R<-vec_R[-which(is.na(vec_R))]
sub_vec_m<-vec_m[-which(is.na(vec_m))]

#No confounding apparent here at the knot level, maybe at the overall level?
scatter_m_r_1<-ggplot()+
  geom_point(aes(x=log(sub_vec_m),y=log(sub_vec_R)))+
  theme_bw()+xlab("Log Instantaneous Natural Mortality")+
  ylab("Log Recruit Density (kg per km^2)")
ggsave(filename=paste0(fig_dir,"scatter_m_r_1.png"),plot=scatter_m_r_1,
       width=6,height=4)

check_R<-matrix(nrow=nrep,ncol=NY)
check_m<-matrix(nrow=nrep,ncol=NY)
for (i in which(!is.na(SigmaO_B))) {
  check_m[i,]<-Report_list[[i]]$mean_m[-(NY+1)]
  check_R[i,]<-Report_list[[i]]$totR
}
check_vec_R<-as.vector(check_R)
check_vec_m<-as.vector(check_m)
sub_check_R<-check_R[-which(is.na(check_R))]
sub_check_m<-check_m[-which(is.na(check_m))]

#Setting 2:

params$log_qI<-rep(NA,data$n_h)
set.seed(122)
for (i in 1:length(params$log_qI)){
  if (data$h_s[i]==0) params$log_qI[i]<-log(rbeta(1,data$prior_pars[1],data$prior_pars[2]))
  if (data$h_s[i]==1) params$log_qI[i]<-log(rbeta(1,data$prior_pars[3],data$prior_pars[4]))
  if (data$h_s[i]==2) params$log_qI[i]<-log(rbeta(1,data$prior_pars[5],data$prior_pars[6]))
}

random<-c("omega_B","omega_R","omega_m")

maps<-list(log_H_input_R=c(factor(NA),factor(NA)))

obj <- MakeADFun( data=data, parameters=params,random=random, map = maps , DLL="SEBDAM")

fit_params$log_qI<-rep(log(0.3),data$n_s)

for (i in 1:nrep) {
  tryCatch({
    simdata <- obj $ simulate(complete=T)
    sim_save[[i]]<-simdata
    fit_params$log_B0<-log(max(exp(simdata$logI),na.rm=T))+1
    obj2 <- MakeADFun ( simdata , fit_params , random=random, map=maps )
    Opt2 <- try(optimx::optimr(obj2$par,obj2$fn,obj2$gr,
                               control=list(maxit=100000,maxeval=100000),method="nlminb"),T)
    while (Opt2$message=="iteration limit reached without convergence (10)") {
      obj2$par<-obj2$env$last.par[which(names(obj2$env$last.par) %in% non_r)]
      Opt2 <- try(optimx::optimr(obj2$par,obj2$fn,obj2$gr,
                                 control=list(maxit=100000,maxeval=100000),method="nlminb"),T)
    }
    rep_list[[i]] <- sdreport(obj2,bias.correct=F)
    Report_list[[i]] = obj2$report()
    if (!is.null(rep_list[[i]])){
      msg_list[[i]] <- Opt2$message
    }
    if ((i %% 5)==0){save(rep_list,sim_save,msg_list,Report_list,file=paste0("interrupt_sehbam_",i,"_good.RData"))}
  }, error=function(e) {})
}

converge<-0
false<-0
x_conv<-0
singular<-0
null_vec<-nrep
no_conv<-0
chosen_ones<-c()
for (i in 1:nrep){
  if(is.null(msg_list[[i]])) no_conv<-no_conv+1
  if (!is.null(msg_list[[i]])){null_vec<-null_vec-1
  if (msg_list[[i]]=="relative convergence (4)") {converge<-converge+1
  chosen_ones<-c(chosen_ones,i)}
  if (msg_list[[i]]=="false convergence (8)") {false<-false+1}
  if(msg_list[[i]]=="singular convergence (7)"){singular<-singular+1}
  if(msg_list[[i]]=="both X-convergence and relative convergence (5)"){x_conv<-x_conv+1
  chosen_ones<-c(chosen_ones,i)}
  if(msg_list[[i]]=="nlminb failed"){no_conv=no_conv+1}
  }
}
conv_frame_2<-data.frame(converge=converge,false=false,singular=singular,x=x_conv,no_conv=no_conv)

sigma_epsilon<-c(rep(NA,nrep))
sigma_upsilon<-c(rep(NA,nrep))
S<-matrix(ncol=n_hab,nrow=nrep)
qI<-matrix(rep(NA,nrep*n_hab),ncol=n_hab)
qR<-matrix(ncol=n_hab,nrow=nrep)
mu_I<-c(rep(NA,nrep))
mu_IR<-matrix(ncol=n_hab,nrow=nrep)
log_var_I<-c(rep(NA,nrep))
var_IR<-matrix(ncol=n_hab,nrow=nrep)
B0<-c(rep(NA,nrep))
R0<-c(rep(NA,nrep))
m0<-c(rep(NA,nrep))
kappa_B<-c(rep(NA,nrep))
kappa_R<-c(rep(NA,nrep))
kappa_m<-c(rep(NA,nrep))
tau_B<-c(rep(NA,nrep))
tau_m<-c(rep(NA,nrep))
tau_R<-c(rep(NA,nrep))
log_H_input_B1<-c(rep(NA,nrep))
log_H_input_B2<-c(rep(NA,nrep))
log_H_input_m1<-c(rep(NA,nrep))
log_H_input_m2<-c(rep(NA,nrep))
SigmaO_B<-rep(NA,nrep)
SigmaO_R<-rep(NA,nrep)
SigmaO_m<-rep(NA,nrep)
#to start the parameters
for (i in chosen_ones){
  sigma_epsilon[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="sigma_epsilon")]
  sigma_upsilon[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="sigma_upsilon")]
  S[i,]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="S")]
  qR[i,]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="qR")]
  mu_I[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="mu_I")]
  mu_IR[i,]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="mu_IR")]
  log_var_I[i]<-log(rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="var_I")])
  var_IR[i,]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="var_IR")]
  B0[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="B0")]
  R0[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="R0")]
  m0[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="m0")]
  kappa_B[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="kappa_B")]
  kappa_R[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="kappa_R")]
  kappa_m[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="kappa_m")]
  tau_B[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="tau_B")]
  tau_R[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="tau_R")]
  tau_m[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="tau_m")]
  log_H_input_B1[i]<-Report_list[[i]]$log_H_input_B[1]
  log_H_input_B2[i]<-Report_list[[i]]$log_H_input_B[2]
  log_H_input_m1[i]<-Report_list[[i]]$log_H_input_m[1]
  log_H_input_m2[i]<-Report_list[[i]]$log_H_input_m[2]
  SigmaO_B[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="SigmaO_B")]
  SigmaO_R[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="SigmaO_R")]
  SigmaO_m[i]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="SigmaO_m")]
  # }
}

for (i in chosen_ones){
  for (j in 1:n_hab){
    qI[i,j]<-rep_list[[i]]$value[which(names(rep_list[[i]]$value)=="qI")][j]
  }
}

qI_diffs<-matrix(rep(NA,nrep*n_hab),ncol=n_hab)
qR_diffs<-matrix(ncol=n_hab,nrow=nrep)
mu_IR_diffs<-matrix(ncol=n_hab,nrow=nrep)
log_var_IR_diffs<-matrix(ncol=n_hab,nrow=nrep)
S_diffs<-matrix(ncol=n_hab,nrow=nrep)
for (i in chosen_ones){
  qI_diffs[i,]<-qI[i,]-exp(sim_save[[i]]$log_qI)
  qR_diffs[i,]<-qR[i,]-sim_save[[i]]$qR
  mu_IR_diffs[i,]<-mu_IR[i,]-params$mu_IR
  log_var_IR_diffs[i,]<-log(var_IR[i,])-params$log_var_IR
  S_diffs[i,]<-S[i,]-sim_save[[i]]$S
}

all_params<-data.frame(sigma_epsilon,sigma_upsilon,
                       mu_I,log_var_I,
                       B0,R0,m0,kappa_B,
                       kappa_R,kappa_m,tau_B,tau_m,tau_R,
                       log_H_input_B1,log_H_input_B2,
                       log_H_input_m1,log_H_input_m2,
                       SigmaO_B,SigmaO_R,SigmaO_m)

hist_plot_2<-data.frame(value=c(sigma_epsilon,sigma_upsilon,
                                mu_I,log_var_I,B0,R0,m0,
                                kappa_B,kappa_R,kappa_m,
                                tau_B,tau_m,tau_R,
                                log_H_input_B1,log_H_input_B2,
                                log_H_input_m1,log_H_input_m2,
                                SigmaO_B,SigmaO_R,SigmaO_m),
                        parameter=rep(c("sigma_epsilon","sigma_upsilon",
                                        "mu_I","log_var_I","B0",
                                        "R0","m0",
                                        "kappa_B",
                                        "kappa_R",
                                        "kappa_m","tau_B",
                                        "tau_m","tauR",
                                        "H_input_B1","H_input_B2",
                                        "H_input_m1","H_input_m2",
                                        "SigmaO_B","SigmaO_R","SigmaO_m")
                                      ,each=nrep),
                        true=rep(c(exp(params$log_sigma_epsilon),
                                   exp(params$log_sigma_upsilon),params$mu_I,
                                   params$log_var_I,
                                   exp(params$log_B0),
                                   exp(params$log_R0),exp(params$log_m0),
                                   exp(params$log_kappa_B),exp(params$log_kappa_R),
                                   exp(params$log_kappa_m),
                                   exp(params$log_tau_B),exp(params$log_tau_m),exp(params$log_tau_R),
                                   params$log_H_input_B[1],params$log_H_input_B[2],
                                   params$log_H_input_m[1],params$log_H_input_m[2],
                                   1/((exp(params$log_tau_B))*sqrt(4*pi)*exp(params$log_kappa_B)),
                                   1/((exp(params$log_tau_R))*sqrt(4*pi)*exp(params$log_kappa_R)),
                                   1/((exp(params$log_tau_m))*sqrt(4*pi)*exp(params$log_kappa_m))),each=nrep))

hist_plot_2$parameter<-factor(hist_plot_2$parameter,
                              labels=c(expression(B [0]),expression("H" ["input"]^"B1"),
                                       expression("H" ["input"]^"B2"),
                                       expression("H" ["input"]^"m1"),expression("H" ["input"]^"m2"),
                                       expression(kappa ["B"]),expression(kappa ["m"]),
                                       expression(kappa ["R"]),expression(log (sigma["I"]^2)),
                                       expression("m" [0]),
                                       expression(mu["I"]),
                                       expression("R" [0]),
                                       expression(sigma [epsilon]),
                                       expression(sigma[upsilon]),
                                       expression(sigma[B]^Omega),
                                       expression(sigma[m]^Omega),
                                       expression(sigma[R]^Omega),
                                       expression(tau ["B"]),
                                       expression(tau ["m"]),expression(tau["R"])))

par_hist_set_2<-ggplot(hist_plot_2)+geom_histogram(aes(x=value),col="black",fill="grey")+
  geom_vline(aes(xintercept = true,col="True Value"))+
  theme(strip.background = element_blank(),legend.position = "none")+
  facet_wrap(~parameter,scales="free",labeller=label_parsed)+
  theme_bw()+
  xlab("Estimated Value")+ylab("Frequency")
ggsave(filename=paste0(fig_dir,"par_hist_2.png"),plot=par_hist_set_2,
       width=9,height=8)

plot_qI<-data.frame(value=as.vector(qI_diffs),
                    param=rep("qI",n_hab*nrep))
plot_qR<-data.frame(value=as.vector(qR_diffs),
                    param=rep("qR",n_hab*nrep))
plot_mu_IR<-data.frame(value=as.vector(mu_IR_diffs),
                       param=rep("mu_IR",n_hab*nrep))
plot_log_var_IR<-data.frame(value=as.vector(log_var_IR_diffs),
                            param=rep("log_var_IR",n_hab*nrep))
plot_S<-data.frame(value=as.vector(S_diffs),
                   param=rep("S",n_hab*nrep))

diff_hist_plot_2<-rbind(plot_qI,plot_qR,plot_mu_IR,plot_log_var_IR,plot_S)

diff_hist_2<-ggplot(diff_hist_plot_2)+geom_histogram(aes(x=value),col="black",fill="grey")+
  geom_vline(aes(xintercept = 0),lty="dashed")+
  theme(strip.background = element_blank(),legend.position = "none")+
  facet_wrap(~param,scales="free",labeller=label_parsed)+
  theme_bw()+
  xlab("Estimated Value")+ylab("Frequency")
ggsave(filename=paste0(fig_dir,"diff_hist_2.png"),plot=diff_hist_2,
       width=7,height=6)

plot_qI2<-data.frame(value=as.vector(qI_diffs),
                     hab=rep(1:3,each=nrep))

qI_by_hab_2<-ggplot(plot_qI2)+geom_histogram(aes(x=value),col="black",fill="grey")+
  geom_vline(aes(xintercept = 0),lty="dashed")+
  theme(strip.background = element_blank(),legend.position = "none")+
  facet_wrap(~hab,scales="free")+
  theme_bw()+
  xlab("Estimated Value")+ylab("Frequency")
ggsave(filename=paste0(fig_dir,"qI_by_hab_2.png"),plot=qI_by_hab_2,
       width=6,height=4)

diff_2_B<-matrix(ncol=nrep,nrow=(NY+1))
diff_2_R<-matrix(ncol=nrep,nrow=NY)
diff_2_m<-matrix(ncol=nrep,nrow=(NY+1))
perc_diff_2_m<-matrix(ncol=nrep,nrow=(NY+1))
for (i in chosen_ones) {
  diff_2_B[,i]<-(Report_list[[i]]$totB-sim_save[[i]]$totB)/sim_save[[i]]$totB
  diff_2_R[,i]<-(Report_list[[i]]$totR-sim_save[[i]]$totR)/sim_save[[i]]$totR
  diff_2_m[,i]<-(Report_list[[i]]$mean_m-sim_save[[i]]$mean_m)/sim_save[[i]]$mean_m
}
Bquant_2<-matrix(ncol=3,nrow=(NY+1))
Rquant_2<-matrix(ncol=3,nrow=NY)
mquant_2<-matrix(ncol=3,nrow=(NY+1))
for(i in 1:(NY)){
  Bquant_2[i,]<-quantile(diff_2_B[i,],probs=c(0.25,0.5,0.75),na.rm=T)
  if (i <= NY){Rquant_2[i,]<-quantile(diff_2_R[i,],probs=c(0.25,0.5,0.75),na.rm=T)}
  mquant_2[i,]<-quantile(diff_2_m[i,],probs=c(0.25,0.5,0.75),na.rm=T)
}
Bquant_2<-as.data.frame(Bquant_2)
colnames(Bquant_2)<-c("t","med","n")
Rquant_2<-as.data.frame(Rquant_2)
colnames(Rquant_2)<-c("t","med","n")
mquant_2<-as.data.frame(mquant_2)
colnames(mquant_2)<-c("t","med","n")
Bmean2<-mean(diff_2_B,na.rm=T)
Rmean2<-mean(diff_2_R,na.rm=T)
mmean2<-mean(diff_2_m,na.rm=T)

tot_B_diff_2<-ggplot(Bquant_2)+geom_ribbon(aes(x=1:(NY+1),ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:(NY+1),y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  # xlab("Year")+ylab("Predicted - True Biomass")+
  # ylim(c(-0.02,0.12))+
  theme_bw()+
  xlab("Year")+ylab("(Predicted - True Total Biomass)/True Total Biomass")
ggsave(filename=paste0(fig_dir,"B_diff_year_2.png"),plot=tot_B_diff_2,
       width=6,height=4)

tot_R_diff_2<-ggplot(Rquant_2)+geom_ribbon(aes(x=1:NY,ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:NY,y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  theme_bw()+
  # ylim(c(-0.035,0.09))+
  # xlab("Year")+ylab("Predicted - True Recruitment")
  xlab("Year")+ylab("(Predicted - True Total Recruitment)/True Total Recruitment")
ggsave(filename=paste0(fig_dir,"R_diff_year_2.png"),plot=tot_R_diff_2,
       width=6,height=4)

mean_m_diff_2<-ggplot(mquant_2)+geom_ribbon(aes(x=1:(NY+1),ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:(NY+1),y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  theme_bw()+
  xlab("Year")+ylab("Predicted - True Mean Mortality")
ggsave(filename=paste0(fig_dir,"m_diff_year_2.png"),plot=mean_m_diff_2,
       width=6,height=4)

knot_diff_2_B<-array(dim=c(nrep,NY+1,nknot))
knot_diff_2_R<-array(dim=c(nrep,NY,nknot))
knot_diff_2_m<-array(dim=c(nrep,NY+1,nknot))
for (i in which(!is.na(SigmaO_B))) {
  for (j in 1:(NY+1)){
    for (k in 1:nknot){
      knot_diff_2_B[i,j,k]<-(Report_list[[i]]$B[k,j]-sim_save[[i]]$B[k,j])/sim_save[[i]]$B[k,j]
      knot_diff_2_m[i,j,k]<-(Report_list[[i]]$m[k,j]-sim_save[[i]]$m[k,j])/sim_save[[i]]$m[k,j]
      if (j <= NY) knot_diff_2_R[i,j,k]<-(Report_list[[i]]$R[k,j]-sim_save[[i]]$R[k,j])/sim_save[[i]]$R[k,j]
    }
  }
}
Bquant_2_big<-matrix(ncol=3,nrow=(NY+1))
Rquant_2_big<-matrix(ncol=3,nrow=NY)
mquant_2_big<-matrix(ncol=3,nrow=(NY+1))
for(i in 1:(NY+1)){
  Bquant_2_big[i,]<-quantile(knot_diff_2_B[,i,],probs=c(0.25,0.5,0.75),na.rm=T)
  if (i <= NY){Rquant_2_big[i,]<-quantile(knot_diff_2_R[,i,],probs=c(0.25,0.5,0.75),na.rm=T)}
  mquant_2_big[i,]<-quantile(knot_diff_2_m[,i,],probs=c(0.25,0.5,0.75),na.rm=T)
}
Bquant_2_big<-as.data.frame(Bquant_2_big)
colnames(Bquant_2_big)<-c("t","med","n")
Rquant_2_big<-as.data.frame(Rquant_2_big)
colnames(Rquant_2_big)<-c("t","med","n")
mquant_2_big<-as.data.frame(mquant_2_big)
colnames(mquant_2_big)<-c("t","med","n")
Bmean_big2<-mean(knot_diff_2_B,na.rm=T)
Rmean_big2<-mean(knot_diff_2_R,na.rm=T)
mmean_big2<-mean(knot_diff_2_m,na.rm=T)

#Look for confounding
knot_2_R<-array(dim=c(nrep,NY,nknot))
knot_2_m<-array(dim=c(nrep,NY,nknot))
for (i in which(!is.na(SigmaO_B))) {
  for (j in 1:(NY)){
    for (k in 1:nknot){
      knot_2_m[i,j,k]<-Report_list[[i]]$m[k,j]
      knot_2_R[i,j,k]<-Report_list[[i]]$R[k,j]
    }
  }
}
vec_2_R<-as.vector(knot_2_R)
vec_2_m<-as.vector(knot_2_m)
sub_vec_2_R<-vec_2_R[-which(is.na(vec_2_R))]
sub_vec_2_m<-vec_2_m[-which(is.na(vec_2_m))]

#No confounding apparent here at the knot level, maybe at the overall level?
scatter_m_r_2<-ggplot()+
  geom_point(aes(x=log(sub_vec_2_m),y=log(sub_vec_2_R)))+
  theme_bw()+xlab("Log Instantaneous Natural Mortality")+
  ylab("Log Recruit Density (kg per km^2)")
ggsave(filename=paste0(fig_dir,"scatter_m_r_2.png"),plot=scatter_m_r_2,
       width=6,height=4)

#Temporally over all knots
knot_B_diff_year_2<-ggplot(Bquant_2_big)+geom_ribbon(aes(x=1:(NY+1),ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:(NY+1),y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  theme_bw()+
  # xlab("Year")+ylab("Predicted - True Biomass")+
  # ylim(c(-0.02,0.12))+
  xlab("Year")+ylab("(Predicted - Simulated Biomass Density)/Simulated Biomass Density")
ggsave(filename=paste0(fig_dir,"B_knot_diff_year_2.png"),plot=knot_B_diff_year_2,
       width=6,height=4)

knot_R_diff_year_2<-ggplot(Rquant_2_big)+geom_ribbon(aes(x=1:NY,ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:NY,y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  theme_bw()+
  # ylim(c(-0.035,0.09))+
  # xlab("Year")+ylab("Predicted - Simulated Recruitment")
  xlab("Year")+ylab("(Predicted - Simulated Recruitment)/Simulated Recruitment")
ggsave(filename=paste0(fig_dir,"R_knot_diff_year_2.png"),plot=knot_R_diff_year_2,
       width=6,height=4)

knot_m_diff_year_2<-ggplot(mquant_2_big)+geom_ribbon(aes(x=1:(NY+1),ymin=t,ymax=n),alpha=0.5,col="grey")+
  geom_line(aes(x=1:(NY+1),y=med))+
  geom_hline(aes(yintercept=0), lty="dashed")+
  theme_bw()+
  xlab("Year")+ylab("Predicted - Simulated Mean Mortality")
ggsave(filename=paste0(fig_dir,"m_knot_diff_year_2.png"),plot=knot_m_diff_year_2,
       width=6,height=4)

#Combined:
#Overall
B_box<-ggplot()+stat_boxplot(geom="errorbar",aes(y=as.vector(diff_B),x="Setting 1"))+
  geom_boxplot(aes(y=as.vector(diff_B),x="Setting 1"),outlier.shape = NA)+
  stat_boxplot(geom="errorbar",aes(y=as.vector(diff_2_B),x="Setting 2"))+
  geom_boxplot(aes(y=as.vector(diff_2_B),x="Setting 2"),outlier.shape = NA)+
  coord_cartesian(y=c(-0.5,0.5))+
  geom_point(aes(x="Setting 1",y=Bmean))+
  geom_point(aes(x="Setting 2",y=Bmean2))+
  theme_bw()+
  geom_hline(aes(yintercept=0),lty="dashed")+
  xlab("Experiment")+ylab("(Predicted - Simulated Biomass Density)/Simulated Biomass Density")
ggsave(filename=paste0(fig_dir,"boxplot_B_both.png"),plot=B_box,
       width=8,height=6)

R_box<-ggplot()+stat_boxplot(geom="errorbar",aes(y=as.vector(diff_R),x="Setting 1"))+
  geom_boxplot(aes(y=as.vector(diff_R),x="Setting 1"),outlier.shape = NA)+
  stat_boxplot(geom="errorbar",aes(y=as.vector(diff_2_R),x="Setting 2"))+
  geom_boxplot(aes(y=as.vector(diff_2_R),x="Setting 2"),outlier.shape = NA)+
  coord_cartesian(y=c(-0.5,1))+
  geom_point(aes(x="Setting 1",y=Rmean))+
  geom_point(aes(x="Setting 2",y=Rmean2))+
  theme_bw()+
  geom_hline(aes(yintercept=0),lty="dashed")+
  xlab("Experiment")+ylab("(Predicted - Simulated Total Recruitment)/Simulated Total Recruitment")
ggsave(filename=paste0(fig_dir,"boxplot_R_both.png"),plot=R_box,
       width=8,height=6)

m_box<-ggplot()+stat_boxplot(geom="errorbar",aes(y=as.vector(diff_m),x="Setting 1"))+
  geom_boxplot(aes(y=as.vector(diff_m),x="Setting 1"),outlier.shape = NA)+
  stat_boxplot(geom="errorbar",aes(y=as.vector(diff_2_m),x="Setting 2"))+
  geom_boxplot(aes(y=as.vector(diff_2_m),x="Setting 2"),outlier.shape = NA)+
  coord_cartesian(y=c(-0.75,1.25))+
  geom_point(aes(x="Setting 1",y=mmean))+
  geom_point(aes(x="Setting 2",y=mmean2))+
  theme_bw()+
  geom_hline(aes(yintercept=0),lty="dashed")+
  xlab("Experiment")+ylab("Predicted - Simulated Mean Natural Mortality")
ggsave(filename=paste0(fig_dir,"boxplot_m_both.png"),plot=m_box,
       width=8,height=6)

#Knot based:
B_box_knot<-ggplot()+stat_boxplot(geom="errorbar",aes(y=as.vector(knot_diff_B),x="Setting 1"))+
  geom_boxplot(aes(y=as.vector(knot_diff_B),x="Setting 1"),outlier.shape = NA)+
  stat_boxplot(geom="errorbar",aes(y=as.vector(knot_diff_2_B),x="Setting 2"))+
  geom_boxplot(aes(y=as.vector(knot_diff_2_B),x="Setting 2"),outlier.shape = NA)+
  coord_cartesian(y=c(-0.9,1))+
  geom_point(aes(x="Setting 1",y=Bmean_big))+
  geom_point(aes(x="Setting 2",y=Bmean_big2))+
  theme_bw()+
  geom_hline(aes(yintercept=0),lty="dashed")+
  xlab("Experiment")+ylab("(Predicted - Simulated Biomass Density)/Simulated Biomass Density")
ggsave(filename=paste0(fig_dir,"boxplot_B_knot_both.png"),plot=B_box_knot,
       width=8,height=6)

R_box_knot<-ggplot()+stat_boxplot(geom="errorbar",aes(y=as.vector(knot_diff_R),x="Setting 1"))+
  geom_boxplot(aes(y=as.vector(knot_diff_R),x="Setting 1"),outlier.shape = NA)+
  stat_boxplot(geom="errorbar",aes(y=as.vector(knot_diff_2_R),x="Setting 2"))+
  geom_boxplot(aes(y=as.vector(knot_diff_2_R),x="Setting 2"),outlier.shape = NA)+
  coord_cartesian(y=c(-1,2.5))+
  geom_point(aes(x="Setting 1",y=Rmean_big))+
  geom_point(aes(x="Setting 2",y=Rmean_big2))+
  theme_bw()+
  geom_hline(aes(yintercept=0),lty="dashed")+
  xlab("Experiment")+ylab("(Predicted - Simulated Total Recruitment)/Simulated Total Recruitment")
ggsave(filename=paste0(fig_dir,"boxplot_R_knot_both.png"),plot=R_box_knot,
       width=8,height=6)

m_box_knot<-ggplot()+stat_boxplot(geom="errorbar",aes(y=as.vector(knot_diff_m),x="Setting 1"))+
  geom_boxplot(aes(y=as.vector(knot_diff_m),x="Setting 1"),outlier.shape = NA)+
  stat_boxplot(geom="errorbar",aes(y=as.vector(knot_diff_2_m),x="Setting 2"))+
  geom_boxplot(aes(y=as.vector(knot_diff_2_m),x="Setting 2"),outlier.shape = NA)+
  coord_cartesian(y=c(-1,1.5))+
  geom_point(aes(x="Setting 1",y=mmean_big))+
  geom_point(aes(x="Setting 2",y=mmean_big2))+
  theme_bw()+
  geom_hline(aes(yintercept=0),lty="dashed")+
  xlab("Experiment")+ylab("Predicted - Simulated Mean Natural Mortality")
ggsave(filename=paste0(fig_dir,"boxplot_m_knot_both.png"),plot=m_box_knot,
       width=8,height=6)

#Combining both parameter estimates plots
hist_plot$setting<-rep("1",length(hist_plot$value))
hist_plot_2$setting<-rep("2",length(hist_plot_2$value))
hist_plot_both<-rbind(hist_plot,hist_plot_2)

par_hist_both<-ggplot(hist_plot_both)+geom_histogram(aes(x=value,col=setting,fill=setting),alpha=0.2)+
  geom_vline(aes(xintercept = true),col="black",lty="dashed")+
  theme(strip.background = element_blank(),legend.position = "none")+
  facet_wrap(~parameter,scales="free",labeller=label_parsed)+
  theme_bw()+
  scale_color_discrete(name="Setting",labels=c("1","2"))+
  scale_fill_discrete(name="Setting",labels=c("1","2"))+
  xlab("Value")+ylab("Frequency")
ggsave(filename=paste0(fig_dir,"par_hist_both.png"),plot=par_hist_both,
       width=9,height=8)

diff_hist_plot$setting<-rep("1",length(diff_hist_plot$value))
diff_hist_plot$param<-factor(diff_hist_plot$param,
                             labels=c(expression(log (sigma ["IR"])),
                                      expression(mu ["IR"]),
                                      expression("q"^"I"),
                                      expression("q"^"R"),
                                      expression("S")))
diff_hist_plot_2$setting<-rep("2",length(diff_hist_plot_2$value))
diff_hist_plot_2$param<-factor(diff_hist_plot_2$param,
                               labels=c(expression(log (sigma ["IR"])),
                                        expression(mu ["IR"]),
                                        expression("q"^"I"),
                                        expression("q"^"R"),
                                        expression("S")))
diff_hist_plot_both<-rbind(diff_hist_plot,diff_hist_plot_2)

diff_hist_both<-ggplot(diff_hist_plot_both)+geom_histogram(aes(x=value,col=setting,fill=setting),alpha=0.2)+
  geom_vline(aes(xintercept = 0),lty="dashed")+
  theme(strip.background = element_blank(),legend.position = "none")+
  facet_wrap(~param,scales="free",labeller=label_parsed)+
  theme_bw()+
  scale_color_discrete(name="Setting",labels=c("1","2"))+
  scale_fill_discrete(name="Setting",labels=c("1","2"))+
  xlab("Value")+ylab("Frequency")
ggsave(filename=paste0(fig_dir,"diff_hist_both.png"),plot=diff_hist_both,
       width=7,height=6)
