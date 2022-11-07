#Plots from combining SEHBAM, TLM, and SEBDAM outputs
library(ggplot2)
library(sf)

nknot<-60

#Comparative qIs

# load("SEBDAM_BF_fit_60_knots.RData")
load("SEBDAM_BF_fit_60_knots_habitat_prior_converged.RData")

sebdam_rep<-rep
sebdam_Report<-Report

sebdam_qIs<-data.frame(qI=rep$value[which(names(rep$value)=="qI")],
                       sd=rep$sd[which(names(rep$value)=="qI")],
                       source=rep("SEAM",60),
                       knotID=unique(obj$env$data$s_i))

load("manuscript_SEHBAM_fit.RData")

sehbam_rep<-rep
sehbam_Report<-Report

sehbam_qIs<-data.frame(qI=rep$value[which(names(rep$value)=="qI")],
                       sd=rep$sd[which(names(rep$value)=="qI")],
                       source=rep("SEHBAM",60),
                       knotID=unique(obj$env$data$s_i))

all_qIs<-rbind(sebdam_qIs,sehbam_qIs)

#Means
mean(subset(all_qIs,source=="SEAM")$qI)
mean(subset(all_qIs,source=="SEHBAM")$qI)
length(which(sehbam_qIs$qI<sebdam_qIs$qI))/60

fig_dir<-"C:/Users/mcdonaldra/Documents/Github/SEHBAM-paper/Figures/"

jitter_guy<-position_jitter(width=0.5,seed=23)

qI_comp_plot<-ggplot(data=all_qIs,aes(x=knotID,y=qI,col=source,shape=source,
                                      ymin=qI-sd,ymax=qI+sd))+
  geom_point()+
  geom_errorbar(width=0.65)+
  geom_hline(aes(yintercept=0.3078),lty="twodash")+
  geom_hline(aes(yintercept=mean(sebdam_qIs$qI)),col="red",lty="dashed")+
  geom_hline(aes(yintercept=mean(sehbam_qIs$qI)),col="blue",lty="longdash")+
  theme_bw()+
  scale_color_manual(name="Model",values=c("red","blue"))+
  scale_shape_discrete(name="Model")+
  xlab("Knot")+ylab("Value")
# ggsave(filename=paste0(fig_dir,"qI_comp_plot.png"),
# plot=qI_comp_plot,height=8,width=10)
ggsave(filename=paste0(fig_dir,"qI_comp_plot_supp.png"),
       plot=qI_comp_plot,height=8,width=10)

# load("TLM_BF_fit_manuscript.RData")
load("TLM_BF_fit_manuscript_diff_prior.RData")

tlm_rep<-rep
tlm_Report<-Report

NY<-length(tlm_Report$B)-1

#Total bio

bio_tots<-data.frame(tlm=tlm_Report$B,
                     sebdam=sebdam_Report$totB,
                     sehbam=sehbam_Report$totB,
                     tlm_se=tlm_rep$sd[which(names(tlm_rep$value)=="log_B")],
                     sebdam_se=sebdam_rep$sd[which(names(sebdam_rep$value)=="log_totB")],
                     sehbam_se=sehbam_rep$sd[which(names(sehbam_rep$value)=="log_totB")])

rec_tots<-data.frame(tlm=tlm_Report$R[-20],
                     sebdam=sebdam_Report$totR,
                     sehbam=sehbam_Report$totR,
                     tlm_se=tlm_rep$sd[which(names(tlm_rep$value)=="log_R")][-20],
                     sebdam_se=sebdam_rep$sd[which(names(sebdam_rep$value)=="log_totR")],
                     sehbam_se=sehbam_rep$sd[which(names(sehbam_rep$value)=="log_totR")])

m_tots<-data.frame(tlm=tlm_Report$m,
                   sebdam=sebdam_Report$mean_m,
                   sehbam=sehbam_Report$mean_m,
                   tlm_se=tlm_rep$sd[which(names(tlm_rep$value)=="log_m")],
                   sebdam_se=sebdam_rep$sd[which(names(sebdam_rep$value)=="log_mean_m")],
                   sehbam_se=sehbam_rep$sd[which(names(sehbam_rep$value)=="log_mean_m")])

tot_bio_plot<-ggplot(data=bio_tots[-20,])+
  geom_point(aes(x=2001:2019,y=tlm,col="TLM",shape="TLM"))+
  geom_line(aes(x=2001:2019,y=tlm,col="TLM"))+
  geom_ribbon(aes(x=2001:2019,ymax=exp(log(tlm)+2*tlm_se),
                  ymin=exp(log(tlm)-2*tlm_se),col="TLM",fill="TLM"),alpha=0.2)+
  geom_point(aes(x=2001:2019,y=sebdam,col="SEAM",shape="SEAM"))+
  geom_line(aes(x=2001:2019,y=sebdam,col="SEAM"))+
  geom_ribbon(aes(x=2001:2019,ymax=exp(log(sebdam)+2*sebdam_se),
                  ymin=exp(log(sebdam)-2*sebdam_se),col="SEAM",fill="SEAM"),alpha=0.2)+
  geom_point(aes(x=2001:2019,y=sehbam,col="SEHBAM",shape="SEHBAM"))+
  geom_line(aes(x=2001:2019,y=sehbam,col="SEHBAM"))+
  geom_ribbon(aes(x=2001:2019,ymax=exp(log(sehbam)+2*sehbam_se),
                  ymin=exp(log(sehbam)-2*sehbam_se),col="SEHBAM",fill="SEHBAM"),alpha=0.2)+
  scale_fill_viridis_d(name="Model")+
  scale_shape_discrete(name="Model")+
  scale_color_viridis_d(name="Model")+
  theme_bw()+xlab("Year")+ylab("Total Biomass (metric tonnes)")
ggsave(filename=paste0(fig_dir,"tot_B_plot.png"),
       plot=tot_bio_plot,height=7,width=9)
# ggsave(filename=paste0(fig_dir,"tot_B_plot_supp.png"),
#        plot=tot_bio_plot,height=7,width=9)

#Diff in uncertainties
mean(bio_tots$sehbam_se/bio_tots$sehbam)
mean(bio_tots$sebdam_se/bio_tots$sebdam)
mean(bio_tots$tlm_se/bio_tots$tlm)

tot_rec_plot<-ggplot(data=rec_tots)+
  geom_point(aes(x=2001:2019,y=tlm,col="TLM",shape="TLM"))+
  geom_line(aes(x=2001:2019,y=tlm,col="TLM"))+
  geom_ribbon(aes(x=2001:2019,ymax=exp(log(tlm)+2*tlm_se),
                  ymin=exp(log(tlm)-2*tlm_se),col="TLM",fill="TLM"),alpha=0.2)+
  geom_point(aes(x=2001:2019,y=sebdam,col="SEAM",shape="SEAM"))+
  geom_line(aes(x=2001:2019,y=sebdam,col="SEAM"))+
  geom_ribbon(aes(x=2001:2019,ymax=exp(log(sebdam)+2*sebdam_se),
                  ymin=exp(log(sebdam)-2*sebdam_se),col="SEAM",fill="SEAM"),alpha=0.2)+
  geom_point(aes(x=2001:2019,y=sehbam,col="SEHBAM",shape="SEHBAM"))+
  geom_line(aes(x=2001:2019,y=sehbam,col="SEHBAM"))+
  geom_ribbon(aes(x=2001:2019,ymax=exp(log(sehbam)+2*sehbam_se),
                  ymin=exp(log(sehbam)-2*sehbam_se),col="SEHBAM",fill="SEHBAM"),alpha=0.2)+
  scale_fill_viridis_d(name="Model")+
  scale_shape_discrete(name="Model")+
  scale_color_viridis_d(name="Model")+
  theme_bw()+xlab("Year")+ylab("Total Recruitment (metric tonnes)")
# ggsave(filename=paste0(fig_dir,"tot_R_plot.png"),
#        plot=tot_rec_plot,height=7,width=9)
# ggsave(filename=paste0(fig_dir,"tot_R_plot_supp.png"),
#        plot=tot_rec_plot,height=7,width=9)

tot_m_plot<-ggplot(data=m_tots[-20,])+
  geom_point(aes(x=2001:2019,y=tlm,col="TLM",shape="TLM"))+
  geom_line(aes(x=2001:2019,y=tlm,col="TLM"))+
  geom_ribbon(aes(x=2001:2019,ymax=exp(log(tlm)+2*tlm_se),
                  ymin=exp(log(tlm)-2*tlm_se),col="TLM",fill="TLM"),alpha=0.2)+
  geom_point(aes(x=2001:2019,y=sebdam,col="SEAM",shape="SEAM"))+
  geom_line(aes(x=2001:2019,y=sebdam,col="SEAM"))+
  geom_ribbon(aes(x=2001:2019,ymax=exp(log(sebdam)+2*sebdam_se),
                  ymin=exp(log(sebdam)-2*sebdam_se),col="SEAM",fill="SEAM"),alpha=0.2)+
  geom_point(aes(x=2001:2019,y=sehbam,col="SEHBAM",shape="SEHBAM"))+
  geom_line(aes(x=2001:2019,y=sehbam,col="SEHBAM"))+
  geom_ribbon(aes(x=2001:2019,ymax=exp(log(sehbam)+2*sehbam_se),
                  ymin=exp(log(sehbam)-2*sehbam_se),col="SEHBAM",fill="SEHBAM"),alpha=0.2)+
  scale_fill_viridis_d(name="Model")+
  scale_shape_discrete(name="Model")+
  scale_color_viridis_d(name="Model")+
  theme_bw()+xlab("Year")+ylab("Mean Instantaneous Natural Mortality")
# ggsave(filename=paste0(fig_dir,"tot_m_plot.png"),
#        plot=tot_m_plot,height=7,width=9)
# ggsave(filename=paste0(fig_dir,"tot_m_plot_supp.png"),
#        plot=tot_m_plot,height=7,width=9)

#Spatial plots
load("spat_plot_grids.RData")

temp_bio<-data.frame(knotID=rep(1:60,NY+1),
                     Year=rep(2001:2020,each=nknot),
                     B=as.vector(sehbam_Report$B))
Bio_spat<-left_join(temp_bio,gridded_bound_shape,by="knotID") %>% st_as_sf()

temp_rec<-data.frame(knotID=rep(1:60,NY),
                     Year=rep(2001:2019,each=nknot),
                     R=as.vector(sehbam_Report$R))
Rec_spat<-left_join(temp_rec,gridded_bound_shape,by="knotID") %>% st_as_sf()

temp_m<-data.frame(knotID=rep(1:60,NY+1),
                   Year=rep(2001:2020,each=nknot),
                   m=as.vector(sehbam_Report$m))
m_spat<-left_join(temp_m,gridded_bound_shape,by="knotID") %>% st_as_sf()

plot_map<-st_transform(sf_good_map,crs=32619)
plot_map$geometry<-plot_map$geometry/1000
st_crs(plot_map)<-NA

spatial_B_plot<-ggplot()+
  geom_sf(data=Bio_spat,aes(fill=B),col=NA)+
  facet_wrap(~Year)+
  scale_fill_viridis_c(name="Predicted Biomass \nDensity (kg/km^2)")+
  geom_sf(data=plot_map)+theme_bw()+
  coord_sf(xlim=c(690,890),ylim=c(4870,5080))+
  xlab("Easting")+ylab("Northing")
# ggsave(filename=paste0(fig_dir,"spat_B_plot.png"),
#        plot=spatial_B_plot,height=12,width=12)

log_spatial_B_plot<-ggplot()+
  geom_sf(data=Bio_spat,aes(fill=log(B)),col=NA)+
  facet_wrap(~Year)+
  scale_fill_viridis_c(name="Log Predicted Biomass \nDensity (kg/km^2)")+
  geom_sf(data=plot_map)+theme_bw()+
  coord_sf(xlim=c(690,890),ylim=c(4870,5080))+
  xlab("Easting")+ylab("Northing")
# ggsave(filename=paste0(fig_dir,"log_spat_B_plot.png"),
#        plot=log_spatial_B_plot,height=12,width=12)

spatial_R_plot<-ggplot()+
  geom_sf(data=Rec_spat,aes(fill=R),col=NA)+
  facet_wrap(~Year)+
  scale_fill_viridis_c(name="Predicted Recruit \nDensity (kg/km^2)")+
  geom_sf(data=plot_map)+theme_bw()+
  coord_sf(xlim=c(690,890),ylim=c(4870,5080))+
  xlab("Easting")+ylab("Northing")
# ggsave(filename=paste0(fig_dir,"spat_R_plot.png"),
#        plot=spatial_R_plot,height=12,width=12)

log_spatial_R_plot<-ggplot()+
  geom_sf(data=Rec_spat,aes(fill=log(R)),col=NA)+
  facet_wrap(~Year)+
  scale_fill_viridis_c(name="Log Predicted Recruit \nDensity (kg/km^2)")+
  geom_sf(data=plot_map)+theme_bw()+
  coord_sf(xlim=c(690,890),ylim=c(4870,5080))+
  xlab("Easting")+ylab("Northing")
# ggsave(filename=paste0(fig_dir,"log_spat_R_plot.png"),
#        plot=log_spatial_R_plot,height=12,width=12)

spatial_m_plot<-ggplot()+
  geom_sf(data=m_spat,aes(fill=m),col=NA)+
  facet_wrap(~Year)+
  scale_fill_viridis_c(name="Predicted Instantaneous \nNatural Mortality")+
  geom_sf(data=plot_map)+theme_bw()+
  coord_sf(xlim=c(690,890),ylim=c(4870,5080))+
  xlab("Easting")+ylab("Northing")
# ggsave(filename=paste0(fig_dir,"spat_m_plot.png"),
#        plot=spatial_m_plot,height=12,width=12)

log_spatial_m_plot<-ggplot()+
  geom_sf(data=m_spat,aes(fill=log(m)),col=NA)+
  facet_wrap(~Year)+
  scale_fill_viridis_c(name="Predicted Log Instantaneous \nNatural Mortality")+
  geom_sf(data=plot_map)+theme_bw()+
  coord_sf(xlim=c(690,890),ylim=c(4870,5080))+
  xlab("Easting")+ylab("Northing")
# ggsave(filename=paste0(fig_dir,"log_spat_m_plot.png"),
#        plot=log_spatial_m_plot,height=12,width=12)

all_qIs$knotID<-all_qIs$knotID+1
qIs_spat<-left_join(all_qIs,gridded_bound_shape,by="knotID") %>% st_as_sf()

spatial_qI_plot<-ggplot()+
  geom_sf(data=qIs_spat,aes(fill=qI),col=NA)+
  facet_wrap(~source)+
  scale_fill_distiller(name="Commercial Size \nCatchability",palette="Spectral")+
  geom_sf(data=plot_map)+theme_bw()+
  coord_sf(xlim=c(690,890),ylim=c(4870,5080))+
  xlab("Easting")+ylab("Northing")
# ggsave(filename=paste0(fig_dir,"spat_qI_plot.png"),
#        plot=spatial_qI_plot,height=12,width=12)
# ggsave(filename=paste0(fig_dir,"spat_qI_plot_supp.png"),
#        plot=spatial_qI_plot,height=12,width=12)

sebdam_qIs_diffs<-sebdam_qIs
sebdam_qIs_diffs$qI<-(sebdam_qIs_diffs$qI-mean(sebdam_qIs_diffs$qI))/sd(sebdam_qIs_diffs$qI)

sehbam_qIs_diffs<-sehbam_qIs
sehbam_qIs_diffs$qI<-(sehbam_qIs_diffs$qI-mean(sehbam_qIs_diffs$qI))/sd(sehbam_qIs_diffs$qI)

all_qI_diffs<-rbind(sebdam_qIs_diffs,sehbam_qIs_diffs)
all_qI_diffs$knotID<-all_qI_diffs$knotID+1

qIs_diff_spat<-left_join(all_qI_diffs,gridded_bound_shape,by="knotID") %>% st_as_sf()

spatial_qI_diff_plot<-ggplot()+
  geom_sf(data=qIs_diff_spat,aes(fill=qI),col=NA)+
  facet_wrap(~source)+
  scale_fill_distiller(name="Standardized Commercial \nSize Catchability",palette="Spectral")+
  geom_sf(data=plot_map)+theme_bw()+
  coord_sf(xlim=c(690,890),ylim=c(4870,5080))+
  xlab("Easting")+ylab("Northing")
# ggsave(filename=paste0(fig_dir,"spat_qI_diff_plot.png"),
#        plot=spatial_qI_diff_plot,height=12,width=12)
# ggsave(filename=paste0(fig_dir,"spat_qI_diff_plot_supp.png"),
#        plot=spatial_qI_diff_plot,height=12,width=12)


