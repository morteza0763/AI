
##==*Packages===========================================
setwd("/media/vorteza/F/PMC/Codes/New_SD/")

# source("https://neuroconductor.org/neurocLite.R")
# neuro_install('aal', release = "current", release_repo = "github") 
# neuro_install('MNITemplate', release = "current", release_repo = "github") 
# remotes::install_github("TKoscik/nifti.io")
suppressPackageStartupMessages(
  {
    library(spatstat)
    library(MASS)
    # library(gpuR)
    # library(curand)
    # library(gpuRcore)
    library(splines)
    library(mvtnorm)
    library(tmvtnorm)
    library(LaplacesDemon)
    library(lattice)
    library(truncdist);library(sn)

    library(future)
    library(future.apply)
    library(doFuture)
    #library(pbdMPI)

    library(progressr)
    library(beepr)

    library(coda);library(ggmcmc)

    library(brainR);options(rgl.printRglwidget = F)
    library(nifti.io)
    library(neurobase)
    library(aal)
    library(MNITemplate)
    library(AnalyzeFMRI)
    library(tidyverse)

  }
)

#= *read brain templates===================
template_adress="atlas/Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.nii"
sch_mni152=readNIfTI(fname = template_adress)
sch_mni152.h=f.read.header(template_adress)

template <- sch_mni152#readNIfTI(system.file("MNI152_T1_2mm_brain.nii.gz", package="brainR")
                      #, reorient=FALSE)
dtemp=dim(template)
#======= *Load and clean data ===================
library(readxl)
SD_covar <- read_excel("SD_Coordinates_april.24.2018-merged_Tal.xlsx",
                       sheet = "demo_edit")
SD <- read_excel("SD_Coordinates_april.24.2018-merged_Tal.xlsx",
                 sheet = "Coordinates")
SD<-na.exclude(SD)
#--------------- Convert Tal to MNI --------------------
Tal=SD%>%select(X:Z,space)%>%
  filter(space=="Talairach")%>%
  select(X:Z)
write.table(Tal,"Tal.txt",row.names = F,col.names = F)
tal_to_mni=read.table("mni_convrted_spm.txt")
SD[SD$space=="Talairach",c("X","Y","Z")]<-tal_to_mni

#--------- Convert MNI to (normalized) XYZ -----------
# Converting back and form to MNI space
#origin <- c(x = 91, y = 126, z = 72)
#mni2xyz <- function(x) sweep(x, 2, as.double(origin), "+")

foci_xyz=SD%>%select(X,Y,Z)%>%
  as.matrix()%>%t()%>%
  xyz2ijk(L=sch_mni152.h)
foci_xyz<-t(foci_xyz[[1]]  );colnames(foci_xyz)<-c("X","Y","Z")
foci_xyz=tibble::as_tibble(foci_xyz)
# Normalizes----------------
#colin=readNIfTI("colin27_t1_tal_lin.nii")
foci_norm <- sweep(foci_xyz, 2, dim(template), "/") #dim(colin)
colnames(foci_norm)<-paste(c("X","Y","Z"),"norm",sep="_")

foci_norm_pp3<-pp3(foci_norm[,"X_norm"]
                   ,foci_norm[,"Y_norm"]
                   ,foci_norm[,"Z_norm"]
                   ,domain=box3())
plot(foci_norm_pp3)
save(foci_xyz,file="foci_xyz.Rdata")
save(foci_norm,file="foci_norm.Rdata")
save(foci_norm_pp3,file="foci_norm_pp3.Rdata")
library(tibble)
SD=add_column(SD,foci_norm,.after = "Z",.name_repair = "unique")
save(SD,file = "SD.Rdata")
#--------------- Join SD and SD_covar ---------------
library(tidyverse);library(ggplot2)
study_covar0=SD%>%
  #note:experiment is local variable
  select(experiment_author,space:modality,-experiment)%>%
  group_by(experiment_author)%>%
  slice_head(n=1)%>%ungroup()

study_covar1=SD_covar%>%
  select(c(3:4,6,8,10:15))

study_foci_count<-SD%>%count(experiment_author)

joint_covar<-left_join(study_covar0,study_covar1
                    ,by="experiment_author")
joint_covar<-left_join(joint_covar,study_foci_count
                    ,by="experiment_author")
#--- here add new covariates---------------------------
joint_covar<-joint_covar%>%
  rowwise()%>%
  mutate(
    sample_size=ifelse(recruiting_subjets_method=="case-control"
                       ,sum(n_SD,n_RS)
                       ,min(n_RS,n_SD))
  )%>%
  ungroup()%>%
  mutate(
    female_percent=Number_of_female_edited/sample_size
    ,hours_SD=relevel(factor(hours_SD),ref = "under24")
    )

joint_covar%>%
  ggplot(aes(reorder(experiment_author,-n),n))+
  geom_bar(stat = "identity")+theme_light()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

save(joint_covar,file="joint_covar.Rdata")
#------------ study specific response(PP3)----------------------
pts3_mni=by(SD,
            INDICES = SD$experiment_author
            ,FUN =function(sd){
              pts3_mni=pp3(sd$X,sd$Y,sd$Z,domain=box3())
              },simplify = T)
pts3_norm=by(SD,
            INDICES = SD$experiment_author
            ,FUN =function(sd){
              pts3_norm=pp3(sd$X_norm,sd$Y_norm,sd$Z_norm,domain=box3())
            },simplify = T)
names(pts3_mni)
##-------- Hyperframe (pts3_. and joint_covar)--------------
joint_SD_h=cbind.hyperframe(id=1:length(pts3_mni)
                            ,pts3_mni=pts3_mni
                            ,pts3_norm=pts3_norm
                            ,joint_covar)
save(joint_SD_h,file = "joint_SD_h.RData")

#=======*Load cleaned data===========================
{
  load("foci_xyz.Rdata")
  load("foci_norm.Rdata")
  load("foci_norm_pp3.Rdata")
  load("SD.Rdata")
  load("joint_covar.Rdata")
  load("joint_SD_h.RData")
  X_rep=joint_SD_h$pts3_norm

}


#---- *EDA--------------------------------------
#-----Explore Covariates------------------------
library(gtsummary)
### note that "expriment" column is not correct
### its value changes inside a study------------
covar_table=joint_covar%>%
  select(n,Age_edited,female_percent,sample_size
         ,space:recruiting_subjets_method
         ,-n_RS,-n_SD
         ,hours_SD
         ,Task_Category)%>%
  tbl_summary(
    statistic =list(
      # all_continuous()~"{mean} ({p25}-{p75})"
      # n_SD~"{median}({p25},{p75})",
      # n_RS~"{median}({p25},{p75})"
    ),
    label = list(n~"foci per study"),
    digits = list(all_continuous() ~ c(2, 0,0))
  )%>%
  bold_labels()
covar_table

covar_table%>%as_flex_table()%>%
  flextable::save_as_docx(path = "cavar_table.docx")
#--------------------------------
docx=officer::read_docx("cavar_table.docx")
docx=SD%>%count(experiment)%>%
  rename("Number of foci"=n)%>%
  flextable::flextable()
  officer::body_add_table(x = docx)
  
expr_count=SD%>%
  select(experiment)

#-----Explore response -------------------------
## foci per study------------------------------
docx=joint_covar%>%
  ggplot(aes(reorder(experiment_author,-n),n))+
  geom_bar(stat = "identity")+theme_light()+
  xlab("Study")+ylab("# foci")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))%>%
  officer::body_add_gg(x=docx)

joint_covar%>%
  pull(n)%>%
  psych::describe()
## foci by covariates-----------------------
joint_covar%>%
  dplyr::select(n,hours_SD)%>%
  mutate(sqrt.n=sqrt(n))%>%
  tbl_summary(by=hours_SD,
              statistic = all_continuous()~"{mean} ({sd})",
              digits = n~c(2,2))%>%
  add_p()%>%
  modify_spanning_header(starts_with("stat_") ~ "**Duration of SD**")%>%
  as_flex_table()%>%
  flextable::save_as_docx(path = "foci_by_hoursSD.docx")


joint_covar%>%
  mutate(modality=replace(modality,
    modality=="rsfMRI"|modality=="tfMRI","fMRI"))%>%
  select(n,modality)%>%
  mutate(sqrt.n=sqrt(n))%>%
  tbl_summary(by=modality,
              statistic = all_continuous()~"{mean} ({sd})",
              digits = n~c(2,2))%>%
  add_p()

##=== plots================================

pl_hours_SD<-joint_covar%>%
  dplyr::select(n,hours_SD)%>%
  mutate(sqrt_n0=sqrt(n))%>%
  group_by(hours_SD)%>%
  summarise(sqrt_n=mean(sqrt_n0),sd=sd(sqrt_n0))%>%
  ggplot(aes(fill=hours_SD,x=hours_SD,y=sqrt_n))+
  geom_col()+
  geom_errorbar(aes(ymax=sqrt_n+1.96*sd/sqrt(29)
                    ,ymin=sqrt_n-1.96*sd/sqrt(29))
                ,width=.2)+
  ylab(expression(
    sqrt("# foci per study")
  ))+
  scale_fill_discrete(name="SD duration")+
  xlab("SD duration")+
  theme_light()+
  theme(legend.position = "top")

  
pl_Age=joint_covar%>%
  ggplot(aes(x=Age_edited,y=sqrt(n)))+
  geom_point(size=2)+ylab(expression(
    sqrt("# foci per study")
    ))+
  xlab("Age (years)")+
  geom_smooth(inherit.aes = F
              ,data = joint_covar%>%
                filter(Age_edited<31)
              ,aes(x=Age_edited,y=sqrt(n))
              ,method = lm
              )+
  theme_light()

joint_covar%>%
  filter(Age_edited<31)->xxx
cor.test(sqrt(xxx$n),xxx$Age_edited
         ,method = "spearman")

pl_sample_size=joint_covar%>%
  ggplot(aes(x=sample_size,y=sqrt(n)))+
  geom_point(size=2)+
  geom_smooth(
    method = lm
    )+
  ylab(expression(
    sqrt("# foci per study")
  ))+
  xlab("Sample size per study")+
  theme_light()
cor.test(sqrt(joint_covar$n),joint_covar$n_SD
         ,method = "spearman")



c3=cor.test(sqrt(joint_covar$n)
            ,joint_covar$female_percent
            ,method = "spearman")
c3=paste("rho=",round(c3$estimate,2)
         ,",","p-value=",round(c3$p.value,2))

pl_female<-joint_covar%>%
  ggplot(aes(x=female_percent,y=sqrt(n)))+
  geom_point(size=2)+
  geom_smooth(method=lm)+ylab(expression(
    sqrt("# foci per study")
  ))+
  xlab("percentage of females (%)")+theme_light()
geom_text(x=.25,y=7.5,size=4
            ,colour="black"
            ,fontface="plain"
            ,parse=F
            ,label=c3
            ,inherit.aes = F
            )+
  theme_light()
library(gridExtra)
all_pl<-marrangeGrob(grobs = list(pl_hours_SD,pl_Age
                          ,pl_sample_size,pl_female)
             ,ncol = 2,nrow = 2)

## 3d brain plot ===================================
##=======overall 3d plot ===========================
contour3d(template, level = 0# 4500 for aal
          ,alpha = 0.3, draw = TRUE)
im.area1<-as.data.frame(foci_norm_pp3)
rgl.points(x = im.area1$x*dtemp[1],
           y = im.area1$y*dtemp[2],
           z = im.area1$z*dtemp[3],size = 5,alpha=.8,
           col="darkred")

# spheres3d(x = im.area1$x*dtemp[1],
#           y = im.area1$y*dtemp[2],
#           z = im.area1$z*dtemp[3],radius = 4,alpha=.4,
#           col="orange")
# spheres3d(x = im.area1$x*dtemp[1],
#           y = im.area1$y*dtemp[2],
#           z = im.area1$z*dtemp[3],radius = 1,alpha=.8,
#           col="black")

# text3d(x = dtemp[1]/2, y = dtemp[2]/2, z = dtemp[3] * 0.98, text = "Top")
# text3d(x = dtemp[1] * 0.98, y = dtemp[2]/2, z = dtemp[3]/2, text = "Right")

##==== Explore interaction range=====================
int_range_perc=sapply(X_rep,function(.x) {
  sum(has.close(.x,2*.05))/npoints(.x)
  sum(has.close(.x,2*.05))
})
summary(int_range_perc)
boxplot(int_range_perc)

sum(int_range_perc)/sum(sapply(X_rep,npoints))
##======= **Selected study plot =====================
library(ggseg3d)
library(ggsegSchaefer)
require(plotly)

# mfrow3d(byrow = T,1,4,sharedMouse = T)
# for(i in 1:29){
#   if(i%%4==0){
#     open3d()
#     mfrow3d(byrow = T,1,4,sharedMouse = T)
#   }
#     
#   next3d()
#   pl=joint_SD_h$pts3_norm[i]
#   rgl.points(coords(pl[[1]]),col="red",size=8)
#   title3d(paste(names(pl),"-",i))
# }

library(flextable)
j=c(10,16,25)
joint_SD_h[j,c("experiment_author","n",
               "sample_size","female_percent"
               ,"Age_edited","hours_SD")]%>%
  as.data.frame()%>%
  flextable()%>%
  bold(part = "header")%>%
  save_as_docx(path="selected_studies.docx")
  
p_fig2<-vector("list",length(j))
names(p_fig2)<-joint_SD_h$experiment_author[j]

# j=j[1]
someData = schaefer17_3d %>% 
  filter(surf == "inflated" & hemi == "right") %>% 
  unnest(ggseg_3d) %>% 
  ungroup() %>% 
  select(region)
c=1
for(jj in j){
  p_fig2[[c]]<- ggseg3d(
    .data = someData%>%
      slice(1)%>% mutate(region=NA, mycol="lightgray")
    ,atlas = "schaefer17_3d"
    ,hemisphere = "right"
    ,colour="mycol",na.alpha = .3
  )%>%
    # add_glassbrain(hemisphere = "right",colour = "lightgray"
    #                 ,opacity = .3)%>%
    remove_axes()%>%
    pan_camera("right lateral")%>%
    add_markers(inherit = T,
                x = joint_SD_h$pts3_mni[[jj]]$data$x,
                y = joint_SD_h$pts3_mni[[jj]]$data$y,
                z = joint_SD_h$pts3_mni[[jj]]$data$z
                # ,color=I("red")
                , type = 'scatter'
                , mode = 'markers'
                ,size=1,name=jj
                ,showlegend = FALSE
                ,marker=list(
                  # cauto=F,cmin=-60,cmax=60,
                  color=joint_SD_h$pts3_mni[[jj]]$data$x
                  ,colorscale="RedBlue"
                  ,colorbar=list(title="X-axis\n(MNI)"
                                 ,width=2
                                 ,tickfont  = list(size = 17,face="bold"))
                  ,line=list(width=2,color="black")
                  ,name=jj,
                  opacity=.6,size=4))%>%
    add_annotations(
      text =paste("foci count=",joint_SD_h$n[jj])
      ,showarrow = F
      , xref='paper', yref='paper',y=.1
      ,yanchor="bottom",valign="bottom"
      ,font=list(size=17)
      )%>%
    layout(title=joint_SD_h$experiment_author[jj]
           # ,width=400,height=300
           )

  c=c+1
  
}
c=1
p_fig2[[3]]%>%
  orca(file =paste0(joint_SD_h$experiment_author[j[3]],".pdf" )
       ,format = "pdf"#,scale=3
       ,width=400
       ,height=300
  )

# doc <- htmltools::tagList(
#   div(p_fig2[[1]], style = "float:left"),#;width:20%;"),
#   div(p_fig2[[2]],style = "float:left"),#;width:20%;"),
#   div(p_fig2[[3]], style = "float:left"),#;width:20%;")
# )
# 
# htmltools::save_html(html = doc, file = "widgets.html")

        
### classic codes=====
# {
#   rgl.bg(color=c("white","white"))
#   par3d(mouseMode="trackball",windowRect=c(  67,83,1360,796))
#   rgl.viewpoint( theta = 0, phi = 15, fov = 60, zoom = 5,
#                  scale = par3d("scale")
#                  , interactive = TRUE,
#                  type = "modelviewpoint")
#   mfrow3d(byrow = T,1,4,sharedMouse = T)
#   contour3d(template, level = 4500,alpha = 0.3, draw =T ,add = T)
#   pch3d(x = im.area1$x*dtemp[1],
#         y = im.area1$y*dtemp[2],
#         z = im.area1$z*dtemp[3],radius = 2,alpha=.8,pch = 16
#         ,col="darkred")
# 
# 

#   selected.studies=which(
#     sapply(
#       joint_SD_h$pts3_norm,npoints)>30)[2:4]
#   names(joint_SD_h$pts3_norm)[selected.studies]
# 
#   for( j in selected.studies){
#     cat("plot study",j,"\n")
#     next3d()
#     # rgl.bg(color=c("white","white"))
#     # par3d(mouseMode="trackball"
#     #       ,windowRect=c(  67,83,1360,796))
#     # rgl.viewpoint( theta = 0, phi = 15, fov = 60, zoom = 5,
#     #                scale = par3d("scale"), interactive = TRUE,
#     #                type = "modelviewpoint")
#     contour3d(template, level = 4500,alpha = 0.3,add=F,draw = TRUE)
#     pch3d(x = joint_SD_h$pts3_norm[[j]]$data$x*dtemp[1],
#           y = joint_SD_h$pts3_norm[[j]]$data$y*dtemp[2],
#           z = joint_SD_h$pts3_norm[[j]]$data$z*dtemp[3]
#           ,radius =3,alpha=1,pch = 16
#           ,col="red")
# 
#   }
# }
# 
# rgl.postscript(filename = "seplected_plot.pdf",fmt = "pdf")
#writeWebGL(filename = "selected studies 3d view.html")

#======* MODELING=======================================
#======*Load replicated volume==========================
#### for update volume please run shpere_volume.R from batch
{
  resol=c(91,109,91)
  load(paste("Area3d_rep_grid_0.05_res",resol[1],".RData",sep=""))
  Area3d_rep_2e6_r0.05<-Area3d_rep_grid_0.05_res91
  r=0.05;win=bounding.box3(foci_norm_pp3);ncycles=1
}

#======*Load SPP Rnd Gen ===============================
source("3dpp.R")  # resol also set to '91 109 91' via '3dpp.R'

#======*g_rep (Replicated un normalized AI)=============
##----- test initial value for g_rep  -----------------
model.matrix(~hours_SD+(1|id),data = joint_SD_h)


Ltheta<-c(Lbeta0=1.60,Lbeta_h=log(1))#,Lbeta_A=log(1))
Ltheta_local<<-c(Lb_x=0,Lb_y=0,Lb_z=0)
Leta<-c(Leta0=1.36,Leta_h=log(1))
Ltheta;Leta

int=rep(1,nrow(joint_SD_h))
h=as.numeric(joint_SD_h$hours_SD)-1
B_fix=cbind(int,h)
B_rnd=cbind(int)

alphi=c(-.404,-.404);blphi=c(2,2)
nu=10;a=2
sigma_Lphi=rhuangwand(nu,rep(a,2))#rinvwishart(2,100*diag(2))
sigma_Lphi

cov2cor(sigma_Lphi)

Lphi12=rtmvnorm(n = nrow(joint_SD_h),mean = c(0,0),sigma =sigma_Lphi
                ,lower = alphi,upper =blphi )

xt=seq(alphi[1],blphi[1],length.out = 20)
plot(xt,dtrunc(xt,"norm",mean=0,sd=sqrt(sigma_Lphi[1,1])
               ,a = alphi[1],b=blphi[1]),type="l",ylim=c(0,1),ylab="")
lines(xt,dtrunc(xt,"norm",mean=0,sd=sqrt(sigma_Lphi[2,2])
                ,a = alphi[1],b=blphi[1]),type="l",col="blue")

E_fix<-B_fix
E_rnd<-B_rnd
# win=bounding.box3(foci_norm_pp3)

#future:::ClusterRegistry("stop")
#plan(muticore,workers=3)

# system.time({
#   Area3d_rep<-future_lapply(joint_SD_h$pts3_norm,function(pts3_i){
#     sphere_volume_mc(n=2e6,S.centers = pts3_i,S.r = .05)
#   }
#   )
# })

X_rep=joint_SD_h$pts3_norm
g_rep=function(X_rep
               ,B_fix,B_rnd,Ltheta
               ,sigma_Lphi,Lphi12=NULL
               ,E_fix,E_rnd,Leta
               ,r=.05,log=TRUE){

  if(is.null(Lphi12)) {
    Lphi12=rtmvnorm(n = nrow(joint_SD_h),mean = c(0,0),sigma =sigma_Lphi
                    ,lower = alphi,upper =blphi )
  }

  if(!"Area3d_rep_2e6_r0.05" %in% ls(envir = .GlobalEnv)){
    Area3d_rep_2e6_r0.05<-future_lapply(X_rep,function(pts3_i){
      sphere_volume_mc(n=1e6,S.centers = pts3_i,S.r=r)
    }
    )
  }

  n_rep=unlist(anylapply(X_rep,npoints))
  #--- Fist order trend for for homogeneous AI----------
  Lbeta_rep<-n_rep*(B_fix%*%Ltheta+B_rnd*Lphi12[,1])

  #--- Fist order trend for for inhomogeneous AI--------
  # beta_i=function(x,y,z,i,log=T){
  #   Lbeta_ii<-Lbeta_rep[i,]/n_rep[i]+
  #     Ltheta_local["Lb_x"]*x+
  #     Ltheta_local["Lb_y"]*y+
  #     Ltheta_local["Lb_z"]*z
  #   out=ifelse(log==T,Lbeta_ii,exp(Lbeta_ii))
  #   names(out)=ifelse(log==T,paste("Lbeta",i,sep="")
  #                     ,paste("beta",i,sep=""))
  #   return(out)
  #
  # }
  #
  # ii<<-1
  # L_beta_rep_inhom=lapply(X_rep,function(X_i){
  #   beta_i_inhom=mapply(beta_i,X_i$data$x,X_i$data$y,X_i$data$z,ii)
  #   return(beta_i_inhom)
  #   ii<<-ii+1
  # })

  #--- Higher order trend for AI------------------------
  LETA_rep<-E_fix%*%Leta+E_rnd*Lphi12[,2]

  C_rep=numeric(length(n_rep))
  names(C_rep)<-names(n_rep)
  for(i in 1:length(n_rep)){
    C_rep[i]=Area3d_rep_2e6_r0.05[[i]]/(4/3*pi*r^3)-n_rep[i]

  }
  high_order_rep<- (-C_rep)*LETA_rep

  #--- replicated g_rep----------------------------
  finalout=sum(Lbeta_rep)+sum(high_order_rep)
  #----in case of inhom AI-------------------------
  # out=sum(sapply(L_beta_rep_inhom,sum))-sum(high_order_rep)

  #---------------- final output-----------------
  names(finalout)<-"ln"
  if(log==FALSE){
    finalout=exp(finalout); names(finalout)<-"exp"
  }
  return(finalout)
}


g_rep(X_rep = X_rep
      ,B_fix = B_fix,B_rnd = B_rnd,Ltheta = Ltheta
      ,sigma_Lphi = sigma_Lphi,Lphi12 = NULL
      ,E_fix = E_fix,E_rnd = E_rnd,Leta = Leta
      ,r=.05,log=TRUE)



############## Replicated AI ########################
assign("Ltheta_local",Ltheta_local,envir = .GlobalEnv)
handlers("progress")

rags_rep_AreaInter.3d=function(X_rep=joint_SD_h$pts3_norm
                               ,B_fix,B_rnd,Ltheta
                               ,sigma_Lphi,Lphi12
                               ,E_fix,E_rnd,Leta
                               ,r=.05
                               ,win=box3(),ncycles=ncycles
                               ,log=TRUE){

  n_rep<-unlist(anylapply(X_rep,npoints))
  Lbeta_rep<-n_rep*(B_fix%*%Ltheta+B_rnd*Lphi12[,1])

  ###----------please unblock if beta_f_rep is blocked

  # beta_i<-function(x,y,z,myrep){
  #   Lbeta_ii<-
  #     Lbeta_rep[myrep,]/n_rep[myrep]#+
  #   #   Ltheta_local["Lb_x"]*x+
  #   #   Ltheta_local["Lb_y"]*y+
  #   #   Ltheta_local["Lb_z"]*z
  #
  #   out=exp(Lbeta_ii)
  #   names(out)<-paste("beta",myrep,sep="")
  #   return(out)
  #
  # }
  beta_f_rep=mapply(# a function that return function!
    function(myrep0){
      function(x,y,z,myrep=myrep0){
        Lbeta_ii<-
          Lbeta_rep[myrep,]/n_rep[myrep]#+
        #   Ltheta_local["Lb_x"]*x+
        #   Ltheta_local["Lb_y"]*y+
        #   Ltheta_local["Lb_z"]*z

        out=exp(Lbeta_ii)
        names(out)<-paste("beta",myrep,sep="")
        return(out)

      }
    },myrep0=1:length(X_rep)
  )

  LETA_rep<-E_fix%*%Leta+E_rnd*Lphi12[,2]

  sim_rep<-vector("list",length(X_rep))
  # cat("\nsimulating replicated AI\n")

  with_progress(
    {
      pb_AI=progressor(along = 1:length(X_rep),auto_finish = F)
      sim_rep<-future_lapply( 1:length(X_rep)
                              ,future.seed = TRUE
                              ,FUN = function(i){
        pb_AI(sprintf("rep AI=%g",i))
        vlmax<-c(exp(Lbeta_rep[i,]/n_rep[i])*exp(LETA_rep[i,])
                 ,LETA_rep[i,]/(4/3*pi * r^3)) # in case of homogeneous AI
        #assign("myrep",i,envir = .GlobalEnv)

        return(
          ragsAreaInter.3d(beta = beta_f_rep[[i]]#beta_i,myrep=myrep
                           ,eta = exp(LETA_rep[i,])
                           ,r=r,win = win
                           ,ncycles = ncycles
                           ,vlmax=vlmax
          )
        )# end return
      })#end future.lapply

    }) # end bp_AI
  names(sim_rep)<-names(n_rep)
  out_rep=hyperframe(foci_norm_pp3=X_rep,simulate=sim_rep)
  return(out_rep)
}

future:::ClusterRegistry("stop")
plan(multisession,workers=3)
registerDoFuture()
handlers(handler_void())
system.time({
  simulated_3dRAI=rags_rep_AreaInter.3d(X_rep=joint_SD_h$pts3_norm
                                        ,B_fix = B_fix,B_rnd = B_rnd,Ltheta = Ltheta
                                        ,sigma_Lphi =sigma_Lphi,Lphi12 = Lphi12
                                        ,E_fix = E_fix,E_rnd = E_rnd,Leta = Leta
                                        ,win = box3(),ncycles = 1
                                        ,r=.05,log=TRUE)
}
)

sum(sapply(simulated_3dRAI$simulate,npoints))


#plot(im.foci.ashape,transparency = 0.2,col="gold")0
#rgl.points(as.matrix(foci_norm_pp3),pch3=2,size = 5,col="orange")
simulate_rep_AI=do.call(rbind,sapply(simulated_3dRAI$simulate,as.matrix.ppx))
rgl.points(simulate_rep_AI,pch3=20,size = 5,col="red")

simpp3=pp3(simulate_rep_AI[,1],simulate_rep_AI[,2],simulate_rep_AI[,3],domain=box3())
plot(pcf3est(simpp3),xlim=c(0,.3),ylim=c(-1,8))
plot(pcf3est(foci_norm_pp3),add=T,col="blue")

foci_xyz_pp3<-pp3(foci_xyz$X,foci_xyz$Y,foci_xyz$Z,
                  domain=box3(
                    range(foci_xyz$X),range(foci_xyz$Y),
                    range(foci_xyz$Z)
                  ))
plot(pcf3est(foci_xyz_pp3))

# =====*Initialize Adaptive Metropolis-within-Gibbs Algorithm ===============
{
  ##----- test initial value ---------
  alphi=log(c(.2,.5))
  blphi=log(c(4,2))

  pr.lower=c(Lbeta0=log(12),Lbeta_h=log(.8),Leta0=log(1*1/.8*1/.5),Leta_h=log(.8))
  pr.upper=c(Lbeta0=log(12+4),Lbeta_h=log(1.2),Leta0=log(1*1/.8*1/.5+5),Leta_h=log(1.2))

  Ltheta<-c(Lbeta0=(pr.lower[[1]]+pr.upper[[1]])/2#
            ,Lbeta_h=0)#,Lbeta_A=log(1))
  Ltheta_local<<-c(Lb_x=0,Lb_y=0,Lb_z=0)
  Leta<-c(Leta0=(pr.lower[[3]]+pr.upper[[3]])/2
        ,Leta_h=0)
  Ltheta;Leta

  int=rep(1,nrow(joint_SD_h))
  h=as.numeric(as.factor(joint_SD_h$hours_SD))-1
  B_fix=cbind(int,h)
  B_rnd=cbind(int)

  nu=10
  a=c(2,5)
  sigma_Lphi= rhuangwand(nu,a)#rhuangwand(20,rep(1e4,2))#rinvwishart(2,100*diag(2))
  sigma_Lphi
  cov2cor(sigma_Lphi)
  Lphi12=rtmvnorm(n = nrow(joint_SD_h),mean = c(0,0)
                  ,sigma =sigma_Lphi
                  ,lower = alphi,upper =blphi )
  hist(
    exp(
      Lphi12[,1])
    ,ylim=c(0,20)
    ,col=adjustcolor(2,alpha.f=.3))
  hist(
    exp(
      Lphi12[,2])
    ,add=T
    ,col=adjustcolor(4,alpha.f=.3))


  xt=seq(min(alphi)
         ,max(blphi),length.out = 1000)
  plot(xt,dtrunc(xt,"norm",mean=0,sd=sqrt(sigma_Lphi[1,1])
                 ,a = alphi[1],b=blphi[1]),type="l"
       ,ylim=c(0,2))
  lines(xt,dtrunc(xt,"norm",mean=0,sd=sqrt(sigma_Lphi[2,2])
                  ,a = alphi[2],b=blphi[2]),type="l",col="blue")

  E_fix<-B_fix
  E_rnd<-B_rnd
  win=box3()#ounding.box3(foci_norm_pp3)
  ncycles=2
  #-------------------MCMC inits-------------------------
  mcmc.len=5e4
  PAR<- matrix(NA,mcmc.len,length(Ltheta)+length(Leta))
  PAR<-cbind(as.hyperframe(PAR),sigma_Lphi)
  names(PAR)<-c(names(Ltheta),names(Leta),"sigma_Lphi")
  head(PAR)

  par=PAR
  head(par)
  ##---- set initial values--------------------
  par$Lbeta0[1]<-Ltheta[1];par$Lbeta_h[1]<-Ltheta[2]
  par$Leta0[1]<-Leta[1];par$Leta_h[1]<-Leta[2]
  par$sigma_Lphi[[1]]<-sigma_Lphi

  #par<-par[1:mcmc.len]

  Lphi1=matrix(NA,mcmc.len,nrow(joint_SD_h)); colnames(Lphi1)<-row.names(joint_SD_h)
  #Lphi1=rbind(Lphi1[1:i,],matrix(NA,mcmc.len,nrow(joint_SD_h)))
  Lphi2=matrix(NA,mcmc.len,nrow(joint_SD_h)); colnames(Lphi2)<-row.names(joint_SD_h)
  #Lphi2=rbind(Lphi2[1:i,],matrix(NA,mcmc.len,nrow(joint_SD_h)))

  Y_t=vector(mode = "list",length = mcmc.len)
  #Y_t<-Y_t[1:mcmc.len]
  #handlers("progress")
  # pb<- txtProgressBar(title ="Normalizing constant" ,1, 25, style=3)
  par(mfrow=c(4,2),mar=c(4.1, 4.1, 3.1, 3.5),pch=16, xpd=TRUE)
  i=1;accept.counter=0; accept.status=F;r_m=0
  ##---- Proposal initial value----------------
  lambda<-as.data.frame(par)
  lambda[1,]<-c(1,1,2,1)/2#rep(2,ncol(lambda))
  lambda[1,]<-c(2,2,2,2)/1.5
  # lambda<-rbind(lambda
  #               ,matrix(NA,mcmc.len,ncol(lambda)
  #                       ,dimnames = list(NULL,colnames(lambda)))
  #               )

  mu=as.data.frame(par)
  mu[1,]<-mu[1,]+.005
  # mu<-rbind(mu
  #               ,matrix(NA,mcmc.len,ncol(mu)
  #                       ,dimnames = list(NULL,colnames(mu)))
  #               )

  matcor_kk<-0;k=0
  #i0=1.2,e=.95
  gain<-function(i,i0=1.2,e=.95){#e=.95
    if(i0<=1) stop("i0 must be i0> 1")
    if(e<=.5 | e>1) stop("e must be: .5<e<=1")
    i0/max(i0,i^e)

  }

  gain=Vectorize(gain)
  par(mfrow=c(1,1),mar=c(4.1, 4.1, 4.1, 3.5),pch=16, xpd=TRUE)

  plot(1:500,gain(1:500),t="l",ylim=c(0,1))
  lines(1:500,gain(1:500,i0=1.2))
  lines(1:500,gain(1:500,i0=2,e=.95),col="red")
  lines(1:500,gain(1:500,i0=1.2,e=.55),col="blue")
  lines(1:500,gain(1:500,i0=2.5,e=.85),col="orange")
  alpha_star=.50
#Lbeta0=0

  # pr.lower=c(Lbeta0=1.2,Lbeta_h=-.6
  #            ,Leta0=1,Leta_h=-.5)
  # pr.upper=c(Lbeta0=4,Lbeta_h=.6
  #            ,Leta0=7,Leta_h=.6)
  # names(pr.lower)<-colnames(par)
  # names(pr.upper)<-colnames(par)

  ##----prior options-----
  prior.lower=pr.lower
  prior.upper=pr.upper

  matcor.prior=diag(length(pr.lower))
  colnames(matcor.prior)=rownames(matcor.prior)=names(pr.lower)

  diag(matcor.prior)<-50* c(Lbeta0=1,Lbeta_h=1
                            ,Leta0=1,Leta_h=1)
  matcor.prior[upper.tri(matcor.prior)]<-matcor.prior[lower.tri(matcor.prior)]<-0
  matcor.prior["Lbeta0","Leta0"]<-matcor.prior["Leta0","Lbeta0"]<- -40#0
  cov2cor(matcor.prior)
  # det(matcor.prior);solve(matcor.prior)
  #
  # o=Matrix::nearPD(matcor.prior)
  # cov2cor(o$mat)
  # det(o$mat)
  # matcor.prior=o$mat
  matcor<-diag(4)#matcor.prior/50
  colnames(matcor)=row.names(matcor)=colnames(matcor.prior)
  matcor["Lbeta0","Leta0"]<-matcor["Leta0","Lbeta0"]<- -.5
  cov2cor(matcor)

  show_plot=TRUE
}

rbind(pr.lower,par=par[1,1:4,drop=T,strip=F],pr.upper)
alphi;blphi
#====*Parallel setting =================
handlers(handler_void())
ncycles=1

future:::ClusterRegistry("stop")
plan(multicore,workers=4)
registerDoFuture()

#======*Start Adaptive Metropolis-within-Gibbs Algorithm ===============
repeat{

  #----1: Draw v from a proposal distribution Q(par[i,], v)--------------
  if(is.nan(r_m)){
    k=k
  }else{
    k=sample(x = names(par),size = 1,replace = T)
  }

  if(k=="sigma_Lphi"){
    # if(any(abs(par[i,5,T])>10)){
    #   par[i,5]$sigma_Lphi[[1]]=par[i,5,T]/max(abs(par[i,5,T]))
    #   }
    v_k=rhuangwand(nu,a)
  }else{
    matcor_kk=lambda[i,k]*diag(matcor)[k]
    v_k=rtrunc(1,mean=par[i,k,T],sd=sqrt(matcor_kk),spec = "norm"
               ,a=pr.lower[k],b=pr.upper[k])
  }

  v=par[i,];v[,k]<-hyperframe(v_k)
  #----2a: Estimate the normalizing constant ratio R(par[i,],v)----------
  sigma_Lphi<-par[i,"sigma_Lphi",T]
  Lphi12<- rtmvnorm(nrow(joint_SD_h),c(0,0),sigma = sigma_Lphi
                    ,lower =alphi
                    ,upper = blphi)
  Ltheta<-unlist(par[i,c("Lbeta0","Lbeta_h"),T])
  Leta<-unlist(par[i,c("Leta0","Leta_h"),T])
  #----------------------------------------------------------------------#

  #if(accept.status==T | i==1){
  cat("simulating auxillary AI for iter",i,"\n")
  # y_t<-foreach(mc_iter =1:25,.export = c("Ltheta_local","n_rep","Lbeta_rep")) %dopar% {
  # with_progress({
  # pb_mc<-progressor(along = 1:1,auto_finish = F)
  y_t<-lapply(1:1,function(mc_iter){ #future_lapply
    # pb_mc(sprintf("mc_iter=%g",mc_iter))
    return(rags_rep_AreaInter.3d(X_rep=joint_SD_h$pts3_norm # could be spatial covariate or original data=foci_norm_pp3
                                 ,B_fix = B_fix,B_rnd = B_rnd,Ltheta = Ltheta
                                 ,sigma_Lphi = sigma_Lphi,Lphi12 = Lphi12
                                 ,E_fix = E_fix,E_rnd =E_rnd ,Leta = Leta
                                 ,win = win,ncycles = ncycles
                                 ,r=r,log=TRUE)
    )
  }) #end future_lapply

  # }) # end bp_mc
  # }# end if

  sigma_Lphi_v<-v[,"sigma_Lphi",T]
  Lphi12_v=rtmvnorm(nrow(joint_SD_h),mean = c(0,0),sigma =sigma_Lphi_v
                    ,lower = alphi
                    ,upper = blphi
  )

  g_y_v=unlist(lapply(y_t,function(sim){
    g_rep(X_rep=sim$simulate
          ,B_fix,B_rnd,Ltheta = unlist(v[,c("Lbeta0","Lbeta_h"),T])
          ,sigma_Lphi=sigma_Lphi_v,Lphi12 =Lphi12_v
          ,E_fix,E_rnd,Leta = unlist(v[,c("Leta0","Leta_h"),T])
          ,r=r,log=TRUE)

  }))
  g_y_par=unlist(lapply(y_t,function(sim){
    g_rep(X_rep=sim$simulate
          ,B_fix,B_rnd,Ltheta = unlist(par[i,c("Lbeta0","Lbeta_h"),T])
          ,sigma_Lphi = sigma_Lphi,Lphi12=Lphi12
          ,E_fix,E_rnd,Leta = unlist(par[i,c("Leta0","Leta_h"),T])
          ,r=r,log=TRUE)

  }))

  #R_m=mean(exp(g_y_v-g_y_par))
  #gdif=g_y_v-g_y_par
  #lR_m=max(gdif)-log(10)
  #lR_m=max(gdif)+log(sum(exp(gdif-max(gdif))))-log(10)

  lR_m=g_y_v-g_y_par
  #----2b:Hastings ratio-------------------
  lr_m1=
    #(1/R_m)*
    #exp(
    -(lR_m)+
    g_rep(X_rep=joint_SD_h$pts3_norm
          ,B_fix,B_rnd,Ltheta = unlist(v[,c("Lbeta0","Lbeta_h"),T])
          ,sigma_Lphi = sigma_Lphi_v,Lphi12 = Lphi12_v
          ,E_fix,E_rnd,Leta = unlist(v[,c("Leta0","Leta_h"),T])
          ,r=r,log=TRUE)-g_rep(X_rep=joint_SD_h$pts3_norm
                               ,B_fix,B_rnd,Ltheta = unlist(par[i,c("Lbeta0","Lbeta_h"),T])
                               ,sigma_Lphi = sigma_Lphi,Lphi12 = Lphi12
                               ,E_fix,E_rnd,Leta = unlist(par[i,c("Leta0","Leta_h"),T])
                               ,r=r,log=TRUE)+
    dtmvnorm(x = unlist(v[,1:4,T]),mean = unlist(par[1,1:4,T]),sigma = matcor.prior,log = T
             ,lower=prior.lower
             ,upper=prior.upper)+dhuangwand(v[,5,T],nu,a,A=rep(1e6,2),log = T)-
    (dtmvnorm(unlist(par[i,1:4,T]),mean = unlist(par[1,1:4,T]),sigma = matcor.prior,log = T
              ,lower=prior.lower
              ,upper=prior.upper)+dhuangwand(par[i,5,T],nu,a,A=rep(1e6,2),log = T))
  #)#end r_m1
  if(k=="sigma_Lphi"){
    lr_m2=#exp(
      dhuangwand(v[,k,T],nu,a,A=rep(1e6,2),log = T)-
      dhuangwand(par[i,k,T],nu,a,A=rep(1e6,2),log = T)
    #)
  }else{
    lr_m2=#exp(
      dtrunc(v[,k,T],mean = par[i,k,T],sd = sqrt(matcor_kk),spec = "norm"
             ,a=pr.lower[k]
             ,b=pr.upper[k]
             ,log=T)-dtrunc(par[i,k,T],mean = v[,k,T],sd = sqrt(matcor_kk),spec = "norm"
                            ,a=pr.lower[k]
                            ,b=pr.upper[k]
                            ,log=T)
    #)
  }

  r_m=exp(lr_m1+lr_m2)

  #---- Accept/Reject-------------------------
  u=runif(1)
  if(!is.nan(r_m)){
    if(u<=min(1,r_m)){
      accept.status=T
      par[i+1,]<- v
      # PAR[accept.counter,]<-v
      accept.counter<-accept.counter+1
      Lphi1[i+1,]=Lphi12_v[,1]
      Lphi2[i+1,]=Lphi12_v[,2]
      Y_t[[i+1]]<-y_t[[1]]$simulate
    }else{
      accept.status=F
      par[i+1,]<-par[i,]

      Y_t[[i+1]]<-Y_t[[i]]
      Lphi1[i+1,]<-Lphi12[,1]
      Lphi2[i+1,]<-Lphi12[,2]
    }
    ##----- Print information-----
    cat("\n iter",i,",r_m=",r_m
        ,",accept.counter=",accept.counter
        ,",accept.ratio=",accept.counter/i
        ,"\n","var[,",k,"]=",matcor_kk
        ,"\n","v=(",paste(names(v),round(unlist(v[,1:4,T]),3),sep="=",collapse = ", "),")"
        ,"\n----------------------------\n")
    if(accept.counter==mcmc.len ) break;

    ##----- Adaptive part---------
    if(k!="sigma_Lphi"){
      lambda[i+1,k]<-exp(
        log(lambda[i,k])+gain(i+1)*(min(1,r_m)-alpha_star)
      )

      lambda[i+1,!colnames(lambda)%in%k]<-lambda[i,!colnames(lambda)%in%k]
      mu[i+1,]<-mu[i,]+gain(i+1)*(unlist(par[i+1,1:4,T])-mu[i,])
      matcor<-matcor+ gain(i+1)*( t(unlist(par[i+1,1:4,T])-mu[i,])%*%(as.matrix(unlist(par[i+1,1:4,T])-mu[i,])) - matcor )
    }else{
      lambda[i+1,]<-lambda[i,]
      mu[i+1,]<-mu[i,]
    }
    
    tracedata=cbind(
      as.data.frame(par[1:i,1:4,T])
      ,sigma_Lphi1=sapply(par[,5]$sigma_Lphi[1:i],function(m){m[1,1]})
      ,sigma_Lphi2=sapply(par[,5]$sigma_Lphi[1:i],function(m){m[2,2]})
      ,Cov_Lphi12=sapply(par[,5]$sigma_Lphi[1:i],function(m){m[1,2]})
    )

    if(show_plot==TRUE){
      layout(matrix(c(1:8),4,2,T))
       traceplot(coda::mcmc(tracedata),col="darkred")
       plot(1:5,1:5,t="n",axes=F,ylab="",xlab="")
       text(1.2,3.5,adj = 0
         ,paste("v[",k,"]","=",round(v[,k,T],3),collapse = "\n")
         ,cex=1.2,col=ifelse(accept.status==T,"darkgreen","red"))

    }
    
    i=i+1

  }else{
    if(k=="sigma_Lphi"){
      cat("\n sigma_Lphi is NAN"
          ,"\n","---------------------------\n")
    }else{
      cat("\n r_m is",r_m,",is Frozen at iter",i,",try for different proposal"
          ,"\n","var[",k,"]=",matcor_kk
          ,"\n","v=",v[,k,T]
          ,"\n","----------------------------\n")
    }

  }#end nan

}

## *Post MCMC================================
# last i = 8141
dir.create(as.character(i))

tempy=na.exclude(Y_t[1:i])
save(tempy,file = paste(i,"//",i,"_Y_t.Rdata",sep=""))
templ1=na.exclude(Lphi1)
save(templ1,file=paste(i,"//",i,"_Lphi1.Rdata",sep=""))
templ2=na.exclude(Lphi2)
save(templ2,file=paste(i,"//",i,"_Lphi2.Rdata",sep=""))
tempp=na.exclude(par[1:i,])
save(tempp,file=paste(i,"//",i,"_par.Rdata",sep=""))
templam=na.exclude(lambda)
save(templam,file=paste(i,"//",i,"_lambda.Rdata",sep=""))
tempmu=na.exclude(mu)
save(tempmu,file=paste(i,"//",i,"_mu.Rdata",sep=""))

#save(dev,file=paste(i,"_dev.Rdata",sep=""))
save(matcor,file=paste(i,"//",i,"_matcor.Rdata",sep=""))
save(matcor_kk,file=paste(i,"//",i,"_matcor_kk.Rdata",sep=""))

save(k,file=paste(i,"//",i,"_k.Rdata",sep=""))

#=== *Load MCMC===================================================
{
  i=13511#8141
  # load(file = paste(i,"//",i,"_Y_t.Rdata",sep=""))
  # Y_t<-tempy
  load(file=paste(i,"//",i,"_Lphi1.Rdata",sep=""))
  Lphi1<-LPHI1<-templ1
  load(file=paste(i,"//",i,"_Lphi2.Rdata",sep=""))
  Lphi2<-LPHI2<-templ2
  load(file=paste(i,"//",i,"_par.Rdata",sep=""))
  par=tempp
  load(file=paste(i,"//",i,"_lambda.Rdata",sep=""))
  lambda=templam
  load(file=paste(i,"//",i,"_mu.Rdata",sep=""))
  mu=tempmu
  #save(dev,file=paste(i,"_dev.Rdata",sep=""))
  load(file=paste(i,"//",i,"_matcor.Rdata",sep=""))
  load(file=paste(i,"//",i,"_matcor_kk.Rdata",sep=""))
  load(file=paste(i,"//",i,"_k.Rdata",sep=""))

  alphi=log(c(.2,.5))
  blphi=log(c(4,2))

  pr.lower=c(Lbeta0=log(12),Lbeta_h=log(.8),Leta0=log(1*1/.8*1/.5),Leta_h=log(.8))
  pr.upper=c(Lbeta0=log(12+4),Lbeta_h=log(1.2),Leta0=log(1*1/.8*1/.5+5),Leta_h=log(1.2))
}

##==== *Gather All posterior paramters (fixed and random)===========
# Lphi1<-na.exclude(LPHI1)
# Lphi2<-na.exclude(LPHI2)
{
  colnames(Lphi1)<-paste("ran1",names(X_rep),sep="_")
  colnames(Lphi2)<-paste("ran2",names(X_rep),sep="_")
  burnin0=2500;thin0=2
  lagpar=anylapply(par[1:i,], #par
  function(sub,burnin=burnin0,thin=thin0){
    sub[seq(burnin,length(sub),by=thin)]
    }
    )
  sigma_Lphi_vcov=do.call(rbind,
                        lapply(lagpar$sigma_Lphi,as.numeric))[,-2]
  colnames(sigma_Lphi_vcov)<-c("V_Lphi1","Cov_Lphi12","V_Lphi2")
  par1=cbind(as.data.frame(lagpar[-5]),sigma_Lphi_vcov)
  colnames(par1)<-paste("f",colnames(par1),sep="_")

  lag_Lphi12=cbind(Lphi1[seq(burnin0,i,by=thin0),]
                 ,Lphi2[seq(burnin0,i,by=thin0),])

  par1=cbind(par1,lag_Lphi12)#cbind(par1[seq(burnin0,i,by=thin0),],lag_Lphi12)
  colnames(par1)
  par1=na.exclude(par1)
  par1=par1%>%mutate(f_V_Lphi1=replace(f_V_Lphi1,f_V_Lphi1>=3,NA))%>%drop_na()
}


## *MCMC Diag===================================================================
library(ggmcmc);library(coda)
par1.ggs=ggmcmc::ggs(mcmc(par1))
ggs_autocorrelation(par1.ggs,nLags = 100,family = "f_")
ggs_traceplot(par1.ggs,family = "f_")+facet_wrap(~Parameter,nrow = 4,ncol = 2,scales = "free_y")
ggs_density(par1.ggs,family = "f_")+facet_wrap(~Parameter,nrow = 4,ncol = 2,scales = "free")

ggs_caterpillar(par1.ggs
                 ,thick_ci =c(.1/2,1-.1/2)
                ,line = 0,family = "f_L" )
(HPD=ggmcmc::ci(par1.ggs,thick_ci = c(.1/2,1-.1/2)))

ggs_geweke(par1.ggs,family = "f_")
# ggs_ppsd()
ggs_caterpillar(par1.ggs
                ,thick_ci =c(.25,.75)
                ,thin_ci = c(.25,.75)
                ,line = 0,family = "ran1" )+
                xlim(log(c(.2,4)))

library(bayestestR)
mapest=map_estimate(par1)

##====== set boundary for theta ===========================
# alphi=log(c(.5,.5))
# blphi=log(c(1.7,1.7))
#
# pr.lower=c(Lbeta0=log(4+2),Lbeta_h=log(.8),Leta0=log(1*1/.8*1/.5+2),Leta_h=log(.8))
# pr.upper=c(Lbeta0=log(4+4+1),Lbeta_h=log(1.2),Leta0=log(1*1/.8*1/.5+4+2),Leta_h=log(1.2))

prior.lower=pr.lower
prior.upper=pr.upper
alphi;blphi

prior.l=c(prior.lower,V_Lphi1=0
          ,Cov_Lphi12=min(par1$f_Cov_Lphi12)
          ,V_Lphi2=0
          ,rep(alphi,each=length(X_rep))
)
prior.u=c(prior.upper,V_Lphi1=max(par1$f_V_Lphi1,na.rm=T)
          ,Cov_Lphi12=max(par1$f_Cov_Lphi12)
          ,V_Lphi2=max(par1$f_V_Lphi2)
          ,rep(blphi,each=length(X_rep))
)
names(prior.l)<-names(prior.u)<-colnames(par1)
##===== Define function of theta: =====================
##===== Exp(theta );Cif; Productivity;Morphology=======
prior.l.exp=exp(prior.l)
prior.u.exp=exp(prior.u)

round(cbind(prior.l.exp[1:7],prior.u.exp[1:7]),3)

par1.exp=apply(par1,2,exp)
par1.exp.ggs=ggmcmc::ggs(mcmc(par1.exp))
ggs_caterpillar(par1.exp.ggs
                ,thick_ci =c(.05,.95)
                ,line = 1,family = "f_L" )
ggmcmc::ci(par1.exp.ggs
           ,thick_ci =c(.05,.95))

ggs_density(par1.exp.ggs,family = "f_L")+
  facet_wrap(~Parameter,ncol=3,scales = "free")
ggs_density(par1.exp.ggs,family = "ran1")+
  facet_wrap(~Parameter,ncol=3)


##======== Bayesian Estimator under=====
##scale mean Loss and MAP

A=prior.l.exp;B=prior.u.exp
par1.fun=par1.exp
par1.fun0=par1.fun
library(tidyverse)

par1.fun<-par1.fun%>%
  as.data.frame()%>%
  select("f_Leta0",contains("ran2_"))

par1.fun.t=apply(par1.fun[,-1],2,function(co){
  par1.fun[,1]*co
})


E_par2=apply(par1.fun,2,function(p){ mean(p^2,na.rm=T) })
E_par=apply(par1.fun,2,mean,na.rm=T)

d_iq=numeric(ncol(par1.fun));names(d_iq)=colnames(par1.fun)
for( i_d in 1:length(d_iq)){

  p1num=(B[i_d]*A[i_d])-E_par2[i_d]
  p2num=(E_par2[i_d]-B[i_d]*A[i_d])^2
  p3num=(B[i_d]+A[i_d]-2*E_par[i_d])
  p4num=(2*B[i_d]*A[i_d]*E_par[i_d])-(B[i_d]+A[i_d])*E_par2[i_d]
  dnum=(B[i_d]+A[i_d]-2*E_par[i_d])

  d_iq[i_d]=( p1num+sqrt(p2num-(p3num*p4num)) )/ dnum

  #d_iq[i_d]=sqrt(E_par2[i_d])
}

#par1.fun.ggs=ggmcmc::ggs(mcmc(par1.fun))
#d_iq_HPD=ggmcmc::ci(par1.fun.ggs)

CredInt=bayestestR::ci(as.data.frame(par1.fun)
                       ,ci=.89,method="HDI")
Point_Est=point_estimate(as.data.frame(par1.fun)
                         ,centrality = "all"
                         ,dispersion = T)
Point_CredInt=Point_Est%>%
  #arrange(order(names(d_iq)))%>%
  tibble::add_column(CredInt[,-1]
                     ,.before="Median"
                     ,d_iq=d_iq
  )

ggplot(Point_CredInt[grep("f_",Point_CredInt$Parameter),]
       ,aes(y=reorder(Parameter,d_iq),x=d_iq))+
  geom_linerange(mapping = aes(xmin=CI_low,xmax=CI_high),size=.5)+
  geom_point(size=2,col="red")+
  geom_vline( xintercept=1,lwd=.5,lty=2)+
  scale_y_discrete(name="Parameters")
  #geom_point(aes(x=Median),size=2)+
  #geom_point(aes(x=MAP),col="royalblue",size=2)
#+geom_linerange(mapping = aes(xmin=Low,xmax=High),size=1.3)
ggsave("d_iq.jpeg",device="jpeg",dpi=400)

Point_CredInt[grep("f_",Point_CredInt$Parameter),]%>%write.csv(file="d_iq.csv")

Point_CredInt[grep("ran1",Point_CredInt$Parameter),]

ggplot(Point_CredInt[grep("ran2_",Point_CredInt$Parameter),]
       ,aes(y=reorder(Parameter,d_iq) #d_iq)
            ,x=d_iq
            ))+
  # geom_linerange(mapping = aes(xmin=CI_low,xmax=CI_high),size=1)+
  # geom_vline(xintercept = 1,lty=2,col="blue")+
  geom_point(size=2,col="red")

  #+xlim(c(.5,2)
# +geom_point(aes(x=Median),size=2)+
# geom_point(aes(x=MAP),col="royalblue",size=2)
 
## random effects density plot=============================
ggplot(Point_CredInt[grep("ran2_",Point_CredInt$Parameter),]
        ,aes(x=d_iq))+
  geom_density(aes(fill="red"),alpha=.4)+
  geom_density(
    data = Point_CredInt[grep("ran1_",Point_CredInt$Parameter),],
    aes(x=d_iq,fill="blue"),alpha=.4
  )+
  xlab(expression(d_iq))+
  scale_fill_discrete(name="Random effect"
                      ,labels=c(expression(phi[2])
                                ,expression(phi[1]))
                      )+
  theme_light(base_size = 16)+
  theme(legend.position = "top")


#======* GOF =============================
##------------- Predictive Residual----------------------
dim_aal_img=dim(aal_image())
dim_aal_img
labs = aal_get_labels()


### N_obs================================================
N_obs_s<-lapply(X_rep,function(X_rep0){
  c_obs=unmark(X_rep0)
  c_obs=as.matrix(c_obs)
  xyz.c_obs=sweep(c_obs, 2,dim_aal_img , "*")
  lab.c_obs=aal_lookup(xyz.c_obs)
  newNobs_k=table(lab.c_obs$name)

  Nobs_k=rep(0,length(labs$name))
  names(Nobs_k)=labs$name
  Nobs_k[names(Nobs_k)%in%names(newNobs_k)]=newNobs_k
  return(Nobs_k)
})

### N_pred ===============================================
K=length(labs$name)
S=length(X_rep)
burnin_pred=2000#7000
last_pred=2100#length(Y_t)
newY_t=Y_t[seq(burnin0,last_pred,by=thin0)]
(L=length(newY_t))

Npred_list<-vector("list",length = length(X_rep))
names(Npred_list)<-names(X_rep)

library(doFuture)
future:::ClusterRegistry("stop")
plan(muticore,workers=15)
#plan(multisession,workers=3)
registerDoFuture()
handlers(handler_progress())
with_progress({
  pb_N_pred=progressor(along = 1:S)

  Npred_list=foreach(s = 1:S)%dopar%{
    pb_N_pred(sprintf("N_pred_k = %g",s))

    foreach(l = 1:L,.export ="Npred_list" )%dopar%{

      c_pred=unmark(newY_t[[l]][[s]])
      c_pred=as.matrix(c_pred)
      xyz.c_pred=sweep(c_pred, 2,dim_aal_img , "*")
      lab.c_pred=aal_lookup(xyz.c_pred)
      newNpred_k=table(lab.c_pred$name)

      Npred_k=rep(0,length(labs$name))
      names(Npred_k)=labs$name
      Npred_k[names(Npred_k)%in%names(newNpred_k)]=newNpred_k

      return(Npred_k)
      #Npred_list[[s]][[l]]<-Npred_k

    } #end foreach l
  }##end foreach s

}) #end progress
Npred_list[[25]][50]

#============= R_pred===============================
R_Pred_list<-vector("list",length = length(X_rep))
names(R_Pred_list)<-names(X_rep)
R_Pred_K_list<-vector("list",length = K)

for(s in 1:S){
  for(k in 1:K){
    for(l in 1:L){
      R_Pred_K_list[[k]][l]<- N_obs_s[[s]][k]-Npred_list[[s]][[l]][k]

    }
  }
}
for(s in 1:S){
  for(k in 1:K){
    R_Pred_list[[s]][[k]]<-R_Pred_K_list[[k]]
  }
}

R_Pred_list[[S]]
##============ Coverage Provability ==========================
R_Pred_CI<-vector("list",length = length(X_rep))
for(s in 1:S){
  for(k in 1:K){

    R_Pred_CI[[s]][[k]]=quantile(R_Pred_list[[s]][[k]]
                                 ,probs = c(.05,1-.05))
  }
}

CP=numeric(length(X_rep));names(CP)=names(X_rep)
for(s in 1:S){
  for(k in 1:K){
    CP[s]=mean(sign(R_Pred_CI[[s]][[k]]) <=0)

  }
}

as_tibble(CP)%>%filter(value!=1)

##========== show N_pred and N_obs=========================
temp_aal=aal_image()
contour3d(temp_aal,level = 4500,alpha = .3)
k=sig_labels
#===== *Classic GOF================================
study=10
subwin=bounding.box3(X_rep[[study]])
par(mfrow=c(1,1))
plot(X_rep[[study]])
npoints(X_rep[[study]])
burnin0=1;thin0=1
tempyt=Y_t[seq(burnin0,i,by=thin0)]
sim_n_combine=lapply(tempyt,function(y_t){
  sim_n_i=unlist(sapply(y_t,npoints.pp3))
  sim_n=sum(sim_n_i)

  sim_n_subwin=npoints.pp3(y_t[[study]][subwin])

  return(c(sim_n=sim_n,sim_n_subwin=sim_n_subwin))
})

sim_n_combine=do.call(rbind,sim_n_combine)

sim_n_combine1=ggmcmc::ggs(mcmc(sim_n_combine))
ggs_autocorrelation(sim_n_combine1,nLags = 100)
ggs_traceplot(sim_n_combine1)
ggs_histogram(sim_n_combine1)

ggs_caterpillar(sim_n_combine1,thick_ci =c(.05,.95) )+
  facet_wrap(~Parameter,nrow = 2,scales = "free")
ggmcmc::ci(sim_n_combine1,thick_ci = c(.05,.95))

ggs_geweke(sim_n_combine1)



#===Run Cif via Batch ====================================
# system("Rep_Cif/Rep_cif.sh")
# load("Rep_Cif/REP_CIF_real_64_sim1000.RData")
load(paste("Rep_Cif/REP_CIF_real_91_sim"
           ,5559#nrow(par1)
           ,".RData",sep=""))
#===*Bayesian Cif=========================================

cc=vector("list",length(X_rep));names(cc)=names(X_rep)

for(i in 1:length(REP_CIF)){
  for(j in 1:length(REP_CIF[[i]])){
    cc[[j]][[i]]=REP_CIF[[i]][[j]]
  }

}

all_study_cif=lapply(cc,function(study){
  do.call(rbind,lapply(study,marks))
})

minCIF=min(sapply(all_study_cif,min,simplify=T))
maxCIF=max(sapply(all_study_cif,max,simplify=T))
all_study_Bayes_cif=lapply(all_study_cif,function(onecif){
  #point_estimate(exp(as.data.frame(onecif))
  #                ,centrality="map")[,2]

  #bayestestR::ci(as.data.frame(onecif)
  #                 ,ci=.89,method="HDI")[,3]

  #  apply(onecif,2,function(p){
  #  sqrt(mean(exp(p)^2))
  #  })

  onecif.fun=  exp(onecif)
  A=apply(onecif.fun,2,function(x) 0 ) #min)
  B=apply(onecif.fun,2,function(x) exp(maxCIF))#max)
  E_par2=apply(onecif.fun,2,function(p){ mean(p^2) })
  E_par=apply(onecif.fun,2,mean)

  d_iq=numeric(ncol(onecif.fun));names(d_iq)=colnames(onecif.fun)
  for( j1 in 1:length(d_iq)){

    p1num=(B[j1]*A[j1])-E_par2[j1]
    p2num=(E_par2[j1]-B[j1]*A[j1])^2
    p3num=(B[j1]+A[j1]-2*E_par[j1])
    p4num=(2*B[j1]*A[j1]*E_par[j1])-(B[j1]+A[j1])*E_par2[j1]
    dnum=(B[j1]+A[j1]-2*E_par[j1])

    d_iq[j1]=( p1num+sqrt(p2num-(p3num*p4num)) )/ dnum

  }

  d_iq
})

REP_CIF_pp3=lapply(REP_CIF[[1]],unmark)
for(j2 in 1:length(X_rep)){
  marks(REP_CIF_pp3[[j2]])<-all_study_Bayes_cif[[j2]]
}
REP_CIF_pp3[1]

#====== All study Bayes Cif Credible============
all_study_Bayes_cif_cred=lapply(all_study_cif,function(onecif){

  onecif.fun=  onecif
  l.cred.cif=bayestestR::ci(as.data.frame(onecif.fun)
                  ,ci=.95,method="HDI")[,3] # HDI
  u.cred.cif=bayestestR::ci(as.data.frame(onecif.fun)
                  ,ci=.95,method="HDI")[,4]
  cred.cif.sig=l.cred.cif>=0
  
  return( cred.cif.sig)
})

#====== Ready for Cif plot =====================
r=.05
template <- sch_mni152#readNIfTI(system.file("MNI152_T1_2mm_brain.nii.gz"
                       #           , package="brainR")
                      #, reorient=FALSE)
(dtemp=dim(template))
# ------------------scaled r to voxel --------------------
(radii_v=r*dtemp[1])  #--> should be > 1(voxel)
(radii_1v=1)
# ------------------voxel to mm (if each voxel be 2 mm)----
template@pixdim[2]
(radii_m<-radii_v*template@pixdim[2])

#-----colour -----------
X_rep=joint_SD_h$pts3_norm
library(viridis)
rep_alpha<-rep_col<-vector("list",length(X_rep))
names(rep_col)<-names(REP_CIF_pp3)

combine_rep_cif=unlist(lapply(REP_CIF_pp3,marks))
hist(combine_rep_cif)
minc=min(combine_rep_cif)
maxc=max(combine_rep_cif)

rep_col=lapply(REP_CIF_pp3,function(cif){
  tblcif=as_tibble(marks(cif))%>%count(value)
  mincol=(min(marks(cif))-minc)/(maxc-minc)
  maxcol=(max(marks(cif))-minc)/(maxc-minc)
  hc=plasma(nrow(tblcif)
            ,alpha=seq(mincol,maxcol,len=nrow(tblcif))
            ,begin = mincol
            ,end = maxcol
            ,direction =1)
  hc=rep(hc,tblcif$n)
  
  hc[rank(marks(cif))]
})


#= Plot CIF============================
#= Combined PLot=========================
{
  #open3d()
  options(rgl.printRglwidget = FALSE)
  dtemp=dim(template)

  par3d(mouseMode="trackball",windowRect=c(  67,83,1360/2,796))
  rgl.bg(color=c("black","black"))
  rgl.viewpoint( theta = 0, phi = -215, fov = 60, zoom = 1,
                 scale = par3d("scale"), interactive = TRUE,
                 type = "modelviewpoint")
  contour3d(template, level = 1#4500
            ,smooth = "grid",material = "shine"
            ,alpha = .3, draw = TRUE,add=T)

  sig_coords_list=vector("list",length = length(REP_CIF_pp3))
  names(sig_coords_list)=names(REP_CIF_pp3)
  
  for( p.index in 1:length(REP_CIF_pp3)){
    cat("plot study",p.index,"\n")

    cif_voxel0=coords(REP_CIF_pp3[[p.index]])
    cif_voxel_mark=marks(REP_CIF_pp3[[p.index]])
    qCIF89=quantile(combine_rep_cif,probs = .95) #quantile(cif_voxel_mark,probs = .89)
    cif_voxel=cif_voxel0[cif_voxel_mark>=qCIF89,]
    
    new_cif_voxel=cbind.data.frame(
      cif_voxel=cif_voxel0[cif_voxel_mark>=qCIF89,],
      cif=cif_voxel_mark[cif_voxel_mark>=qCIF89]
    )
    #cif_voxel=cif_voxel0[all_study_Bayes_cif_cred[[p.index]],]
    #-------- aal tag--------------------------------
    sig_coords_list[[p.index]]<-new_cif_voxel#cif_voxel
    #---------end aal tag -------------------------------
    if(nrow(cif_voxel)>0){
      rgl.spheres(
        x =cif_voxel[,1]*dtemp[1]
        ,y =cif_voxel[,2]*dtemp[2]
        ,z =cif_voxel[,3]*dtemp[3]
        ,r=radii_v #1
        ,col=rep_col[[p.index]] [cif_voxel_mark>=qCIF89]
        #,col=rep_col[[p.index]] [all_study_Bayes_cif_cred[[p.index]]]
        #,alpha=.4
        ,add=T
      )

      # rgl.points(
      #   x =unlist(X_rep[[p.index]]$data$x)*dtemp[1]
      #   ,y =unlist(X_rep[[p.index]]$data$y)*dtemp[2]
      #   ,z =unlist(X_rep[[p.index]]$data$z)*dtemp[3]
      #   #,r=1
      #   ,col="black"
      #   ,add=T
      # )
    }

  } # end for

  bgplot3d(fields::image.plot(legend.only = T
                              ,horizontal = T
                              , zlim=c(log(minc)
                                        ,log(maxc))
                              , col =plasma(begin = 0
                                            ,n=length(unlist(rep_col))
                                            ,direction =1)
                              , legend.args = list(text="Log (CIF)",col="white",cex=1.5)
                              , axis.args=list(col.axis="white",col.ticks="white",cex.axis=1.5)
                              , add=T)
           ,bg.color =c("black","black"))
}

text3d(x=dtemp[1]/2, y=dtemp[2]/2, z = dtemp[3]*0.98, text="Top",col="white")
text3d(x=-0.98, y=dtemp[2]/2, z = dtemp[3]/2, text="Right",col="white")

sig_coords_df0=do.call(rbind,sig_coords_list)
sig_coords_df0=rownames_to_column(sig_coords_df0)
write.csv(sig_coords_df0,"sig_coords_df.csv",row.names = F)

movie3d(
  movie="3d_filterd_Max_cif_map",
  spin3d( axis = c(0, 0, 1), rpm = 10),
  duration = 8,
  dir = getwd(),
  type = "gif",
  clean = TRUE
)

#writeWebGL(dir="webGL",file.path("webGL", "index.html"))
#htmlwidgets::saveWidget(rglwidget(), "3d_filterd_cif_map.html")
#rgl.bbox(color = c("#333377", "white"), emission = "#333377",
#         specular = "#3333FF", shininess = 5, alpha = 0.8 ,font.size=.5)

#==== Seperate plot========================
{
  options(rgl.printRglwidget = FALSE)
  dtemp=dim(template)
  mfrow3d(byrow = T,2,3,sharedMouse = T)

  foreach( p.index=1:length(REP_CIF_pp3))%do%{
    if(p.index%%7==0){ open3d();mfrow3d(byrow = T,2,3,sharedMouse = T)}

    cat("plot study",p.index,"\n")
    next3d()
    rgl.bg(color=c("black","black"))
    par3d(mouseMode="trackball",windowRect=c(  67,83,1360,796))

    rgl.viewpoint( theta = 0, phi = 15, fov = 60, zoom = 1,
                   scale = par3d("scale"), interactive = TRUE,
                   type = "modelviewpoint")

    # pch3d( x = X_rep[[p.index]]$data$x*dtemp[1]
    #        ,y =X_rep[[p.index]]$data$y*dtemp[2]
    #        ,z =X_rep[[p.index]]$data$z*dtemp[3]
    #        ,radius =2
    #        ,pch = 16
    #        # ,cex=.2
    #        # ,alpha=alpha_rep
    #        ,col=col_cif[[p.index]])
    cif_voxel0=coords(REP_CIF_pp3[[p.index]])
    cif_voxel_mark=marks(REP_CIF_pp3[[p.index]])

    # cif_voxel=cif_voxel0[cif_voxel_mark>0,] 
     cif_voxel=cif_voxel0[all_study_Bayes_cif_cred[[p.index]],]

    if(nrow(cif_voxel)>0){
      rgl.spheres(
        x =cif_voxel[,1]*dtemp[1]
        ,y =cif_voxel[,2]*dtemp[2]
        ,z =cif_voxel[,3]*dtemp[3]
        ,r=radii_v# 1
        # ,col=rep_col[[p.index]][cif_voxel_mark>0]
        ,col=rep_col[[p.index]] [all_study_Bayes_cif_cred[[p.index]]]
        #,alpha=.4
        ,add=T
      )
    }
  
  bgplot3d(fields::image.plot(legend.only = T
                              ,horizontal = T
                              , zlim=c(minc,maxc)
                              , col =plasma(begin = 0
                                            ,n=length(unlist(rep_col))
                                            ,direction =1)
                              , legend.args = list(text=paste("Log (CIF)",names(REP_CIF_pp3[p.index]))
                                                   ,col="white",cex=1)
                              , axis.args=list(col.axis="white",col.ticks="white",cex.axis=1.5)
                              , add=T)
           ,bg.color =c("black","black"))

    # bgplot3d(fields::image.plot(legend.only = T,horizontal = T
    #                             , zlim =range(cif_voxel_mark)+c(-.01,.01)
    #                             , col =plasma(begin = 0
    #                                           ,length(cif_voxel_mark)
    #                                           ,direction =1)
    #                             , legend.args = list(text=paste("CIF\n",names(X_rep)[p.index],sep="")
    #                                                  ,col="white")
    #                             , axis.args=list(col.axis="white",col.ticks="white")
    #                             , add=T)
    #          ,bg.color =c("black","white"))
    contour3d(template
              ,level = 4500
              ,smooth = "grid",material = "shine"
              ,alpha = .2, draw = TRUE,add=T)


    # rgl.bbox(color = c("#333377", "white"), emission = "#333377",
    #          specular = "#3333FF", shininess = 5, alpha = 0.8 ,font.size=.5)


  }
  highlevel(integer()) #
}


#===== Tag brain regions corresponding to indices =============
library(neurobase)
library(aal)
library(MNITemplate)

sig_coords_df=sig_coords_df[,1:3]
plot3d(sig_coords_df)
##----- Find significant lables---------------------------------
study_specific_sig=lapply(sig_coords_list,function(sig_coords){
  if(nrow(sig_coords)>0){
    xyz.foci=sweep(sig_coords, 2,dim(template), "*")
    return(aal_lookup(as.matrix(xyz.foci)))
   
     }
})

#= find significant cif on brain slices--------------------------
save(sig_coords_list,file = "sig_coords_list.Rdata")
for(ss in 10:75 ){
  jpeg(paste("sig_lab_slice//",ss,".jpeg",sep=""))
  oro.nifti::slice(template,z=ss)
  text(0,.95,adj=0,label="R",col="white",cex=2.2)
  lapply(sig_coords_list,function(sig_coords){
  if(nrow(sig_coords)>0){
    xyz.foci=sweep(sig_coords[,1:3], 2, dim(template), "*")
  
    #-----------------
    xyz.foci_zslice=cbind(sig_coords[,1:3],slice_z=ceiling(xyz.foci[,3]))
    o=xyz.foci_zslice%>%
      filter(slice_z==ss)%>%select(-3:-4)
    points(o[,1],o[,2],pch=16,col=2,cex=0.05*91/2)# 2 for 2mm
    ss0=world.to.mni(cbind(0,0,ss),nii.file="MNI152_T1_2mm_brain.nii")
    title(main=paste("z=",ss0[3]),col.main="white")
  }
    
  })
  dev.off()
}
system("convert -delay 30 -loop 0 sig_lab_slice/*.jpeg sig_ab_slice_sch.gif")

#= sginificant table area=======================
study_specific_sig_df=do.call(rbind,study_specific_sig)
sig_labels=as.data.frame(table(study_specific_sig_df$name))
sig_labels
ggplot(sig_labels,aes(x=reorder(Var1,Freq),y=Freq))+
  geom_bar(stat="identity",)+coord_flip()+
  theme_light()
ggsave(file="sig_labels.jpeg",device="jpeg",dpi=500)
#- save significant lables-------------------------
write.csv(sig_labels,"sig_labels.csv")

#sig_labels=sig_labels%>%dplyr::filter(Freq>=1)
#== Highlight sig labels on Brain ===================
img = aal_image()
cut <- 4500
dtemp <- dim(template)
# All of the sections you can label
labs = aal_get_labels()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  #open3d()
  options(rgl.printRglwidget = FALSE)
  dtemp=dim(template)
  par3d(mouseMode="trackball",windowRect=c(  67,83,1360/2,796))
  rgl.bg(color=c("black","black"))
  rgl.viewpoint( theta = 145, phi = 16, fov = 60, zoom = 1,
                 scale = par3d("scale"), interactive = TRUE,
                 type = "modelviewpoint")
  contour3d(template, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3]
            , level = cut
            , alpha = 0.4, draw = TRUE)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sig_lab_col_code=punif(sig_labels$Freq,min(sig_labels$Freq),max(sig_labels$Freq))  

  sig_lab_col= plasma(length(table(sig_lab_col_code)))
  sig_lab_col=rep(sig_lab_col,table(sig_lab_col_code))
  sig_lab_col=sig_lab_col[rank(sig_lab_col_code)]
  #rainbow(nrow(sig_labels) )

  for (i in 1:nrow(sig_labels)){
    sig_lab_index=labs$index[labs$name%in%sig_labels$Var1[i]]
    # mask significant Areas
    mask_sig=remake_img(vec = img %in% sig_lab_index, img = img)
    #mask_sig_ind=which(mask_sig==1,arr.ind = T)
    #dim(mask_sig_ind)
    # this would be the ``activation'' or surface you want to render -----
   
    contour3d(mask_sig, level = c(0.3)
              , alpha = c(0.4)
              , add = TRUE,smooth = 80
              , color=sig_lab_col[i]
              )
     mask_ave_loc=apply(which(mask_sig@.Data==1,arr.ind=T),2,mean)
     text3d(x=-3, y=dtemp[2]/2+1, z = dtemp[3]/2+1
            , text=as.character(sig_labels$Var1[i]),col="red")
     arrow3d(
        c(x=-3, y=dtemp[2]/2, z = dtemp[3]/2)
       ,c(mask_ave_loc[1],mask_ave_loc[2],mask_ave_loc[3])
       ,col=2)
    cat(as.character(sig_labels$Var1[i]),"\n")

   }
}
   bgplot3d(fields::image.plot(legend.only = T
                              ,horizontal = T
                              , zlim=range(sig_labels$Freq)
                              , col =plasma(length(sig_lab_col_code)
                                              ,begin=min(sig_lab_col_code)
                                              ,end=max(sig_lab_col_code)
                                              # ,alpha=sort(sig_lab_col_code)
                                            )
                              , legend.args = list(text="#points per region",col="white",cex=1.5)
                              , axis.args=list(col.axis="white",col.ticks="white",cex.axis=1.5)
                              , add=T)
           ,bg.color =c("black","black"))
}

#------------- add text------------------------------------------------
text3d(x=dtemp[1]/2, y=dtemp[2]/2, z = dtemp[3]*0.98, text="Top",col="white")
text3d(x=-0.98, y=dtemp[2]/2, z = dtemp[3]/2, text="Right",col="white")

#writeWebGL(filename = "activation_regions.html")

#rglwidget()
movie3d(
  movie="activation_regions",
  spin3d( axis = c(0, 0, 1), rpm = 10),
  duration = 8,
  dir = getwd(),
  type = "gif",
  clean = TRUE
)

#=========Visualizing arbitrary vertices==================

library(ggseg)
ggseg(colour="black", size=.7
      ,adapt_scales = FALSE
      ,position = "stacked"
      , mapping=aes(fill=region)) +
  theme_void()

someData =data.frame(area =sig_labels$Var1,p =1:nrow(sig_labels))# Figure 4
ggseg(.data=someData, mapping=aes(fill=p))+
  labs(title="A nice plot title", fill="p-value")+
  scalefillgradient(low="firebrick",high="goldenrod")

#=======================================================
img = aal_image()

#template = readMNI(res = "2mm")
cut <- 4500
dtemp <- dim(template)
# All of the sections you can label
labs = aal_get_labels()
# Pick the region of the brain you would like to highlight - in this case the hippocamus_L
hippocampus = labs$index[grep("Parietal_S", labs$name)]
mask = remake_img(vec = img %in% hippocampus, img = img)
mask.ind=which(mask==1,arr.ind = T)
dim(mask.ind)

### this would be the ``activation'' or surface you want to render
contour3d(template, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], level = cut, alpha = 0.1, draw = TRUE)
contour3d(mask, level = c(0.5), alpha = c(0.5), add = TRUE, color=c("red") )
### add text
text3d(x=dtemp[1]/2, y=dtemp[2]/2, z = dtemp[3]*0.98, text="Top")
text3d(x=-0.98, y=dtemp[2]/2, z = dtemp[3]/2, text="Right")
rglwidget()


## END. ===================================================================

###########################################################################
mod01 <- list(cif="strauss",par=list(beta=400,gamma=0.2,r=.05),
              w=c(0,1,0,1))
X1.strauss <- rmh(model=mod01,start=list(n.start=ns),
                  control=list(nrep=nr,nverb=nv))
plot(X1.strauss)
# X1.strauss=as.data.frame(X1.strauss)



mod02 <- list(cif="",par=list(beta=400,eta=1.5,r=.05),
              w=c(0,1,0,1))
X1.geyer <- ragsAreaInter(200,2,.05,win=owin())
# plot(X1.geyer)
X1.geyer=as.data.frame(X1.geyer)
plot(
  table(cutree(hclust(dist(X1.strauss)),h = .05))
)
plot(
  table(cutree(hclust(dist(X1.geyer)),h = .05))
)
plot(
  table(cutree(hclust(dist(im.foci)),h = .05))
)

sum(table(cutree(hclust(dist(im.foci)),h = .05))>1)
sum(table(cutree(hclust(dist(X1.strauss)),h = .05))>1)
sum(table(cutree(hclust(dist(X1.geyer)),h = .05))>1)

plot(hclust(dist(X1.strauss)),hang=-1)
plot(hclust(dist(X1.geyer)),hang=-1)
plot(hclust(dist(im.foci)),hang=-1)








