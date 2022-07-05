library(future.apply)
# library(alphashape3d)
library(spatstat)
library(progressr)

load("joint_SD_h.RData")

# load("hh_sim.RData")
future:::ClusterRegistry("stop")
plan(multicore,workers=10)
# plan(multisession,workers=4)

# sphere_volume_mc<-function(n=2e6
#                            ,brain=im.foci.ashape
#                            ,S.centers=im.foci3,S.r=.05,
#                            win=box3(),...){
#   cube=runifpoint3(n,domain = win,...)
#   #---------- start parallel inashape3d------------------------------------
#   NPOINT_cube=n
#   CUT=30; if(NPOINT_cube<=300) CUT=3
#
#   in_brain<- do.call(c,future_lapply(split.ppx(cube,cut(1:NPOINT_cube,CUT))
#                                      ,FUN = function(pm){
#                                        in_brain=inashape3d(brain,points = as.matrix.ppx(pm[[1]]))
#
#                                        return(in_brain)
#                                      })# end lapply
#   )# end do.call
#
#   in_union<-has.close.pp3(cube,S.r/2,S.centers)
#
#   #---------- end parallel inashape3d-------------------------------------
#
#   mcvol=mean(in_brain & in_union)*volume(win)
#   return(mcvol)
# }
#
# system.time({
#    Area3d_rep<-future_lapply(hh_sim$Pts3,function(pts3_i){
#      cat("1 end\t")
#      return(
#        sphere_volume_mc(n=6e6,win = box3()
#                         ,S.centers =pts3_i
#                         ,S.r = 0.07)
#      )
#    }
#    )
# })
#
# ## with n=1e7-->
# Area3d_rep_5e6_r07<-Area3d_rep
# save(Area3d_rep_5e6_r07,file = "Area3d_rep_5e6_r0.07.RData")
#

############## Grid based Volume ##########################

#========= Uncomment if you need change the resolution other than 64,128,264========

# allvoxel=264^3#prod(dtemp)
# z<-y<-x<-seq(0,1,len=allvoxel^(1/3))
# (pric=1/allvoxel^(1/3))
# cube=expand.grid(x,y,z)
# CUT=30
# cube_split=split.data.frame(cube,cut(1:nrow(cube),CUT))
#  with_progress({
#    pb_split=progressor(along=1:CUT)
#    cub3_split=future_lapply(cube_split,function(d){
#      pb_split(sprintf("%g",1))
#      pp3(x=d$Var1,y=d$Var2,z=d$Var3
#          ,domain=box3(range(d$Var1),range(d$Var2),range(d$Var2))
#      )
#
#    })
#  })
# save(cub3_split,file="cub3_split_264.RData")

#================for test======================
# p=pp3(.5,.5,.5,domain=box3())
# a=has.close.pp3(cube3,.05,im.foci3)
# sum(a)/npoints.pp3(cube3)
# 4/3*pi*((.05)^3)

#===============run grid MC_vol===============

source("3dpp.R")
load(paste("cub3_split_",resol[1],".RData",sep="")) # resol is set by default via source("3dpp.R")
#load("cub3_split_128.RData")
#load("cub3_split_264.RData")
# system.time(
#   sphere_volume_mc()
# )

cat("start sphere_volume_mc with",resol,"resolution\n")
system.time({
  with_progress({
    pb_study<-progressor(along = 1:length(joint_SD_h$pts3_norm))
    Area3d_rep<-future_lapply(future.seed = TRUE
      ,joint_SD_h$pts3_norm
      ,function(pts3_i){
      pb_study(sprintf("%g",1))
      return(
        sphere_volume_mc(S.r = .05
                         ,S.centers = pts3_i,cube.pp3 = cub3_split)
      )
    }
    )
  })

})


Area3d_rep_grid_0.05_res91<-Area3d_rep
save(Area3d_rep_grid_0.05_res91,file = paste("Area3d_rep_grid_0.05_res",resol[1],".RData",sep=""))
cat("Replicated Area saved in:",paste("Area3d_rep_grid_0.05_res",resol[1],".RData",sep=""),"\n")
