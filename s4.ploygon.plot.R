
args <- commandArgs(TRUE)
polygon<-args[1]
prefix <-args[2]
output_dir<-args[3]
source('/zfssz3/BC_COM_P7/F17FTSECWLJ2619/MANiggR-01/scRNA/Spatial/ovary/bin/ST.subfunction.xiu.R')
source('/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/yanghuanjie/project/FM/01.ST/gem.to.rds.r')


all <- read.table(polygon,header=T)
complete_data<-NULL
  for (i in unique(all$CellID))
    {
     single=all[all$CellID==i, c("x","y", "CellID","OO.genes","GC.genes","GC.C.genes","GC.M.genes","GC.3.genes","TC.genes","smc.genes","endo.genes","Macrophage.genes"),]
     single.chull<- single[chull(single[single$CellID== i, c("x", "y")]),]
     complete_data<-rbind(complete_data,data.frame(single.chull)) 
    }
