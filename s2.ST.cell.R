
##############################################my file ###################################################


args <- commandArgs(TRUE)
gemfile<-args[1]
prefix<-args[2]
output_dir<-args[3]
source('/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/yanghuanjie/project/ovary/bin/ST.subfunction.cell.R')


gem=gemfile;
seurat_spatialObj=gem_2_rds(gem,prefix,20)
#（0）设置QC阈值：
min_gene <- 1
max_gene <- 20000
min_count <-50
max_count <- 40000
mincell <-5
mt_cut <- 30
cv_method <- 'vst'; # vst, mvp, disp
top_cv_select <- 3000 #高可变基因
use_pcs_max <- 50; #聚类需要用到前几个PCs
resolution <- 0.5; #聚类的分辨率，值越高，类数目越多，0.4-1.2范围


#（1）进行QC质控评估
seurat_spatialObj=calc_percenterage_feature_set(seurat_spatialObj)
options(repr.plot.width = 16,repr.plot.height =4)
pdf(paste0(output_dir,prefix,".raw.VlnPlot.pdf"),w=10,h=6)
VlnPlot(seurat_spatialObj, features = c("nCount_Spatial","nFeature_Spatial","percent.mt","percent.ribo","percent.hb"), 
               cols="blue",pt.size=0,ncol = 5,fill.by="blue") 
dev.off()

##################################################################
st_dat=seurat_spatialObj
#st_dat=filter_st_dat(seurat_spatialObj, minFeature=min_gene, maxFeature=max_gene, minCount=min_count, maxCount=max_count, minCell=mincell, mtRatio=mt_cut)

#barplot(seurat_spatialObj@assays$Spatial@counts,output_dir,paste0(prefix,"raw"))
#barplot(st_dat@assays$Spatial@counts,output_dir,paste0(prefix,"filtered"))


pdf(paste0(output_dir,prefix,".filtered.VlnPlot.pdf"),w=10,h=6)
VlnPlot(st_dat, features = c("nCount_Spatial","nFeature_Spatial","percent.mt","percent.ribo","percent.hb"), 
               cols="blue",pt.size=0,ncol = 5,fill.by="blue")
dev.off() 


dim(seurat_spatialObj) #原始的
dim(st_dat) #过滤后的
seurat_spatialObj@images$image = new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = seurat_spatialObj@meta.data[, c('y', 'x')]
    )
st_dat@images$image = new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = st_dat@meta.data[, c('y', 'x')]
    )
saveRDS(seurat_spatialObj,file=paste0(output_dir,prefix,".raw.rds"))
saveRDS(st_dat,file=paste0(output_dir,prefix,".ST.rds"))

pdf(paste0(output_dir,prefix,".iPlot.Seurat.pdf"),w=8,h=8)
p1=SpatialFeaturePlot(seurat_spatialObj, features="nFeature_Spatial",stroke = NA, alpha = c(0.7, 1),images="image")
p2=SpatialFeaturePlot(seurat_spatialObj, features="nCount_Spatial",stroke = NA, alpha = c(0.7, 1),images="image")
p3=SpatialFeaturePlot(st_dat, features="nFeature_Spatial", stroke = NA, alpha = c(0.7, 1),images="image")
p4=SpatialFeaturePlot(st_dat, features="nCount_Spatial", stroke = NA, alpha = c(0.7, 1),images="image")
multiplot(p1,p2,p3,p4,cols=2)
dev.off()

pdf(paste0(output_dir,prefix,".iPlot.pdf"),w=8,h=8)
p1=iPlot(seurat_spatialObj, features="nFeature_Spatial", pt.size = 0.04, cluster_Palette = ColorPalette(50))
p2=iPlot(seurat_spatialObj, features="nCount_Spatial", pt.size = 0.04, cluster_Palette = ColorPalette(50))
p3=iPlot(st_dat, features="nFeature_Spatial", pt.size = 0.04, cluster_Palette = ColorPalette(50))
p4=iPlot(st_dat, features="nCount_Spatial", pt.size = 0.04, cluster_Palette = ColorPalette(50))
multiplot(p1,p2,p3,p4,cols=2)
dev.off()
#
st_dat =Normalization(st_dat,use_pcs_max)
OO.genes<- c("ZP3","TUBB8","FIGLA","DDX4")
GC.genes <- c("AMH","GSTA1","WT1","INHA")
GC.C.genes <- c("TNNI3","AMH","DSP","FST","HSD17B1","SERPINE2","MAGED2","PRKAR2B")
GC.M.genes <-c("INSL3","APOE","S100A6","GSTA1","APOA1","FDX1","CYP17A1") 
GC.3.genes <- c("DCN","LGALS1","LGALS3","ADIRF","IGFBP5","ADAMTS3","ADAMTS1","ZNF331")
TC.genes <-c("DCN","STAR")
smc.genes <- c("ACTA2","MUSTN1","DES")
endo.genes <- c("TM4SF1","VWF","TREK")
Macrophage.genes <-c("CD68","CD14")

allgene=rownames(st_dat@assays$Spatial)
OO.genes <- intersect(OO.genes,allgene)
GC.genes <- intersect(GC.genes,allgene)
GC.C.genes <- intersect(GC.C.genes,allgene)
GC.M.genes <- intersect(GC.M.genes,allgene)
GC.3.genes <- intersect(GC.3.genes,allgene)
TC.genes <- intersect(TC.genes,allgene)
smc.genes <- intersect(smc.genes,allgene)
endo.genes <- intersect(endo.genes,allgene)
Macrophage.genes <- intersect(Macrophage.genes,allgene)
st_dat <- AddModuleScore(st_dat,features = list(OO.genes),name="OO.genes",n.bin = 500,ctrl.size=500)
st_dat <- AddModuleScore(st_dat,features = list(GC.genes),name="GC.genes",n.bin = 500,ctrl.size=500)
st_dat <- AddModuleScore(st_dat,features = list(GC.C.genes),name="GC.C.genes",n.bin = 500,ctrl.size=500)
st_dat <- AddModuleScore(st_dat,features = list(GC.M.genes),name="GC.M.genes",n.bin = 500,ctrl.size=500)
st_dat <- AddModuleScore(st_dat,features = list(GC.3.genes),name="GC.3.genes",n.bin = 500,ctrl.size=500)
st_dat <- AddModuleScore(st_dat,features = list(TC.genes),name="TC.genes",n.bin = 500,ctrl.size=500)
st_dat <- AddModuleScore(st_dat,features = list(smc.genes),name="smc.genes",n.bin = 500,ctrl.size=500)
st_dat <- AddModuleScore(st_dat,features = list(endo.genes),name="endo.genes",n.bin = 500,ctrl.size=500)
st_dat <- AddModuleScore(st_dat,features = list(Macrophage.genes),name="Macrophage.genes",n.bin = 500,ctrl.size=500)
saveRDS(st_dat,file=paste0(output_dir,prefix,".ST.Normalization.rds"))
write.table(st_dat@meta.data,file=paste(output_dir,prefix,".metadata.xls",sep=""),sep="\t",quote=F,row.names=F)

