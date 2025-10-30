#untar("GSE117570_RAW.tar")
getwd()

fs=list.files('./')
fs
library(tidyverse)
library(data.table)
library(Seurat)
sceList = lapply(fs,function(pro){ 
  # pro=samples[1] 
  print(pro) 
  sce=CreateSeuratObject( Read10X(file.path(dir,paste0(pro,"/filtered_feature_bc_matrix"))), 
                          project = pro,
                          min.cells = 5,
                          min.features = 300 ) 
  return(sce)
})

sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = fs)

sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
### 该操作会在metadata数据里面增加一列叫做percent.mt
### 质控数据可视化，使用VlnPlot函数
#pdf(file="normal单细胞质控.pdf",width=10,height=12)
#VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
sce.all <- subset(sce.all, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
sce.all <- NormalizeData(sce.all)
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst", nfeatures = 2000)
sce.all <- ScaleData(sce.all, features = rownames(sce.all))
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all),reduction.name = "pca")
ElbowPlot(sce.all)
library(harmony)
sce.all <- RunHarmony(sce.all,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
sce.all <- FindNeighbors(sce.all, reduction = "harmony", dims = 1:15)
### 设置多个resolution选择合适的resolution
sce.all <- FindClusters(sce.all, resolution = 1)
sce.all <- RunUMAP(sce.all, reduction = "harmony", dims = 1:15,reduction.name = "umap")





DimPlot(sce.all,group.by = "seurat_clusters",label = T,cols = cols)


Idents(sce.all) = sce.all$celltype
levels(sce.all)
markers = FindAllMarkers(sce.all,only.pos = T)
top10 = markers %>% 
  group_by(cluster) %>% 
  top_n(10,wt = avg_log2FC)
