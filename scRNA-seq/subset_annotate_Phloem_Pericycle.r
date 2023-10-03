# load data 
alldata <- readRDS("/Volumes/files/EvoCell_lab/Xue/scRNA-seq_GSE141730/E-GEOD-152766/h5_qc_dr.rds")
# load library
suppressPackageStartupMessages({
    library(Seurat)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(rafalib)
    library(clustree)
})
# determine the k-nearest neighbor graph
alldata <- FindNeighbors(alldata, dims = 1:30, k.param = 60, prune.SNN = 1/15)
# determine the clusters for various resolution
res.used=seq(0.1,1,by=0.1)
for (res in res.used) {
    alldata <- FindClusters(alldata, graph.name = "RNA_snn", resolution = res.used, algorithm = 1)
}
# visualize with clusetre
clustree(alldata)+theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low="grey80",high="red") #or clustree(alldata@meta.data, prefix = "RNA_snn_res.")
# decide use res=0.7
# from this clustree we know res0.7 is good option
sel.clust = "RNA_snn_res.0.7"
alldata <- SetIdent(alldata, value = sel.clust)
table(alldata@active.ident)
# subset data cluster5 where phloem cells locate
clusters_to_subset <- 5
subset_cells <- which(alldata$RNA_snn_res.0.7%in% clusters_to_subset)
cluster_seurat <- subset(alldata,cells=subset_cells)
cluster_seurat
cluster_seurat <- FindNeighbors(cluster_seurat, dims = 1:30, k.param = 60, prune.SNN = 1/15)
res.used=seq(0.1,1,by=0.1)
for (res in res.used) {
    cluster_seurat <- FindClusters(cluster_seurat, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
 clustree(cluster_seurat)+theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low="grey80",high="red")#or clustree(alldata@meta.data, prefix = "RNA_snn_res.")
sel.clust = "RNA_snn_res.0.3"
cluster_seurat <- SetIdent(cluster_seurat, value = sel.clust)
table(cluster_seurat@active.ident)
DimPlot(cluster_seurat,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
 DimPlot(cluster_seurat,
        reduction = "tsne",
        label = TRUE,
        label.size = 6)
## check where the cc and protophloem markers locate
Protophloem<- c("AT1G11915","DOF5.6",'AT1G77700','COPT4','AT1G55760','PIN4','DOF2.4','CLE45','BAM3')
 CC <- c('AT1G04330','AT1G07407', 'SUC2','APL','HIPP21','AT4G14020','AT4G00880','OBP2','AT1G10380')
VlnPlot(cluster_seurat, features = rev(Protophloem), group.by = sel.clust, assay = "RNA")
VlnPlot(cluster_seurat, features = rev(CC), group.by = sel.clust, assay = "RNA")
library(dplyr)
# to be use top_n function 
markers_genes <- FindAllMarkers(cluster_seurat, log2FC.threshold = 0.2, test.use = "wilcox",
                                min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50,
                                assay = "RNA")
markers_genes %>%    group_by(cluster) %>%    top_n(-20, p_val_adj) ->top20

# add subcelltype annotation 
cluster_seurat@meta.data$celltype <- Idents(cluster_seurat)
cluster_seurat$celltype <- recode(cluster_seurat$celltype,"0"="MSE","1"="PSE","2"="Companion cells")
alldata <- readRDS("seurat_with_celltype_annotation_update.rds")
 alldata@meta.data$celltype<- as.character(alldata@meta.data$celltype)
cluster_seurat@meta.data$celltype<- as.character(cluster_seurat@meta.data$celltype)
for (i in seq_len(nrow(cluster_seurat@meta.data))) {
  id <- cluster_seurat@meta.data$cells[i]
  alldata@meta.data[alldata@meta.data$cells == id, "celltype"] <- cluster_seurat@meta.data$celltype[i]
}

# change cell type again as factor
alldata@meta.data$celltype <- as.factor(alldata@meta.data$celltype)
Idents(alldata)<-alldata$celltype
# DEG analysis for CC, PME, PSE
cluster.markers <- FindMarkers(alldata, ident.1 = "MSE", min.pct = 0.25)
write.csv(cluster.markers,file = file.path(save_path,paste0("MSE", "_markers.csv") ))

#Dotplot
dotplot <- DotPlot(alldata, features = rev(CC), group.by = "celltype", assay = "RNA") +
    ylab('') +
    xlab('')
    ggsave(filename = "Companion_dotplot.pdf",
       plot = dotplot, width = 18, height =8)  
       dotplot <- DotPlot(cluster_seurat, features = rev(Protophloem), group.by = "celltype", assay = "RNA") +
    ylab('') +
    xlab('')
    ggsave(filename = "Companion_phloem_subtype_dotplot
    .pdf",
       plot = dotplot, width = 18, height =8)  
############################subset PPP and XPP cells
       alldata <- readRDS("/Volumes/files/EvoCell_lab/Xue/scRNA-seq_GSE141730/E-GEOD-152766/seurat_with_celltype_annotation_update.rds")
       k <- subset(alldata,idents="Pericycle")
       # determine the k-nearest neighbor graph
k <- FindNeighbors(k, dims = 1:30, k.param = 60, prune.SNN = 1/15)
# determine the clusters for various resolution
res.used=seq(0.1,1,by=0.1)
for (res in res.used) {
   k <- FindClusters(k, graph.name = "RNA_snn", resolution = res.used, algorithm = 1)
}
# visualize with clusetre
clustree(k)+theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low="grey80",high="red") #or clustree(alldata@meta.data, prefix = "RNA_snn_res.")

sel.clust = "RNA_snn_res.0.3"
cluster_seurat <- SetIdent(k, value = sel.clust)
table(cluster_seurat@active.ident)

## check where the cc and protophloem markers locate
PPP<- c("NPF2.9","PHB6",'CALS8','AtbZIP6','MSE7','AT3G11930','SBT4.12','AT1G26450','SUC2','APL','NEN4','NAC086','TPPD','ZIP4') # cluster 2
XPP <- c('AT1G02460','AT4G30450', 'DOT1','XTH21','AT2G28780','CYP73A5','AT2G44300','GATA23','CEP5') # cluster 0,1
VlnPlot(cluster_seurat, features = rev(PPP), group.by = sel.clust, assay = "RNA")
VlnPlot(cluster_seurat, features = rev(XPP), group.by = sel.clust, assay = "RNA")
library(dplyr)

# add subcelltype annotation 
cluster_seurat@meta.data$celltype <- Idents(cluster_seurat)
cluster_seurat$celltype <- recode(cluster_seurat$celltype,"0"="XPP","1"="XPP","2"="PPP")

 alldata@meta.data$celltype<- as.character(alldata@meta.data$celltype)
cluster_seurat@meta.data$celltype<- as.character(cluster_seurat@meta.data$celltype)
for (i in seq_len(nrow(cluster_seurat@meta.data))) {
  id <- cluster_seurat@meta.data$cells[i]
  alldata@meta.data[alldata@meta.data$cells == id, "celltype"] <- cluster_seurat@meta.data$celltype[i]
}

# change cell type again as factor
alldata@meta.data$celltype <- as.factor(alldata@meta.data$celltype)
Idents(alldata)<-alldata$celltype
# DEG analysis for CC, PME, PSE
cluster.markers <- FindMarkers(alldata, ident.1 = "PPP", min.pct = 0.25)
write.csv(cluster.markers,file = file.path("/Volumes/files/EvoCell_lab/Xue/scRNA-seq_GSE141730/E-GEOD-152766",paste0("PPP", "_markers.csv") ))
cluster.markers <- FindMarkers(alldata, ident.1 = "XPP", min.pct = 0.25)
write.csv(cluster.markers,file = file.path("/Volumes/files/EvoCell_lab/Xue/scRNA-seq_GSE141730/E-GEOD-152766",paste0("XPP", "_markers.csv") ))