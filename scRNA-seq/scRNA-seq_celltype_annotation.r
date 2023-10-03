# add all required library
# load the data dowload from embl-ebi
filter <- Read10X_h5("/Volumes/files/EvoCell_lab/Xue/scRNA-seq_reference_data/reference_data_aggregated_filtered_gene_bc_matrices.h5", use.names = TRUE, unique.features = TRUE)
# create a seurat object
filter <- CreateSeuratObject(filter)
filter
======================================================================================
# STEP1: quality control 
# check the percentage of mitochondra genes 
filter$mitoRatio <- PercentageFeatureSet(object = filter, pattern = "M")
filter$mitoRatio <- filter@meta.data$mitoRatio / 100
# Add number of genes per UMI for each cell to metadata
filter$log10GenesPerUMI<-log10(filter$nFeature_RNA) /log10(filter$nCount_RNA)
min(filter$log10GenesPerUMI)
max(filter$log10GenesPerUMI)
# add package to be able to use pip function 
# reorganize the metadata
metadata  <-filter@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%  dplyr::rename(sample = orig.ident, nUMI = nCount_RNA,nGene = nFeature_RNA)
filter@meta.data <- metadata
# check some parameters
median(filter@meta.data$nUMI) # 37095
mean(filter@meta.data$nUMI) # 51392.38 
median(filter@meta.data$nGene) #6781
# this means the genes per cell is consistent with literature but not the reads( cell filtering part is done). we need to filter reads at gene level
# check the distribution of avergage count using scater package 
counts <- GetAssayData(object = filter, slot = "counts")
nonzero <- counts > 0
ave.counts <- rowMeans(counts)
# decide which value as treshold
hist(log10(ave.counts), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))
#The filter threshold should cut the distribution at some point along the rectangular component to remove the majority of low-abundance genes.
abline(v=log10(0.006), col="blue", lwd=2, lty=2)
keep <- ave.counts >= 0.006
sum(keep) # 20556
filtered_counts <- counts[keep, ]
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filter@meta.data)
filtered_seurat ~# 20556 genes across 5145 cells 
# nomalization, feature selection, dimension reduction 
data.filt <- filtered_seurat
data.filt <-NormalizeData(data.filt)# log transforme
# check the top genes 
par(mar =c(4, 8, 2, 1))
C <-data.filt@assays$RNA@counts
C <-Matrix::t(Matrix::t(C)/Matrix::colSums(C)) *100
most_expressed <-order(apply(C, 1, median), decreasing =T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex =0.1, las =1, xlab ="% total count per cell",    col =(scales::hue_pal())(20)[20:1], horizontal =TRUE)

======================================================================================================================================
# STEP2 feature selection, Dimension reduction 
data.filt<-FindVariableFeatures (data.filt, verbose =F)
ggData<-data.filt@assays$RNA@meta.features
p1 <- ggplot(ggData, aes(vst.mean, vst.variance, color = vst.variable)) + 
    geom_point() + scale_color_manual(values = c("black", "firebrick")) + 
    geom_line(aes(vst.mean, vst.variance.expected), color = "dodgerblue") + 
    xlab("Average Expression") + ylab("Gene Variance") + 
    scale_x_log10() + scale_y_log10()
  p2 <- ggplot(ggData, aes(vst.mean, vst.variance.standardized, 
                         color = vst.variable)) + 
    geom_point() + scale_color_manual(values = c("black", "firebrick")) + 
    xlab("Average Expression") + ylab("Standardized Variance") + 
    scale_x_log10() +  theme(legend.position = "none")
p2 <- LabelPoints(plot = p2, repel = TRUE, points = VariableFeatures(data.filt)[1:10])
ggsave(p1 + p2, width = 10, height = 14, filename = "/Users/xung0001/Desktop/basicHVG.png")
# dimension reduction 
# scale the data 
data.filt <- ScaleData(data.filt, vars.to.regress = c("mitoRatio", "nGene"),
    assay = "RNA")
 # run PCA
data.filt <- RunUMAP(data.filt, reduction = "pca", dims = 1:30, n.components = 2, n.neighbors = 30,
    n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
# visualize PCA components
plot_grid(ncol = 3, DimPlot(data.filt, reduction = "pca", dim=1:2),DimPlot(data.filt, reduction = "pca", dim=3:4),DimPlot(data.filt, reduction = "pca", dim=5:6))
 VizDimLoadings(data.filt, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)
# we can plot the amount of variance explained by each PC
ElbowPlot(alldata, reduction = "pca", ndims = 50)
# run t-SNE
data.filt <- RunTSNE(data.filt, reduction = "pca", dims = 1:30, 
                   perplexity=30,
                   max_iter=1000,
                   theta=0.5,
                   eta=200,
                   num_threads=0 )
library(cowplot)
# Visualize tsne
p1 <- DimPlot(data.filt, reduction = "tsne", pt.size = 0.1, shuffle = TRUE)+ coord_fixed()
p1
# run uMAP
data.filt<- RunUMAP(data.filt, reduction = "pca", dims = 1:30, n.components = 2, n.neighbors = 30,
      n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
data.filt <- RunUMAP(data.filt, reduction.name = "UMAP10_on_PCA", reduction = "pca",
    dims = 1:30, n.components = 10, n.neighbors = 30, n.epochs = 200, min.dist = 0.3,learning.rate = 1, spread = 1)
    # visualize
    
plot_grid(ncol = 3, DimPlot(data.filt, reduction = "umap") +
              ggplot2::ggtitle(label = "UMAP_on_PCA"), DimPlot(data.filt, reduction = "UMAP10_on_PCA",
                                dims = 1:2) + ggplot2::ggtitle(label = "UMAP10_on_PCA"),
          DimPlot(data.filt, reduction = "UMAP10_on_PCA",dims = 3:4) +
              ggplot2::ggtitle(label = "UMAP10_on_PCA"))
 # visualize some marker genes using featureplot 
p <- c("APL", "AT1G29520", "DOF2.4", "OBP2", "NPF7.3", "AT5G01740","SMXL5","ATHB-8") #phloem/ pericycle/ procabium genes 
plot_list <- list()
for (i in seq_along(p) ){
    plot_list[[i]] <- FeaturePlot(data.filt, reduction = "tsne", dims = 1:2, features = p[i],
                                  ncol = 3, order = T) + NoLegend() + NoAxes() + NoGrid()+ ggplot2::ggtitle(label =p[i])
}
plot_grid(plotlist = plot_list, ncol = 3)
saveRDS(data.filt, "/Volumes/files/EvoCell_lab/Xue/scRNA-seq_reference_data/E-GEOD/h5_qc_dr.rds")

=====================================================================================================================
# STEP3 Clustering 
# load data 
alldata <- readRDS("/Volumes/files/EvoCell_lab/Xue/scRNA-seq_reference_data/E-GEOD/h5_qc_dr.rds")
# determine the k-nearest neighbor graph
alldata <- FindNeighbors(alldata, dims = 1:30, k.param = 60, prune.SNN = 1/15)
# determine the clusters for various resolution
res.used=seq(0.1,1,by=0.2)
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
DimPlot(alldata,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
 DimPlot(alldata,
        reduction = "tsne",
        label = TRUE,
        label.size = 6)
        # identify marker gene
# visualize clusters
p<-list( )
        p <- map(c(levels(Idents(alldata))),function(x){DimPlot(alldata,cells.highlight = CellsByIdentities(object = alldata,idents = x))})
        plot_grid(plotlist = p)

markers_genes <- FindAllMarkers(alldata, log2FC.threshold = 0.2, test.use = "wilcox",
                                min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50,
                                assay = "RNA")
markers_genes %>%    group_by(cluster) %>%    top_n(-20, p_val_adj) ->top20
markers_genes %>%    group_by(cluster) %>%    top_n(-5, p_val_adj) ->top5
# visualization top5 genes using heatmap, violin plot, dotplot 
alldata <- ScaleData(alldata, features = as.character(unique(top5$gene)), assay = "RNA")
DoHeatmap(alldata, features = as.character(unique(top5$gene)), group.by = sel.clust,
    assay = "RNA")
DotPlot(alldata, features = rev(top5$gene), group.by = sel.clust,
        assay = "RNA") + coord_flip()
        VlnPlot(alldata,features = top5$gene,pt.size=0.5,ncol = 5, group.by = sel.clust,assay = "RNA")

 ========================================================================================================= 
 # STEP4: cell type annotation 
# check the markers gene location in the clusters 
#we need to create a list with all the marker gene 
# do this to each individual cluster
for (i in names(marker.list)) {
    # Get the marker genes for the current cell type
    markers <- marker.list[[i]]
    save_path <- '/Volumes/files/EvoCell_lab/Xue/scRNA-seq_reference_data/E-GEOD'
    file_name <- paste0(i, "_DotPlot_gene_plots.pdf")
    VlnPlot <- VlnPlot(alldata, features = rev(markers), group.by = sel.clust, assay = "RNA") +
        ylab('') +
        xlab('')
    ggsave(filename = file.path(save_path, file_name),
           plot = VlnPlot, width = 18, height =8, units = "in", dpi = 300)   
}
# after I see where all the markees express 
# annotate cell type 
  alldata@meta.data$celltype <- Idents(alldata)
alldata$celltype <- recode(alldata$celltype,
"0"="LRC","1"="Epidermis+Atrichoblast","2"="Pericycle","3"="Pericycle","4"="Initials","5"="Phloem","6"="Cortex",
"7"="Procambium","8"="Endodermis","9"="Xylem","10"="LRC","11"="LRC","12"="Dividing","13"="Columella","14"="Trichoblast")

# subset data cluster1
clusters_to_subset <- 1 
subset_cells <- which(alldata$RNA_snn_res.0.7%in% clusters_to_subset)
cluster_seurat <- subset(alldata,cells=subset_cells)
cluster_seurat
cluster_seurat <- FindNeighbors(cluster_seurat, dims = 1:30, k.param = 60, prune.SNN = 1/15)
res.used=seq(0.1,1,by=0.1)
for (res in res.used) {
    cluster_seurat <- FindClusters(cluster_seurat, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
 clustree(cluster_seurat)+theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low="grey80",high="red")#or clustree(alldata@meta.data, prefix = "RNA_snn_res.")
sel.clust = sel.clust = "RNA_snn_res.0.2"
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
markers_genes <- FindAllMarkers(cluster_seurat, log2FC.threshold = 0.2, test.use = "wilcox",
                                min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50,
                                assay = "RNA")
markers_genes %>%    group_by(cluster) %>%    top_n(-20, p_val_adj) ->top20

# load the data from article 
library(readxl)
all_sheet_names <- excel_sheets("/Volumes/files/EvoCell_lab/Xue/scRNA-seq_reference_data/E-GEOD/xue_tables1.xlsx")
# make a list to store the interested cell types
data.list<- list()
p<- c("LRC","Atrichoblast","Trichoblast","Columella","Epidermis" ,"QC","Dividing","Cortex","Endodermis","Initials","Pericycle","Phloem","Procambium","Xylem")
for (i in p) {
 data<- as.data.frame(read_excel("/Volumes/files/EvoCell_lab/Xue/scRNA-seq_reference_data/E-GEOD/xue_tables1.xlsx", sheet = i))
 data.list[[i]]<- data
  names(data.list[[i]])[1] <- "gene"
 # only keep gene symbol
   data.list[[i]]$gene <- substr(data.list[[i]]$gene ,11,nchar(data.list[[i]]$gene))
}
# check the structure
summary(data.list)
## check the top markers they used in the paper…
gene_list <- top20$gene
# create a list to restore the subset data
extracted_data <- list()
final_data <- data.frame()  # Initialize an empty data frame instead of a vector
for (i in gene_list) {
    extracted_data_gene <- lapply(data.list, function(df) {
        subset(df, gene == i)  # Use == for exact string matching
    })
    
    # Check if any data frame in the list is not empty before storing it in 'extracted_data'
    non_empty_data <- Filter(function(df) nrow(df) > 0, extracted_data_gene)
    
    if (length(non_empty_data) > 0) {
        extracted_data[[i]] <- non_empty_data
        genedata <- extracted_data[[i]]
        gene_data <- bind_rows(genedata, .id = "cell_type") 
        gene_data_filtered <- gene_data %>%
            filter(avg_logFC > 0)
        final_data <- bind_rows(final_data, gene_data_filtered)  # Append data frames using bind_rows
    }
}

cluster <- c()
for (i in final_data$gene) { if (i %in% markers_genes$gene) {
    cluster <- rbind(cluster,data.frame(markers_genes[markers_genes$gene == i, ]))
} 
}
final_data <- cbind(final_data, cluster)
final_data$cell_type <- as.factor(final_data$cell_type)
levels(final_data$cell_type)
final_data<- final_data[order(final_data[11]), ]
View(final_data)
# some explore 
## res=0.3, I found subcluster 5 as Columella:626 ,4 as dividing cells 477CELLS,3 as endodermis cells (397cells), 2as Procabium cells:759CELLS  1mightbe epidermis
top20[top20$cluster=="2",]$gene %in% data.list$Procambium[data.list$Procambium$avg_logFC>0,]$gene
# res=0.4, not make dig change dispite we are able to differentiate Columella, dividing cells, endodermis cells, procabium cells 
# subset clusters_to_subset <- c(7,8,10,11,12,13) : res=0.2, 1: procambium cells:726  0:LRC mix, 3:dividing cells 465, 2,endodermis 349 top20identical 4:columella651
length(which(markers_genes[markers_genes$cluster=="1",]$gene%in% data.list$Procambium[data.list$Procambium$avg_logFC>0,]$gene))
# check if most genes are included in the excel file from publications in our generated data  especially top ones 
# or use commands 
subset(data.list$LRC[data.list$LRC$avg_logFC>0,], gene %in% markers_genes[markers_genes$cluster == "2",]$gene)
data.list$Procambium[which(markers_genes[markers_genes$cluster=="2",]$gene%in% data.list$Procambium[data.list$Procambium$avg_logFC>0,]$gene),]
##### we found where the cluster belongs
top20[top20$cluster=="1",]$gene %in% data.list$Atrichoblast[data.list$Atrichoblast$avg_logFC>0,]$gene
 # [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
top20[top20$cluster=="0",]$gene %in% data.list$Epidermis[data.list$Epidermis$avg_logFC>0,]$gene
 # [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
# annotate the subset cell type
cluster_seurat@meta.data$celltype <- Idents(cluster_seurat)
cluster_seurat$celltype <- recode(cluster_seurat$celltype,"0"="Atrichoblast","1"="Epidermis")
# change the cell type of cluster1 in dalldata
 h@meta.data$celltype<- as.character(alldata@meta.data$celltype)
cluster_seurat@meta.data$celltype<- as.character(cluster_seurat@meta.data$celltype)
for (i in seq_len(nrow(cluster_seurat@meta.data))) {
  id <- cluster_seurat@meta.data$cells[i]
  alldata@meta.data[alldata@meta.data$cells == id, "celltype"] <- cluster_seurat@meta.data$celltype[i]
}
# change cell type again as factor
alldata@meta.data$celltype <- as.factor(alldata@meta.data$celltype)
## do the same for cluster 2,10,12 since we found QC markers locate in those clusters
clusters_to_subset <- c(2,10,12)
subset_cells <- which(alldata$RNA_snn_res.0.7%in% clusters_to_subset)
cluster_seurat <- subset(alldata,cells=subset_cells)
cluster_seurat

cluster_seurat <- FindNeighbors(cluster_seurat, dims = 1:30, k.param = 60, prune.SNN = 1/15)
res.used=seq(0.1,1,by=0.1)
for (res in res.used) {
    cluster_seurat <- FindClusters(cluster_seurat, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
 clustree(cluster_seurat)+theme(legend.position = "bottom")+scale_color_brewer(palette = "Set1")+scale_edge_color_continuous(low="grey80",high="red")#or clustree(alldata@meta.data, prefix = "RNA_snn_res.")
sel.clust = sel.clust = "RNA_snn_res.0.7"
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
markers_genes <- FindAllMarkers(cluster_seurat, log2FC.threshold = 0.2, test.use = "wilcox",
                                min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50,
                                assay = "RNA")
markers_genes %>%    group_by(cluster) %>%    top_n(-20, p_val_adj) ->top20
# we found 4 QC, 2 dividing cells ,0,3 pericycle 1:LRC
cluster_seurat@meta.data$celltype <- Idents(cluster_seurat)
cluster_seurat$celltype <- recode(cluster_seurat$celltype,"0"="Pericycle","1"="LRC","2"="Dividing","3"="Pericycle","4"="QC")
# change the cell type of cluster1 in dalldata
 h@meta.data$celltype<- as.character(alldata@meta.data$celltype)
cluster_seurat@meta.data$celltype<- as.character(cluster_seurat@meta.data$celltype)
for (i in seq_len(nrow(cluster_seurat@meta.data))) {
  id <- cluster_seurat@meta.data$cells[i]
  alldata@meta.data[alldata@meta.data$cells == id, "celltype"] <- cluster_seurat@meta.data$celltype[i]
}
# change cell type again as factor
alldata@meta.data$celltype <- as.factor(alldata@meta.data$celltype)
=============================================================
# after we found all the clusters , we identify all markers for each cluster, check if everything is correct
 # now I want to change the order of my clusters in seurat object 
        alldata@meta.data$celltype <- factor(alldata@meta.data$celltype, 
                            levels=c("Trichoblast",
                                 "Atrichoblast",
                                   "Epidermis",
                                "LRC",
                                  "Columella", 
                                "Endodermis",
  "Cortex",  "Procambium", "Pericycle", "Phloem", 
 "Xylem", "Dividing", "Initials", "QC"))

# generate dotplot for all the markers 
for (i in names(marker.list)) {
    # Get the marker genes for the current cell type
    markers <- marker.list[[i]]
    save_path <- '/Volumes/files/EvoCell_lab/Xue/scRNA-seq_reference_data/E-GEOD'
    file_name <- paste0(i, "_DotPlot_gene_plots.pdf")
   dotplot <- DotPlot(alldata, features = rev(markers), group.by = "celltype", assay = "RNA") +
        ylab('') +
        xlab('')
    ggsave(filename = file.path(save_path, file_name),
           plot = dotplot, width = 18, height =8, units = "in", dpi = 300)   
}
# after validate all the information, everything is correct.
DimPlot(alldata, group.by = "celltype", label = TRUE, repel = TRUE) + NoAxes()
saveRDS(alldata, "/Volumes/files/EvoCell_lab/Xue/scRNA-seq_reference_data/seurat_with_celltype_annotation.rds")



























