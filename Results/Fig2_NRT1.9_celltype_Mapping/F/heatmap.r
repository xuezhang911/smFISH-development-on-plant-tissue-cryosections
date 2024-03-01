# NRT1.9 was expressed in 3050 cells with median value median(NRT@assays$RNA@counts) =3 NRT is seurat object including only NRT1.9 gene. 

library(ggplot2)
library(tidydr)
library(Seurat)
setwd('/Volumes/Expansion/zhang_et_al_2024_finaldata/Fig2_NRT1.9_celltype_Mapping/A&B')
alldata <- readRDS("seurat_with_celltype_annotation_update.rds")
# check the factors
table(Idents(alldata))
# subset scRNA-seq data 
k <- c("Epidermis","Cortex","Endodermis","Pericycle","Phloem","Procambium","Xylem")
sub <- subset(alldata,idents=k)
table(Idents(sub))
# reorder celltype
sub@meta.data$celltype<- factor(sub@meta.data$celltype,levels=c("Epidermis","Cortex","Endodermis","Pericycle","Phloem","Procambium","Xylem"))

    k1 <- sub@assays$RNA@data["NPF2.9",]
    k1 <- t(t(k1))
    k1 <- as.data.frame(k1)
    # test if the row.names are equal with metadata
    row.names(k1)==sub@meta.data$cells
    # yes, it's true, create a data.frame
    k1$"celltype" <- sub@meta.data$celltype
    names(k1)[1] <- "expression"
    head(k1)
   
    # extract scRNA-seq
 sum_expression <- aggregate(expression ~ celltype, data = k1, FUN = median)
# extract smFISH data
setwd('/Volumes/Expansion/zhang_et_al_2024_finaldata/Fig2_NRT1.9_celltype_Mapping/E')

# load packages
library(ggplot2)
library(readxl)
library(ggsignif)
library(dplyr)
# load data
data <- read_excel("counts_celltype_MAX_s8_SR_NRT1_CDPK1_CON5-1_SR4_ch1.xlsx")
head(data)
# modify the structure of data 
names(data) <- c("replicates",'celltype',"transcripts")
data$replicates <- as.factor(data$replicates)
data$celltype <- as.factor(data$celltype)
levels(data$celltype)
data$celltype <- factor(data$celltype,levels = c("Epidermis","Cortex","Endodermis","Pericycle","Phloem","Procambium","Xylem"))
# normalize the dataset with cell quantity for each cell type

celltype_counts <- data %>%
    group_by(celltype) %>%
    summarize(total_cells = n())

    data <- data %>%
  left_join(celltype_counts, by = "celltype")
data$norm <- (data$transcripts/data$total_cells)*10
# extract data 
library(dplyr)
expression <-data %>%
 group_by(celltype) %>%
  summarize(total_cells = median(norm))
# create a dataframe including both data 
df <- cbind(sum_expression[,2],expression[,2])
colnames(df) <- c('scRNA-seq','smFISH')
row.names(df) <- c("Epidermis","Cortex","Endodermis","Pericycle","Phloem","Procambium","Xylem")

# draw figure :
library(pheatmap)
library(ComplexHeatmap)
names(df)[2] <- 'cryo-smFISH'

pdf("heatmap.pdf", width = 8, height = 10)
Heatmap(scale(df),cluster_rows = T,cluster_columns = F,col=colorRampPalette(c("#779fd3","lightgrey", "#a13037"))(100),rect_gp = gpar(col='lightgrey',lwd=1.5),name='value',width=unit(5,'cm'),height=unit(10,'cm'),heatmap_legend_param = list(at=c(-1,0,1,2),title='',gp=gpar(fontsize=14)),column_names_rot = 30,cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", df[i, j]), x, y, gp = gpar(fontsize = 12,col='blue'))},row_names_gp = gpar(fontsize = 16),column_names_gp = gpar(fontsize = 16))
    dev.off()