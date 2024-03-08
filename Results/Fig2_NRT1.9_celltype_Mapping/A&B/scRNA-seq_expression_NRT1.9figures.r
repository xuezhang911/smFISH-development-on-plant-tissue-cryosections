# NRT1.9 was expressed in 3050 cells with median value median(NRT@assays$RNA@counts) =3 NRT is seurat object including only NRT1.9 gene. 
library(ggplot2)
library(tidydr)
library(Seurat)
setwd("/Volumes/Expansion/zhang_et_al_2024_finaldata/Fig2_NRT1.9_celltype_Mapping/A&B")
alldata <- readRDS("seurat_with_celltype_annotation_update.rds")
# check the factors
table(Idents(alldata))
# subset scRNA-seq data 
k <- c("Epidermis","Cortex","Endodermis","Pericycle","Phloem","Procambium","Xylem")
sub <- subset(alldata,idents=k)
table(Idents(sub))
# reorder celltype
sub@meta.data$celltype<- factor(sub@meta.data$celltype,levels=c("Epidermis","Cortex","Endodermis","Pericycle","Phloem","Procambium","Xylem"))
# visualize our gene of interest
genes <- c("NPF2.9")

# Initialize an empty list to store data frames for each gene
gene_data <- list()

# Loop over each gene
for (gene in genes) {
    # Extract gene expression data for the current gene
    k <- sub@assays$RNA@data[gene,]
    k <- t(t(k))
    k <- as.data.frame(k)
    
    # Store the data frame in the gene_data list
    gene_data[[gene]] <- k
}
# Combine the data frames for all genes
combined_data <- do.call(cbind, gene_data)
# Add celltype information
combined_data$"celltype" <- sub@meta.data$celltype
# Rename columns
colnames(combined_data) <- c(genes, "celltype")
# Display the head of the combined data frame
head(combined_data)
# add the embeddings information 
head(sub@reductions$tsne@cell.embeddings)
# see if we can directly cbind data 
which(!names(sub@reductions$tsne@cell.embeddings)==names(combined_data))

# combine the tsne information 
 df <- cbind(combined_data,sub@reductions$tsne@cell.embeddings)
# make celltype as factor
df$celltype <- as.factor(df$celltype)
# change name
names(df)[1] <- c('NRT1.9')
# draw figure
ggplot(df, aes(x = tSNE_1, y = tSNE_2, colour = NRT1.9)) +
    geom_point(size = 0.3, alpha = 1) +
    scale_colour_gradientn(colours = c("#779fd3","lightgrey", "#a13037"), limits = c(0, 3),breaks=seq(0,3,0.5), oob = scales::squish)
       # add lable to the plot 
    library(dplyr)
    library(ggrepel)
# calculate the place where to put label 
label.df_2 <- df %>% 
     group_by(celltype) %>% 
  summarize(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2))

p <- ggplot(df, aes(x = tSNE_1, y = tSNE_2, colour = NRT1.9)) +
    geom_point(size =1, alpha = 1) +
    scale_colour_gradientn(colours = c("#779fd3","lightgrey", "#a13037"), limits = c(0, 4),breaks=seq(0,4,1), oob = scales::squish) +
    ggrepel::geom_label_repel(data = label.df_2, aes(x = tSNE_1, y = tSNE_2, label = celltype), 
                              colour = "black", inherit.aes = FALSE,label.size = NA,fill=NA,size=8
    )+ tidydr::theme_dr()+theme(plot.title = element_text(size = 25,hjust=0.5),legend.title= element_text(size = 16,hjust=1,vjust = 1),legend.text= element_text(size = 14,hjust=1,vjust = 1),
                                axis.text = element_blank(),
                                axis.title = element_text(size = 20,vjust=0.5),panel.grid = element_blank())
# save as figure
ggsave(filename = "Fig2_scRNA_feature_plot.pdf", plot = p, width = 8, height = 6)

    # we can also extract information and use ggplot2 to draw violin plot 
# violin plot+ box plot 
p <- ggplot(df, aes(x=celltype, y=NRT1.9)) + 
    geom_violin(aes(fill=celltype),trim = T, scale = "width")+geom_boxplot(width=0.1)+
    scale_fill_manual(values=c("#6BB952", "#EC748B", "#C4A751", "#36ACA2", "#6EB1DE","#B67FB3",'#394A92'))
#scale_fill_manual(values=c("#A68BC2", "#85C680", "#D09B7E", "#79AED2", "#EFB266","#DCCA83",'#EE7677'))
#scale_fill_manual(values=c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2","#8491B4B2",'#91D1C2B2'))

p <- p + theme(panel.grid.major = element_blank(),axis.line = element_line(colour = "black"), 
               panel.background = element_blank(),
               axis.text = element_text(size = 20),
               axis.title = element_text(size =20),
               axis.title.x = element_blank(),  # Remove x-axis title
               axis.text.x = element_text(size = 20, angle = 35, hjust = 1, vjust = 1),  # Set the size and angle of x-axis text
               plot.title = element_text(colour = "black",size = 25, face = "bold",hjust = 0.5, vjust = 0.5),
               legend.text = element_text(size = 16 ),
               panel.grid = element_blank(),legend.title  = element_blank(),
               legend.position = " ")+
    labs(title = expression(bold('scRNA-seq')),
         x = "Cell type", y = "log-transformed expression")
        
## add statistics
a1<-aov(df$NRT1.9~ df$celltype)
summary(a1)
test<-TukeyHSD(a1)
test
# add significance in
m <- as.data.frame(test$`df$celltype`)
m$Significance <- ifelse(m$`p adj` > 0.05, "ns", 
                      ifelse(m$`p adj` <= 0.05 & m$`p adj` > 0.01, "*", 
                      ifelse(m$`p adj` <= 0.01 & m$`p adj` > 0.001, "**", 
                      ifelse(m$`p adj` <= 0.001 & m$`p adj` > 0.0001, "***", "****"))))

# Printing the updated data frame
print(m)
## optional
 library(xlsx)
library(ggsignif)
write.xlsx(m,'/Volumes/Expansion/zhang_et_al_2024_finaldata/Fig2_NRT1.9_celltype_Mapping/A&B/scrna_statistics.xlsx')


p <- p+geom_signif(annotations="***",y_position = 3,xmin=3.9,xmax=4.8,tip_length = 0.01,textsize = 6)
  p <- p+geom_signif(annotations="***",y_position = 3,xmin=5.3,xmax=6,tip_length = 0.01,textsize = 6)
   p <- p+geom_signif(annotations="*",y_position = 3.6,xmin=3.9,xmax=6,tip_length = 0.01,textsize = 6)+ylim(0,4)
   ggsave("Fig2_scRNA-seq_Violin.pdf",p,width = 8, height = 6)