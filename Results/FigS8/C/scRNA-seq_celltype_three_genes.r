library(ggplot2)
library(tidydr)
library(Seurat)
setwd('/Volumes/Expansion/zhang_et_al_2024_finaldata/Results/FigS8/C')
alldata <- readRDS("seurat_with_celltype_annotation_update.rds")
# check the factors
table(Idents(alldata))
# subset scRNA-seq data 
k <- c("Epidermis","Cortex","Endodermis","Pericycle","Phloem","Procambium","Xylem")
sub <- subset(alldata,idents=k)
table(Idents(sub))
# subset gene expression of NRT1.9, PP2A and asSOFL1
g <- VlnPlot(sub,features = "NPF2.9",pt.size=0.5,group.by = "celltype")+ ggplot2::ggtitle(label = expression(bolditalic("NRT1.9"))) +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.text = element_text(size = 18, face = "bold"),
          axis.title = element_text(size = 20,face="bold"),axis.title.x = element_blank(),legend.text = element_text(size = 16, hjust = 1, vjust = 1),)+ labs(
              y = "log-transformed expression")

g2 <- VlnPlot(sub,features = "AT1G26208",pt.size=0.5,group.by = "celltype")+ ggplot2::ggtitle(label = expression(bolditalic("NRT1.9"))) +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.text = element_text(size = 18, face = "bold"),
          axis.title = element_text(size = 20,face="bold"),axis.title.x = element_blank(),legend.text = element_text(size = 16, hjust = 1, vjust = 1),)+ labs(
              y = "log-transformed expression")

g3<- VlnPlot(sub,features = "PP2A3",pt.size=0.5,group.by = "celltype")+ ggplot2::ggtitle(label = expression(bolditalic("NRT1.9"))) +
    theme(plot.title = element_text(size = 24, face = "bold"),
          axis.text = element_text(size = 18, face = "bold"),
          axis.title = element_text(size = 20,face="bold"),axis.title.x = element_blank(),legend.text = element_text(size = 16, hjust = 1, vjust = 1),)+ labs(
              y = "log-transformed expression")
# create data frame
m <- cbind(g$data,g3$data$PP2A3,g2$data$AT1G26208)
names(m) <- c('NRT1.9','celltype','PP2A','asSOFL1')
# change the structure 
l<- melt(m, id.vars = "celltype", variable.name = "gene", value.name = "expression")
l$celltype <- as.character(l$celltype)

l$celltype <- as.factor(l$celltype)
l$celltype <- factor(l$celltype,levels=c("Epidermis", "Cortex", "Endodermis", "Pericycle", "Phloem", "Procambium", "Xylem"))

l$gene <- factor(l$gene,levels=c("NRT1.9","PP2A","asSOFL1"))

# draw figures with ggplot2 
p <- ggplot(l, aes(x=celltype, y=expression)) + 
    geom_violin(aes(fill=celltype),trim = FALSE, scale = "width")+geom_boxplot(width=0.1)+ scale_fill_manual(values=c("#6BB952", "#EC748B", "#C4A751", "#36ACA2", "#6EB1DE","#B67FB3",'#394A92'))+facet_grid(gene~., scales = 'free_y')
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
               panel.grid = element_blank(),legend.title  = element_blank(),strip.background.y = element_blank(),strip.text.y = element_text(size=16,face='bold.italic'),panel.spacing.y = unit(0,'cm'),
               legend.position = " ")+
    labs(title = "scRNA-seq",
         x = "Cell type", y = "log-transformed expression")

     ggsave(filename = "Figs9_scRNA-seq_celltype_violin.pdf", plot = p, width = 8, height = 6.5)
         # stastics analysis
test <- list()
k<- c()
for (i in c('NRT1.9','PP2A','asSOFL1')) {
    print(i)
    a1 <- aov(l[l$gene == i, ]$expression ~ l[l$gene == i, ]$celltype)
    test[[i]] <- TukeyHSD(a1) 
    m <- as.data.frame( test[[i]]$`l[l$gene == i, ]$celltype`)
m$Significance <- ifelse(m$`p adj` > 0.05, "ns", 
                      ifelse(m$`p adj` <= 0.05 & m$`p adj` > 0.01, "*", 
                      ifelse(m$`p adj` <= 0.01 & m$`p adj` > 0.001, "**",
                      ifelse(m$`p adj` <= 0.001 & m$`p adj` > 0.0001, "***",
                       "****"))))
k[[i]]<- m
#save
        write.xlsx(k[[i]],paste0('/Volumes/Expansion/zhang_et_al_2024_finaldata/Results/FigS8/C/scRNA-seq', i,'_statistics.xlsx'))
}
