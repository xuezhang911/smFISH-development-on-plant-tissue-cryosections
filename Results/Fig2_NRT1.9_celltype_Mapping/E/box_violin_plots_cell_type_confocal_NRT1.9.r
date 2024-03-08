setwd("/Volumes/Expansion/zhang_et_al_2024_finaldata/Fig2_NRT1.9_celltype_Mapping/E")

# load packages
library(ggplot2)
library(readxl)
library(ggsignif)
library(dplyr)

# load  the single z-stack confocal data : 
#  data <- read_excel("NRT19_cryosections_by_celltype-sub.xlsx")
# load maxium projection, dapi, and WF mixed mannually quantified data 
# data <- read_excel("development_COUNTS_NRT1.9_celltypes.xlsx")
# Or load max projection  
data <- read_excel("counts_celltype_MAX_s8_SR_NRT1_CDPK1_CON5-1_SR4_ch1.xlsx")
head(data)
# modify the structure of data 
names(data) <- c("replicates",'celltype',"transcripts")
data$replicates <- as.factor(data$replicates)
# levels(data$replicates) <- c("rep1","rep2","rep3","rep4","rep5")
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

# make figures 

a1<-aov(data$norm~ data$celltype)
summary(a1)
test<-TukeyHSD(a1)
test

# Violin plus box plot 
# Violin plus box plot 
p <- ggplot(data, aes(x = celltype, y = norm)) + 
    geom_violin(aes(fill = celltype),trim = F,scale='width',adjust=1) +
    geom_boxplot(width = 0.1)  +
    scale_fill_manual(values = c("#989B7A", "#E2E8AA", "#5D9C76", "#A0D5B2", "#684FA1","#D1EBE7",'#3A4A92'))

p <- p + theme(panel.grid.major = element_blank(), axis.line = element_line(colour = "black"), 
               panel.background = element_blank(),
               axis.text = element_text(size = 20),
               axis.title = element_text(size = 20),
               axis.title.x = element_blank(),
               axis.text.x = element_text(size = 20, angle = 35, hjust = 1, vjust = 1),
               plot.title = element_text(colour = "black", size = 25, face = "bold", hjust = 0.5, vjust = 0.5),
               legend.text = element_text(size = 16),
               legend.title = element_blank(),
               legend.position = '') +
    guides(color = guide_legend(override.aes = list(size = 2)))  # Adjust the size value

p <- p + labs(title = expression(bold('cryo-smFISH')),
              x = "Cell type", y = "Transcript count per cell")

p
## add statistics
a1<-aov(data$norm~ data$celltype)
summary(a1)
test<-TukeyHSD(a1)
test
# add significance in
m <- as.data.frame(test$`data$celltype`)
m$Significance <- ifelse(m$`p adj` > 0.05, "ns", 
                      ifelse(m$`p adj` <= 0.05 & m$`p adj` > 0.01, "*", 
                      ifelse(m$`p adj` <= 0.01 & m$`p adj` > 0.001, "**", 
                      ifelse(m$`p adj` <= 0.001 & m$`p adj` > 0.0001, "***", "****"))))

# Printing the updated data frame
print(m)
## optional
 library(xlsx)
library(ggsignif)
write.xlsx(m,'/Volumes/Expansion/zhang_et_al_2024_smFISH_cryosections/smFISH_in_cryosections/Final_data_each_figure/Fig2_NRT1.9_celltype_Mapping/E/smFISH_statistics.xlsx')

# add statistics for max projection
p <- p+geom_signif(annotations="*",y_position = 15, xmin=4.1,xmax=4.7,tip_length = 0.01,textsize = 6)
p <- p+geom_signif(annotations="ns",y_position = 26, xmin=4.1,xmax=5.8,tip_length = 0.01,textsize = 6)+ylim(0,28)
p <- p+geom_signif(annotations="***",y_position = 15, xmin=5.3,xmax=5.9,tip_length = 0.01,textsize = 6)
ggsave(filename = "Fig2_smFISH_celltype_violin.pdf", plot = p, width = 6, height = 6)