setwd("/Volumes/Expansion/zhang_et_al_2024_finaldata/Fig3_asSOFL1_subcellular/C")
# load packages
library(ggplot2)
library(readxl)
library(dplyr)
library(ggsignif)
library(xlsx)
# load data 
data <- read_excel("counts_per_celltype_root5_con_sr_asSOFL1-2.xlsx")
head(data)
# modify the structure of data 
names(data) <- c("replicates",'celltype',"transcripts")
data$replicates <- as.factor(data$replicates)
levels(data$replicates)
data$celltype <- factor(data$celltype,levels = c("Epidermis","Cortex","Endodermis","Pericycle","Phloem","Procambium","Xylem"))
# draw figure
p <- ggplot(data, aes(x = celltype, y = transcripts)) + 
    geom_violin(aes(fill = celltype),trim = FALSE,scale='width') +  #geom_jitter(shape = 16, size = 2, position = position_jitter(0.2)) +
    geom_boxplot(width = 0.1) +
  
    scale_fill_manual(values=c("#989B7A", "#E2E8AA", "#5D9C76", "#A0D5B2", "#684FA1","#D1EBE7",'#3A4A92'))
# or  scale_color_manual(values = c("#1F77B4B2", "#FF7F0EB2", "#2CA02CB2", "#d62728B2", "#9467bdb2", "#BC564BB2", '#E377C2B2'))

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
    guides(color = guide_legend(override.aes = list(size = 2))) # Adjust the size value

p <- p + labs(title = expression(bolditalic("asSOFL1")),
              x = "Cell type", y = "Transcript count")+ylim(0,4)

p
# add statistics 
a1<-aov(data$transcripts~ data$celltype)
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
#save
write.xlsx(m,'/Volumes/Expansion/zhang_et_al_2024_finaldata/Fig3_asSOFL1_subcellular/C/asSOFL1_statistics.xlsx')

ggsave("Fig3_asSOFL1_boxplot_Violin_root5_sr.pdf",p,width = 5.54, height = 6.6)

