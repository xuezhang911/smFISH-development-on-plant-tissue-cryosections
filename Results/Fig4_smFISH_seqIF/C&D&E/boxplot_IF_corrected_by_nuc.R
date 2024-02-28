# set working directory
setwd("/Volumes/Expansion/zhang_et_al_2024_finaldata/Fig4_smFISH_seqIF/C")
# load libraries
library(readxl)
library(reshape2)
library(dplyr)
library(ggplot2)
#load data
df <- read.csv("RNA_protein_one_cell_mask_includingcelltype_CORRECTED_nuc_area_only.csv",sep=";",header=T)
head(df)
df[is.na(df)] <- 0
# load rna absolute count and protein abundance normalized with cell size 
df <- subset(df,select=c("RNA_transcripts","nuc_mean_prot_intensity","celltype"))
df$nuc_mean_prot_intensity <- df$nuc_mean_prot_intensity*100
names(df) <- c("RNA_counts","Protein_abundance","celltype")
head(df)

### histogram for protein
# check which cells have protein abundance greater than 0
which(df$Protein_abundance>0)
# create a new dataframe
m <- df[which(df$Protein_abundance>0),]
head(m)
# order the cells
m <- m[order(m$Protein_abundance),]
# draw histogram for protein 

p<- ggplot(m,aes(x=Protein_abundance))+geom_histogram(alpha=0.4, color="black", fill="gray", bins = 10)

p<-p+geom_vline(aes(xintercept = median(Protein_abundance)),col="black",size=1,linetype="dashed",)+
  xlab("Protein abundance per cell") + ylab("Frequency")

p<-p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             legend.title = element_blank(),
             plot.title = element_text(colour="black", size =20,face="bold",hjust = 0.5,vjust=0.5 ),
             legend.position = "",
             axis.text=element_text(size=20),
             axis.title=element_text(size=20))+ggtitle("H4Ac")

 ggsave(filename ='hist_H4Ac.pdf', plot = p, width = 6, height = 6)

# draw histogram for RNA
df <- df[order(df$RNA_counts),]
p<- ggplot(df,aes(x=RNA_counts))+geom_histogram(alpha=0.4, color="black", fill="gray", bins = 7)

p<-p+geom_vline(aes(xintercept = median(RNA_counts)),col="black",size=1,linetype="dashed",)+
    xlab("Number of transcripts per cell") + ylab("Frequency")

p<-p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             legend.title = element_blank(),
             plot.title = element_text(colour="black", size =20,face="bold.italic",hjust = 0.48 ),
             legend.position = "",
             axis.text=element_text(size=20),
             axis.title=element_text(size=20))+ggtitle("PP2A")

ggsave(filename ='hist_RNA_PP2A.pdf', plot = p, width = 6, height = 6)

######################################################################################
# normalize RNA and protein by cell size 
# extract cell size 
df <- read.table("Intensity_immunoFilterCells.txt",sep="\t",header=T)
# put cell size in the following file 
data <- read.csv("RNA_protein_one_cell_mask_includingcelltype_CORRECTED_nuc_area_only.csv",sep=";",header=T)
 # add cellsize to data
data$cellarea <- df$AreaShape_Area
 # get normalized RNA counts
 data$noRNA <- data$RNA_transcripts/data$cellarea
data$noRNA <- round(data$noRNA,6)*10000
# assign NA as 0
df[is.na(df)] <- 0
 # make protein normalized
 data$noProtein <- data$Protein_abundance/data$cellarea
data$noProtein<- data$noProtein*100
  # subset the data 
data<- subset(data,select=c("noRNA","noProtein","celltype"))
# rename the data 
names(data) <- c("RNA_counts","Protein_abundance","celltype")

data$celltype <- factor(data$celltype,levels=c("Epidermis","Cortex","Endodermis","Pericycle","Phloem","Procambium","Xylem"))


s <- melt(data, id.vars = "celltype", 
          measure.vars = c("RNA_counts", "Protein_abundance"),
          variable.name = "Measurement", 
          value.name = "Value")
          # change the levels name
levels(s$Measurement)[2] <- "H4Ac"
levels(s$Measurement)[1] <- "PP2A"
### draw plot 
 p <- ggplot(s, aes(x=celltype, y=Value)) + 
    geom_violin(aes(fill=celltype),trim = FALSE, scale = "width")+geom_boxplot(width=0.1)+  scale_fill_manual(values = c("#989B7A", "#E2E8AA", "#5D9C76", "#A0D5B2", "#684FA1","#D1EBE7",'#3A4A92'))+facet_grid(Measurement~.)
#
#

p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
                          panel.background = element_rect(fill=NA,color='black'),
                          axis.text = element_text(size = 20),
                          axis.title = element_text(size =20),
                          axis.title.x = element_blank(),  # Remove x-axis title
                          axis.text.x = element_text(size = 20, angle = 35, hjust = 1, vjust = 1),  # Set the size and angle of x-axis text
                          plot.title = element_text(colour = "black",size = 25, face = "bold",hjust = 0.5, vjust = 0.5),
                          legend.text = element_text(size = 16 ),
                          panel.grid = element_blank(),legend.title  = element_blank(),strip.background.y = element_blank(),strip.text.y = element_blank(),panel.spacing.y = unit(0.2,'cm'),
                          legend.position = " ")+
    labs(title = " ",  
         x = "Cell type", y = "Value")+  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))


##statistics 
# create a list to store the statistics results for both PP2A and H4Ac
data_list <- list()
# loop 
for (variable in c('PP2A','H4Ac')) 
 {a1 <- aov(Value ~ celltype, data = s[s$Measurement==variable,])
summary(a1)
test<-TukeyHSD(a1)
print(test)
h<- as.data.frame(test$celltype[,ncol(test$celltype)])
names(h)<- variable
h$Significance <- ifelse(h[,1] > 0.05, "ns", 
                      ifelse(h[,1] <= 0.05 & h[,1] > 0.01, "*", 
                      ifelse(h[,1] <= 0.01 & h[,1] > 0.001, "**", 
                      ifelse(h[,1] <= 0.001 & h[,1] > 0.0001, "***", '****'))))
                      print(h)
data_list[[variable]]<- h
}
 # combine into a single data frame
    data <-cbind(data_list[[1]],data_list[[2]])
# save data
write.xlsx(data,'/Volumes/Expansion/zhang_et_al_2024_finaldata/Fig4_smFISH_seqIF/C/Fig4_smFISH_seqIF_statistics.xlsx')


## we make use of the following format to draw pvalue in the figure
library(rstatix)
# only would like to use this format, the significance is annova test
stat.test <- s %>%
    group_by(Measurement) %>%
    t_test(Value ~ celltype)
    # add xy position 
stat.test <- stat.test %>% add_xy_position(x = "celltype",  step.increase = 0.12,)
stat.test<- stat.test %>%
    mutate(group1group2 = paste(group2, group1, sep = "-"))
# change p-value and significance with annova 
# check if group name is the same as data 
which(!stat.test$group1group2==rownames(data))
#place p value
stat.test$p.adj <- c(data$PP2A,data$H4Ac)

# Redefine p.adj.signif based on p.adj
stat.test <- stat.test %>%
  mutate(p.adj.signif = ifelse(p.adj > 0.05, 'ns',
                             ifelse(p.adj <= 0.05 & p.adj > 0.01, '*',
                             ifelse(p.adj <= 0.01 & p.adj > 0.001, '**',
                             ifelse(p.adj <= 0.001, '***', '****')))))


## add this to the plot 
library(ggpubr)
p +  
    stat_pvalue_manual(
        stat.test, bracket.nudge.y = 1, hide.ns = TRUE,
        label = "{p.adj.signif}",size=6
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

# but we found one is very close , we adjust mannually
stat.test$y.position[27] <- 14
p <-p +  
    stat_pvalue_manual(
        stat.test, bracket.nudge.y = 1, hide.ns = TRUE,
        label = "{p.adj.signif}",size=6
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(filename = "Fig4_smFISH_IF_Violin_plot.pdf", plot = p, width = 6, height = 6)