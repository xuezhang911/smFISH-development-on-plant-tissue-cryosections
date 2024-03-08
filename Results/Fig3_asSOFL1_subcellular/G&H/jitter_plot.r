setwd("/Volumes/Expansion/zhang_et_al_2024_finaldata/Fig3_asSOFL1_subcellular/G&H")
library(ggplot2)
library(smatr)
library(ggpubr)
library(reshape2)
library(xlsx)
analyze_gene_expression <- function(var) {
  # List the files in the directory
    files <- list.files(pattern = paste0(var, "\\.txt$"))
    # Read each file and store its contents in the list
    data_list <- list()
    for (file in files) {
        data <- read.table(file, header = TRUE)
        data_list[[file]] <- data
    }
    # combine into a single data frame
    data <- do.call(rbind, data_list)
    #filter cells without count
    data <- data[!data$N_thres_Total == 0,]
    data$cyto <- data$N_thres_Total - data$N_thres_Nuc
    #select data for cyto,nuc,toto
    data <- data[, c((ncol(data)-2):ncol(data))]
    # rename the data
    colnames(data) <- c('tot', 'nuc', 'cyto')
    # R loop for count as numeric
    for (x in c('cyto', 'nuc', 'tot')) {
        data[,names(data) == x] <- as.numeric(data[,names(data) == x])
    }
    # calculate the count ratio compare to total count in individual cell
    for (x in c('cyto', 'nuc')) {
        data[,names(data) == x] <- data[,names(data) == x] / data[,names(data) == 'tot'] * 100
    }
    
  # keep only nuc, cyto and transform data 
    m <- melt(data[,!names(data) == 'tot'])
    
    # plotting
    p <- ggplot(m, aes(x = variable, y = value)) + 
        geom_jitter(aes(color = variable), position = position_jitter(0.2), size = 2, alpha = 0.5) +
        scale_color_manual(values = c('#EA5C15', '#002D8E')) +
        labs(x = "", y = "Transcript count%", title = bquote(bolditalic(.(var)))) +
        stat_summary(fun.y = median, geom = "crossbar", size = 0.5, width = 0.5, color = "black") +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
        theme(panel.grid.major = element_blank(), axis.line = element_line(colour = "black"), 
              panel.background = element_blank(),
              axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size = 20, angle = 35, hjust = 1, vjust = 1),
              plot.title = element_text(colour = "black", size = 25, face = "bold", hjust = 0.5, vjust = 0.5),
              legend.text = element_text(size = 16),
              legend.title = element_blank(),
              legend.position = '') 
              #statistics
a1<-aov(m$value~ m$variable)
summary(a1)
test<-TukeyHSD(a1)
test
# add significance in
m <- as.data.frame(test$`m$variable`)
m$Significance <- ifelse(m$`p adj` > 0.05, "ns", 
                      ifelse(m$`p adj` <= 0.05 & m$`p adj` > 0.01, "*", 
                      ifelse(m$`p adj` <= 0.01 & m$`p adj` > 0.001, "**", 
                      ifelse(m$`p adj` <= 0.001 & m$`p adj` > 0.0001, "***", "***"))))
                      print(m)
#save
write.xlsx(m,paste0('/Volumes/Expansion/zhang_et_al_2024_finaldata/Fig3_asSOFL1_subcellular/G&H/',var,'smFISH_statistics.xlsx'))
   # save data 
    ggsave(filename = paste0("Fig3_", var, "_sub.pdf"), plot = p, width = 6, height = 6)
}


# Run the function for both asSOFL1 and PP2A
analyze_gene_expression("asSOFL1")
analyze_gene_expression("PP2A")
