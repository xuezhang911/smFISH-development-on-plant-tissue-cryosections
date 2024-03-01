setwd("/Volumes/Expansion/zhang_et_al_2024_finaldata/Fig3_asSOFL1_subcellular/B")
library(ggplot2)
 data <- read.delim("__FQ_batch_summary_MATURE_230503.txt",header=F,sep="\t")
 names(data) <- data[5,]
data <- data[-c(1:5),]
data$FILE <- as.factor(data$FILE)
levels(data$FILE)
data <- data[data$FILE=="root5_con_sr_asSOFL1_2_aligned_ch1_outline.txt",]
 transcript_counts <- table(data$N_thres_Total)
 data$N_thres_Total <- as.numeric(data$N_thres_Total)
 transcript_counts_df <- data.frame(transcripts = as.numeric(names(transcript_counts)),
                count = as.numeric(transcript_counts))
  p <- ggplot(transcript_counts_df, aes(x = transcripts, y = count)) +geom_vline(aes(xintercept = median(data$N_thres_Total)),col="black",size=1,linetype="dashed",)+ geom_bar(stat = "identity", fill = "gray", color = "black",alpha=0.4,width = 1)+  
    theme(
        panel.grid.major = element_blank(),axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 20, face = "bold.italic", hjust = 0.5,vjust=0.5),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)
    )+ xlab("Number of transcripts per cell") +
    ylab("number of cells") +
    ggtitle(expression(bolditalic("asSOFL1")))

p
     

              ggsave(filename = "Fig3_smFISH_bar_plot_asSOFL1.pdf", plot = p, width = 5, height = 6)