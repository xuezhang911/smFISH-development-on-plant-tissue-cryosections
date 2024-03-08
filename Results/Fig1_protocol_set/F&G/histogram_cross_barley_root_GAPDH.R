setwd("/Volumes/Expansion/zhang_et_al_2024_finaldata/Fig1_protocol_set/F&G")
library(ggplot2)
#################
#transverse
#################


data <- read.delim("__FQ_batch_summary_MATURE_230504_root6.txt")
head(data)
data <- data[order(data$N_thres_Total),]

#histograms

p<- ggplot(data,aes(x=N_thres_Total))+geom_histogram(alpha=0.4, color="black", fill="gray", bins = 10,boundary=0)

p<-p+geom_vline(aes(xintercept = median(N_thres_Total)),col="black",size=1,linetype="dashed",)+
    xlab("Number of transcripts per cell") + ylab("Frequency")


p<-p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             legend.title = element_blank(),
             plot.title = element_text(colour="black", size =20,face="bold.italic",hjust = 0 ),
             legend.text =  element_blank(),
             legend.position = "",
             axis.text=element_text(size=16),
             axis.title=element_text(size=16))

p

ggsave("Fig1_histogram_Barley_root_trans_HvGAPDH.pdf",p,width = 5, height = 3) 



