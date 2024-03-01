setwd('/Volumes/Expansion/zhang_et_al_2024_finaldata/Results/FigS4/B')
library(ggplot2)
 
 for (var in c('leaf','root')) {
       files <- list.files(pattern = paste0(var, "\\.txt$"))
       data <- read.delim(files)
       
p<- ggplot(data,aes(x=N_thres_Total))+geom_histogram(alpha=0.4, color="black", fill="gray", bins = 8)

p<-p+geom_vline(aes(xintercept = median(N_thres_Total)),col="black",size=1,linetype="dashed",)+
    xlab("Number of transcripts per cell") + ylab("Frequency")

plot_title <- paste0(var," ","cross-section")
p<-p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             legend.title = element_blank(),
             plot.title = element_text(colour="black", size =25,hjust = 0.48 ),
             legend.text = element_text(size = 16),
             legend.position = "",
             axis.text=element_text(size=20),
             axis.title=element_text(size=20))+ labs(title = plot_title)

p
ggsave(paste0("barley_",var,"_cross_HvGAPDH.pdf"),plot=p,width=7,height=4)

    }

 
 