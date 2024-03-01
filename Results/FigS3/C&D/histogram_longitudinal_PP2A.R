setwd("/Volumes/Expansion/zhang_et_al_2024_finaldata/Results/FigS3/C&D")

library(ggplot2)

#################
#################


data <- read.delim("__FQ_batch_summary_MATURE_230424_long_leaf_PP2A.txt")
head(data)

#histograms

p<- ggplot(data,aes(x=N_thres_Total))+geom_histogram(alpha=0.4, color="black", fill="gray", bins = 10)

p<-p+geom_vline(aes(xintercept = median(N_thres_Total)),col="black",size=1,linetype="dashed",)+
    xlab("Number of transcripts per cell") + ylab("Frequency")



p<-p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             legend.title = element_blank(),
             plot.title = element_text(colour="black", size =20,face="bold.italic",hjust = 0.5 ),
             legend.text = element_text(size = 20),
             legend.position = "",
             axis.text=element_text(size=20),
             axis.title=element_text(size=20))+ggtitle(expression(bolditalic("PP2A") ~ bold("young leaf longitudinal-section")))

p
ggsave("leaf_long_PP2A.pdf",plot = p,width = 6,height = 6)

 
 for (var in c('leaf','root')) {
       files <- list.files(pattern = paste0(var, "\\.txt$"))
       data <- read.delim(files)
       
p<- ggplot(data,aes(x=N_thres_Total))+geom_histogram(alpha=0.4, color="black", fill="gray", bins = 10)

p<-p+geom_vline(aes(xintercept = median(N_thres_Total)),col="black",size=1,linetype="dashed",)+
    xlab("Number of transcripts per cell") + ylab("Frequency")

plot_title <- paste0(var," ","longitudinal-section")
p<-p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             legend.title = element_blank(),
             plot.title = element_text(colour="black", size =25,hjust = 0.48 ),
             legend.text = element_text(size = 16),
             legend.position = "",
             axis.text=element_text(size=20),
             axis.title=element_text(size=20))+ labs(title = plot_title)

p
ggsave(paste0("Arabidopsis_",var,"_longitudinal_PP2A.pdf"),plot=p,width=7,height=4)

    }




for (var in c('leaf','root')) {
    files <- list.files(pattern = paste0(var, "\\.txt$"))
    for (file in files) {
        data <- read.delim(file)
      p<- ggplot(data,aes(x=N_thres_Total))+geom_histogram(alpha=0.4, color="black", fill="gray", bins = 10)

p<-p+geom_vline(aes(xintercept = median(N_thres_Total)),col="black",size=1,linetype="dashed",)+
    xlab("Number of transcripts per cell") + ylab("Frequency")

plot_title <- paste0(var," ","longitudinal-section")
p<-p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             legend.title = element_blank(),
             plot.title = element_text(colour="black", size =25,hjust = 0.48 ),
             legend.text = element_text(size = 16),
             legend.position = "",
             axis.text=element_text(size=20),
             axis.title=element_text(size=20))+ labs(title = plot_title)

p
ggsave(paste0("Arabidopsis_",var,"_longitudinal_PP2A.pdf"),plot=p,width=6,height=6)
    }
}






