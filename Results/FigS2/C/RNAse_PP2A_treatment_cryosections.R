
setwd("/Volumes/Expansion/zhang_et_al_2024_finaldata/Results/FigS2/C")
library(ggplot2)
library(ggsignif)
######################
####################### PP2A
##################################


data<-read.delim("__FQ_batch_summary_MATURE_230426.txt")
head(data)

data_control<-data[data$FILE=="ROOT12_pp2a_GAPDH_5CON_1-1ch1_aligned_outline.txt",]
COND<-rep("RNAse(-)",dim(data_control)[1])
data_control<-cbind(COND,data_control)
head(data_control)

data_treat<-data[data$FILE=="ROOT12_pp2a_GAPDH_5CON_RNASE-1_aligned_ch1_outline.txt",]
COND<-rep("RNAse(+)",dim(data_treat)[1])
data_treat<-cbind(COND,data_treat)
head(data_treat)

data_join<-rbind(data_control,data_treat)
head(data_join)

log_counts<-log(data_join[,7]+1)
tx_area<-(data_join[,7]+1)/data_join[,4]
log_tx_area<-log(tx_area)
data_log<-cbind(data_join,log_counts,tx_area,log_tx_area)
head(data_log)

#violin plots

p <- ggplot(data, aes(x=data_join$COND, y=data_join$N_thres_Total, fill=data_join$COND)) + 
    geom_violin(trim = TRUE, size=0.5)
p<-p+scale_fill_manual(values=c("gray40", "gray90"))
p<-p+ geom_boxplot(coef=0.95,width=0.3, outlier.shape = NA, fatten=2, notch= FALSE, notchwidth = 2, fill="white")+
    #stat_summary(fun.y=mean,geom="errorbar", aes(ymax = ..y.., ymin = ..y..), size=0.5, width=0.1, linetype = "dashed")+
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 1),
          axis.text = element_text(size = 20),
          axis.title = element_text(colour="black",size = 20),
          axis.text.x=element_text(hjust=0.5),
          plot.title = element_text(colour = "black", size = 25, face = "bold", hjust = 0.5, vjust = 0.5),
          #panel.border = element_blank()
          panel.grid = element_blank(),
          legend.position="none"
    )+
    labs(
         x = element_blank(),  title = expression(italic(PP2A)),
         y = "Transcripts per cell")+
    expand_limits(y = 0)
p


test<-pairwise.t.test(data_join$N_thres_Total,data_join$COND)
p_val = round(test$p.value,digits=4)
# add statistics 
 m1 <- p+geom_signif(annotations = "****",y_position=10,xmin=1.1,xmax=2,tip_length = 0.01,textsize=6)
ggsave("violin_cross_PP2A_RNAse_treatmentV.pdf",plot = m1,width = 4,height = 4)



