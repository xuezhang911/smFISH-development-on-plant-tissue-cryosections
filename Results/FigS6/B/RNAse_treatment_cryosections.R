setwd("/Volumes/Expansion/zhang_et_al_2024_finaldata/Results/FigS6/B")
library(ggplot2)

######################
####################### NRT1.9
##################################


data<-read.delim("__FQ_batch_summary_MATURE_230203_CY5.txt")
head(data)

data_control<-data[data$FILE=="sr_Rroot4_TRANSVERSE_con2l__outline.txt",]
COND<-rep("RNAse(-)",dim(data_control)[1])
data_control<-cbind(COND,data_control)
head(data_control)

data_treat<-data[data$FILE=="sr_Rroot4_TRANSVERSE_Rnase1-2__outline.txt",]
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
p<-p+ geom_boxplot(coef=0.95,width=0.2, outlier.shape = NA, fatten=2, notch= FALSE, notchwidth = 2, fill="white")+
    #stat_summary(fun.y=mean,geom="errorbar", aes(ymax = ..y.., ymin = ..y..), size=0.5, width=0.1, linetype = "dashed")+
    theme(panel.background =element_rect(fill = "white", colour = "black", size = 1),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text.x=element_text(hjust=0.5),
          plot.title = element_text(size = 25,hjust=0.5),
          #panel.border = element_blank()
          panel.grid = element_blank(),
          legend.position="none"
    )+
    labs(title = expression(bolditalic("NRT1.9")),
         x = element_blank(), 
         y = "Transcripts per cell")+
    expand_limits(y = 0)
p


a1<-aov(data_join$N_thres_Total ~ data_log$COND)
summary(a1)
test<-pairwise.t.test(data_join$N_thres_Total,data_join$COND)
p_val = round(test$p.value,digits=4)
p <- p+geom_signif(annotations="****",y_position = 20, xmin=1.2,xmax=2,tip_length = 0.01)

ggsave("violin_cross_NRT1.9_RNAse_treatment.pdf",plot = p,width = 4,height = 4)
