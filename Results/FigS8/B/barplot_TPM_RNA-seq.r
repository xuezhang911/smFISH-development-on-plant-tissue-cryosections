## barplot for TMP in RNA-seq
setwd('/Volumes/Expansion/zhang_et_al_2024_finaldata/Results/FigS8/B')
library(ggplot2)
library(reshape2)
library(readxl)
library(dplyr)
library(ggpubr)
library(rstatix)
# create data_frame or loading data 
df <- data.frame(gene=rep(c("NRT1.9","PP2A","asSOFL1"),each=2),
                 total_relative_counts=c(46,47.27,56.89,61.12,0.32,0.37))

# calculate summary statistics 
df_summary <- df %>% # Specifying the name of the new dataframe `df_summary`
  group_by(gene) %>% # we aim to know the value of two genes
  summarise(            # Function for mutating data into summary stats.
    mean_relative_counts = mean(total_relative_counts), # Creating column to store means
    N = n(),                              # Creating column to store sample size
    SE_total_relative_counts= sd(total_relative_counts) / sqrt(n()) # Creating column to store standard error
  )


df_summary$gene <- factor(df_summary$gene,levels=c('asSOFL1','PP2A','NRT1.9'))
View(df_summary)
# creating barplot with ggplot


p <-   ggplot(df_summary, aes(x = gene,
                              y = mean_relative_counts)) + # Specifies the data we want plotted
  geom_bar(stat="identity",position=position_dodge()) + # Tells R what type of plot we want
  geom_errorbar(aes(
    ymin = mean_relative_counts- SE_total_relative_counts, # Established error bars
    ymax = mean_relative_counts + SE_total_relative_counts,
    width = 0.3) ) +labs(x = element_blank(), 
            y = "TPM counts")+ scale_y_continuous(expand = c(0,0), limits = c(0, 70))
                                                               
p<-p+ggtitle("bulk RNA-seq") +theme_classic(base_size = 16)+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 20, face = "italic"),
          axis.text.y =element_text(size = 20),
          # if we have lengend legend.title = element_blank(),
          axis.title = element_text(colour="black",size = 20),
          plot.title = element_text(colour="black", size =20,face="bold",hjust = 0.5 ))
          library(dplyr)
          library(rstatix)

          stat.test <- df %>%
    t_test(total_relative_counts ~ gene)

    stat.test <- stat.test %>%
    mutate(p.adj.signif = ifelse(p.adj > 0.05, 'ns',
                                 ifelse(p.adj <= 0.05 & p.adj > 0.01, '*',
                                        ifelse(p.adj <= 0.01 & p.adj > 0.001, '**',
                                               ifelse(p.adj <= 0.001, '***', '****')))))
 # add xy position 
stat.test <- stat.test %>% add_xy_position(x = "gene",  step.increase = 0.12,)
    ## add this to the plot 
  
p<-p +  
    stat_pvalue_manual(
        stat.test, bracket.nudge.y = 1, hide.ns = T,
        label = "{p.adj.signif}",size=6
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
        # save data 
        m <- stat.test[,c(2,3,9,10)]
        write.xlsx(m,'/Volumes/Expansion/zhang_et_al_2024_finaldata/Results/FigS8/B/bulk_RNA-seq_statistics.xlsx')

ggsave("total_transcripts_RNA-seq.pdf",p,width = 5, height = 5)

