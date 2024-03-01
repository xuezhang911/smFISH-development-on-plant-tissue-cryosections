## barplot for TMP in RNA-seq
setwd('/Volumes/Expansion/zhang_et_al_2024_finaldata/Results/FigS8/A')
library(ggplot2)
library(reshape2)
library(readxl)
library(dplyr)
library(ggpubr)
library(rstatix)
# LOAD PP2A
data <- read.delim("__FQ_batch_summary_MATURE_230421_PP2A_new.txt")
data <- data.frame(data$FILE,data$N_thres_Total)
names(data) <- c('rep','PP2A')
PP2A<-data
PP2A$rep <- as.factor(PP2A$rep)
# need 3nd replicates
data <- read.delim("__FQ_batch_summary_MATURE_230423_pp2a.txt")
data <- data[data$FILE=="root13_CON_nrt1.9_PP2A_5-2_ch3_outline.txt",]
data <- data.frame(data$FILE,data$N_thres_Total)
names(data) <- c('rep','PP2A')
PP2A <- rbind(PP2A,data)
PP2A$rep <- as.factor(PP2A$rep)

PP2A <- PP2A %>%
    group_by(rep) %>%
    summarize(PP2A = sum(PP2A))
PP2A <- melt(PP2A)

###################################################################### for asSOFL1
data <- read.delim("__FQ_batch_summary_MATURE_230503_asSOFL1_strict.txt")
data <- data.frame(data$FILE,data$N_thres_Total)
names(data) <- c('rep','asSOFL1')
asSOFL1<-data

# need 3nd replicates
data <- read_excel("counts_per_celltype_asSOFL1.xlsx")
head(data)
# modify the structure of data 
names(data) <- c("replicates",'celltype',"transcripts")
data$replicates <- as.factor(data$replicates)
levels(data$replicates)
data <- data[data$replicates=="root5_con_sr_asSOFL1-2",] 
# here, we no need cell type information # subset data 
data <- as.data.frame(data[,c(1,3)])
names(data) <- c('rep','asSOFL1')
asSOFL1 <- rbind(asSOFL1,data)
asSOFL1$rep <- as.factor(asSOFL1$rep)
asSOFL1 <- asSOFL1 %>%
    group_by(rep) %>%
    summarize(asSOFL1 = sum(asSOFL1))
asSOFL1 <- melt(asSOFL1)
###################################################################### for NRT1.9
data <- read_excel("cyto_nuclei_NRT1.9_update.xlsx") 
# only analysis 5-2
data <- data[data$FILE=="root13_CON_nrt1.9_PP2A_5-2_ch3_outline.txt",]
data <- data.frame(data$FILE,data$N_thres_Total)
names(data) <- c('rep','NRT1.9')
NRT1.9<-data
# rep3
data <- read_excel("development_COUNTS_NRT1.9_celltype_update.xlsx") 
 data$FILE <- as.factor(data$FILE)
 levels(data$FILE) 
 ## I decided to mix two development stage 
 # modify the structure of data 
names(data) <- c("replicates",'celltype',"transcripts")
data$replicates <- as.factor(data$replicates)
levels(data$replicates)
# subset data 
 c <- c('root13_CON_SRnrt1.9_PP2A_4',"sr_Rroot4_TRANSVERSE_con2l")
 data <- data[data$replicates%in%c,] 
data <- as.data.frame(data[,c(1,3)])
names(data) <- c('rep','NRT1.9')
NRT1.9 <- rbind(NRT1.9,data)
NRT1.9$rep <- as.factor(NRT1.9$rep)
NRT1.9 <- NRT1.9 %>%
    group_by(rep) %>%
    summarize(NRT1.9 = sum(NRT1.9))
NRT1.9 <- melt(NRT1.9)
################################################################## combine data 
df <- rbind(NRT1.9,PP2A,asSOFL1)
## no need replicates more 
df <- df[,c(2:3)]
# check the structure 
 str(df)

  # calculate summary statistics 
df_summary <- df %>% # Specifying the name of the new dataframe `df_summary`
  group_by(variable) %>% # we aim to know the value of two genes
  summarise(            # Function for mutating data into summary stats.
    mean_relative_counts = mean(value), # Creating column to store means
    N = n(),                              # Creating column to store sample size
    SE_total_relative_counts= sd(value) / sqrt(n()) # Creating column to store standard error
  )

# ggplot2 draw figures 
df_summary$variable <- factor(df_summary$variable,levels=c('asSOFL1','PP2A','NRT1.9'))
View(df_summary)
# creating barplot with ggplot


p <-   ggplot(df_summary, aes(x = variable,
                              y = mean_relative_counts)) + # Specifies the data we want plotted
  geom_bar(stat="identity",position=position_dodge()) + # Tells R what type of plot we want
  geom_errorbar(aes(
    ymin = mean_relative_counts- SE_total_relative_counts, # Established error bars
    ymax = mean_relative_counts + SE_total_relative_counts,
    width = 0.3) ) +labs(x = element_blank(), 
            y = "Transcript abundance")+ scale_y_continuous(expand = c(0,0))
                                                               
p<-p+ggtitle("cryo-smFISH") +theme_classic(base_size = 16)+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 20, face = "italic"),
          axis.text.y =element_text(size = 20),
          # if we have lengend legend.title = element_blank(),
          axis.title = element_text(colour="black",size = 20),
          plot.title = element_text(colour="black", size =20,face="bold",hjust = 0.5 ))
names(df) [1]<- 'gene'

#statistics
          stat.test <- df %>%
    t_test(value ~ gene)
 # add xy position 
stat.test <- stat.test %>% add_xy_position(x = "gene",  step.increase = 0.12,)

 stat.test <- stat.test %>%
    mutate(p.adj.signif = ifelse(p.adj > 0.05, 'ns',
                                 ifelse(p.adj <= 0.05 & p.adj > 0.01, '*',
                                        ifelse(p.adj <= 0.01 & p.adj > 0.001, '**',
                                               ifelse(p.adj <= 0.001, '***', '****')))))

        # save data 
        m <- stat.test[,c(2,3,9,10)]
        write.xlsx(m,'/Volumes/Expansion/zhang_et_al_2024_finaldata/Results/FigS8/A/cryo_smFISH_statistics.xlsx')
    ## add this to the plot 
   
p<-p +  
    stat_pvalue_manual(
        stat.test, bracket.nudge.y = 1, hide.ns = F,
        label = "{p.adj.signif}",size=6
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

      ggsave("total_transcripts_smFISH.pdf",p,width = 5, height = 5)