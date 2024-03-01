setwd('/Volumes/Expansion/zhang_et_al_2024_finaldata/Results/FigS8/D')
library(readxl)
library(reshape2)
library(xlsx)
library(dplyr)
library(rstatix)
################################################################ for PP2A file1: 
 # load data 
data <- read.delim("__FQ_batch_summary_MATURE_230421_PP2A_new.txt")

# only analysis 5-2
data <- data[data$FILE=="root11_CON_PP2A_asSOFL1_2_ch2_scaled_outline.txt",] 

# for cell, we only keep the number 
data$CELL <- as.numeric(substring(data$CELL, 9, nchar(data$CELL)) )

# for 1-2

cell_numbers <- c(
    80,52,46,44,30,18,23,31,35,
    47,40,32,21,79,
59,62,48,28,13,8,6,10,83,82,
                58,64,57,19,3,4,26,75,
               1,2,5,12,25,45,63,65,67,69,68,66,71,53,70,20,7,72,73,74,76,
56,54,50,42,41,36,33,29,15,17,27,34,37,77,43,
               60,51,39,22,16,9,78,11,14,24,38,49,55,61,81

)


# check if it has duplicates and sort it to see if anthing is missing 
which(duplicated(cell_numbers))
sort(cell_numbers)
cell_types <- c(rep('Phloem',9),rep('Xylem',5),rep('Endodermis',10),rep('Cortex',8),rep('Epidermis',21),rep('Procambium',15),rep('Pericycle',15))

# length 
 length(cell_types)


data$CELL%in%cell_numbers


data$celltype <- NA
# Loop through the cell numbers
for (i in 1:length(cell_numbers)) {
    # Find rows where data$CELL matches the current cell number
    matching_rows <- data$CELL == cell_numbers[i]
    
    # Assign the corresponding cell type to the matched rows
    data$celltype[matching_rows] <- cell_types[i]
}
# merge to data 
data$AREA_cyto <- data$AREA_cell-data$AREA_nuc
# DATA normalized by the size of area
data$N_thres_Cyto <- data$N_thres_Total-data$N_thres_Nuc
data$N_Area_tot <- data$N_thres_Total/data$AREA_cell
data$N_Area_nuc <- data$N_thres_Nuc/data$AREA_nuc
data$N_Area_cyto <- data$N_thres_Cyto/data$AREA_cyto

# subset data 
c <- c("FILE",'CELL',"N_Area_tot","N_Area_nuc","N_Area_cyto","N_thres_Total","N_thres_Nuc","N_thres_Cyto",'celltype')
data
data <- data[c]

# create dataframe 
m <- cbind(data$celltype,data$N_thres_Total,data$N_thres_Nuc,data$N_thres_Cyto) 


m <- as.data.frame(m)
names(m) <- c('celltype','total','nuc','cyto')
# change the structure
s <- melt(m, id.vars = "celltype", 
          measure.vars = c("total", "nuc",'cyto'),
          variable.name = "compartment", 
          value.name = "Value")

str(s)
s$Value <- as.numeric(s$Value)
s$celltype <- as.factor(s$celltype)

s$celltype <- factor(s$celltype,levels=c("Epidermis", "Cortex", "Endodermis", "Pericycle", "Phloem", "Procambium", "Xylem"))
# or normalized by the cell quantity optional 
cell_count <- s %>%
    group_by(celltype, compartment) %>%
    summarize(Cell_Count = n())

s<- s %>%
      left_join(cell_count, by =c("celltype","compartment") )
s$norm <- (s$Value/s$Cell_Count)*10
# create data with normalized value 
s <- s[,c(1,5)]
names(s)
names(s)[2] <- 'PP2A'
# rename as df
df <- s # this one so far only include 1-2, we will add 1-3, and 2-4 for PP2A, three rep for asOSFL1 , three rep for NRT1.9
#################### ################################################rep2 for PP2A

data <- read.delim("__FQ_batch_summary_MATURE_230421_PP2A_new.txt")
data <- data[data$FILE=="root11_CON_PP2A_asSOFL1_3_1_ch2_scaled_outline.txt",]
# for cell, we only keep the number 
data$CELL <- as.numeric(substring(data$CELL, 9, nchar(data$CELL)) )
#for 1-3-1 PP2A_SOFL1
cell_numbers <- c(
    27,26,23,51,53,44,
    49,40,34,28,24,79,
31,52,57,15,10,9,12,19,72,59,64,65,
               63,55,35,11,7,73,18,78,47,66,
              75,58,76,70,68,69,67,60,77,13,8,5,2,3,1,4,6,74,22,45,
20,29,33,36,38,39,41,42,46,50,71,80,
              62,61,54,43,30,25,17,14,16,21,32,37,48,56

)


# check if it has duplicates and sort it to see if anthing is missing 
which(duplicated(cell_numbers))
sort(cell_numbers)
cell_types <- c(rep('Phloem',6),rep('Xylem',6),rep('Endodermis',12),rep('Cortex',10),rep('Epidermis',20),rep('Procambium',12),rep('Pericycle',14))
# length 
 length(cell_types)
data$CELL%in%cell_numbers


data$celltype <- NA
# Loop through the cell numbers
for (i in 1:length(cell_numbers)) {
    # Find rows where data$CELL matches the current cell number
    matching_rows <- data$CELL == cell_numbers[i]
    
    # Assign the corresponding cell type to the matched rows
    data$celltype[matching_rows] <- cell_types[i]
}
# merge to data 
data$AREA_cyto <- data$AREA_cell-data$AREA_nuc
# DATA normalized by the size of area
data$N_thres_Cyto <- data$N_thres_Total-data$N_thres_Nuc
data$N_Area_tot <- data$N_thres_Total/data$AREA_cell
data$N_Area_nuc <- data$N_thres_Nuc/data$AREA_nuc
data$N_Area_cyto <- data$N_thres_Cyto/data$AREA_cyto

# subset data 
c <- c("FILE",'CELL',"N_Area_tot","N_Area_nuc","N_Area_cyto","N_thres_Total","N_thres_Nuc","N_thres_Cyto",'celltype')
data
data <- data[c]

# create dataframe 
m <- cbind(data$celltype,data$N_thres_Total,data$N_thres_Nuc,data$N_thres_Cyto) 


m <- as.data.frame(m)
names(m) <- c('celltype','total','nuc','cyto')
# change the structure
s <- melt(m, id.vars = "celltype", 
          measure.vars = c("total", "nuc",'cyto'),
          variable.name = "compartment", 
          value.name = "Value")

str(s)

s$Value <- as.numeric(s$Value)
s$celltype <- as.factor(s$celltype)

s$celltype <- factor(s$celltype,levels=c("Epidermis", "Cortex", "Endodermis", "Pericycle", "Phloem", "Procambium", "Xylem"))
# or normalized by the cell quantity optional 
cell_count <- s %>%
    group_by(celltype, compartment) %>%
    summarize(Cell_Count = n())

s<- s %>%
      left_join(cell_count, by =c("celltype","compartment") )
s$norm <- (s$Value/s$Cell_Count)*10
# create data with normalized value 
s <- s[,c(1,5)]
names(s)
names(s)[2] <- 'PP2A'
# combine data 
df <- rbind(df,s) #here we have two replicates 
############################################################################################## rep3
data <- read.delim("__FQ_batch_summary_MATURE_230423_pp2a.txt")
data <- data[data$FILE=="root13_CON_nrt1.9_PP2A_5-2_ch3_outline.txt",]
# for cell, we only keep the number 
data$CELL <- as.numeric(substring(data$CELL, 9, nchar(data$CELL)) )

cell_numbers <- c(35,34,28,29,31,43,20,
                  21,26,39,53,58,
                  56,38,19,8,10,14,42,61,67,65,
                  69,71,73,66,47,37,11,6,5,7,13,24,48,
                  68,72,76,74,75,70,60,45,36,22,77,78,3,1,2,4,9,79,30,52,63,
 44,46,51,18,15,27,40,49,55,
12,16,17,23,25,32,33,41,50,54,57,59,62,64,80)

cell_types <- c(rep('Phloem',7),rep('Xylem',5),rep('Endodermis',10),rep('Cortex',13),rep('Epidermis',21),rep('Procambium',9),rep('Pericycle',15))

# check if data$CELL is equal to cell_numbers
data$CELL%in%cell_numbers

data$celltype <- NA
# Loop through the cell numbers
for (i in 1:length(cell_numbers)) {
    # Find rows where data$CELL matches the current cell number
    matching_rows <- data$CELL == cell_numbers[i]
    
    # Assign the corresponding cell type to the matched rows
    data$celltype[matching_rows] <- cell_types[i]
}
data$celltype

data<-data[data$AREA_nuc>0,]
data <- as.data.frame(data)
data$AREA_cell <- as.numeric(data$AREA_cell)
data$AREA_nuc <- as.numeric(data$AREA_nuc)


# merge to data 
data$AREA_cyto <- data$AREA_cell-data$AREA_nuc
# DATA normalized by the size of area
data$N_thres_Cyto <- data$N_thres_Total-data$N_thres_Nuc
data$N_Area_tot <- data$N_thres_Total/data$AREA_cell
data$N_Area_nuc <- data$N_thres_Nuc/data$AREA_nuc
data$N_Area_cyto <- data$N_thres_Cyto/data$AREA_cyto

# subset data 
c <- c("FILE",'CELL',"N_Area_tot","N_Area_nuc","N_Area_cyto","N_thres_Total","N_thres_Nuc","N_thres_Cyto",'celltype')
data
data <- data[c]

# create dataframe 
m <- cbind(data$celltype,data$N_thres_Total,data$N_thres_Nuc,data$N_thres_Cyto) 


m <- as.data.frame(m)
names(m) <- c('celltype','total','nuc','cyto')
# change the structure
s <- melt(m, id.vars = "celltype", 
          measure.vars = c("total", "nuc",'cyto'),
          variable.name = "compartment", 
          value.name = "Value")

str(s)
s$Value <- as.numeric(s$Value)
s$celltype <- as.factor(s$celltype)

s$celltype <- factor(s$celltype,levels=c("Epidermis", "Cortex", "Endodermis", "Pericycle", "Phloem", "Procambium", "Xylem"))
# or normalized by the cell quantity optional 
cell_count <- s %>%
    group_by(celltype, compartment) %>%
    summarize(Cell_Count = n())

s<- s %>%
      left_join(cell_count, by =c("celltype","compartment") )
s$norm <- (s$Value/s$Cell_Count)*10
# create data with normalized value 
s <- s[,c(1,5)]
names(s)
names(s)[2] <- 'PP2A'
# combine data 
df <- rbind(df,s) #here we have three replicates 
##################################################################################################### for asSOFL1 rep1 
 data <- read.delim("__FQ_batch_summary_MATURE_230503_asSOFL1_strict.txt")

data <- data[data$FILE=="root11_CON_PP2A_asSOFL1_2_ch1_scaled_outline.txt",] 

# for cell, we only keep the number 
data$CELL <- as.numeric(substring(data$CELL, 9, nchar(data$CELL)) )

# for 1-2

cell_numbers <- c(
    80,52,46,44,30,18,23,31,35,
    47,40,32,21,79,
59,62,48,28,13,8,6,10,83,82,
                58,64,57,19,3,4,26,75,
               1,2,5,12,25,45,63,65,67,69,68,66,71,53,70,20,7,72,73,74,76,
56,54,50,42,41,36,33,29,15,17,27,34,37,77,43,
               60,51,39,22,16,9,78,11,14,24,38,49,55,61,81

)


# check if it has duplicates and sort it to see if anthing is missing 
which(duplicated(cell_numbers))
sort(cell_numbers)
cell_types <- c(rep('Phloem',9),rep('Xylem',5),rep('Endodermis',10),rep('Cortex',8),rep('Epidermis',21),rep('Procambium',15),rep('Pericycle',15))

# length 
 length(cell_types)


data$CELL%in%cell_numbers
data$celltype <- NA
# Loop through the cell numbers
for (i in 1:length(cell_numbers)) {
    # Find rows where data$CELL matches the current cell number
    matching_rows <- data$CELL == cell_numbers[i]
    
    # Assign the corresponding cell type to the matched rows
    data$celltype[matching_rows] <- cell_types[i]
}
# merge to data 
data$AREA_cyto <- data$AREA_cell-data$AREA_nuc
# DATA normalized by the size of area
data$N_thres_Cyto <- data$N_thres_Total-data$N_thres_Nuc
data$N_Area_tot <- data$N_thres_Total/data$AREA_cell
data$N_Area_nuc <- data$N_thres_Nuc/data$AREA_nuc
data$N_Area_cyto <- data$N_thres_Cyto/data$AREA_cyto

# subset data 
c <- c("FILE",'CELL',"N_Area_tot","N_Area_nuc","N_Area_cyto","N_thres_Total","N_thres_Nuc","N_thres_Cyto",'celltype')
data
data <- data[c]

# create dataframe 
m <- cbind(data$celltype,data$N_thres_Total,data$N_thres_Nuc,data$N_thres_Cyto) 


m <- as.data.frame(m)
names(m) <- c('celltype','total','nuc','cyto')
# change the structure
s <- melt(m, id.vars = "celltype", 
          measure.vars = c("total", "nuc",'cyto'),
          variable.name = "compartment", 
          value.name = "Value")

str(s)
s$Value <- as.numeric(s$Value)
s$celltype <- as.factor(s$celltype)

s$celltype <- factor(s$celltype,levels=c("Epidermis", "Cortex", "Endodermis", "Pericycle", "Phloem", "Procambium", "Xylem"))
# or normalized by the cell quantity optional 
cell_count <- s %>%
    group_by(celltype, compartment) %>%
    summarize(Cell_Count = n())

s<- s %>%
      left_join(cell_count, by =c("celltype","compartment") )
s$norm <- (s$Value/s$Cell_Count)*10
# create data with normalized value 
s <- s[,c(1,5)]
names(s)
names(s)[2] <- 'asSOFL1'

## change the stucture , gene as vairable
s <- melt(s, id.vars = "celltype", variable.name = "gene", value.name = "expression")
# do the same for df 
df<- melt(df, id.vars = "celltype", variable.name = "gene", value.name = "expression")
# combine both together 
df <- rbind(df,s)
##################################################################################################### for asSOFL1 rep2
 data <- read.delim("__FQ_batch_summary_MATURE_230503_asSOFL1_strict.txt")

data <- data[data$FILE=="root11_CON_PP2A_asSOFL1_3_1_ch1_scaled_outline.txt",]

# for cell, we only keep the number 
data$CELL <- as.numeric(substring(data$CELL, 9, nchar(data$CELL)) )

#for 1-3-1 PP2A_SOFL1
cell_numbers <- c(
    27,26,23,51,53,44,
    49,40,34,28,24,79,
31,52,57,15,10,9,12,19,72,59,64,65,
               63,55,35,11,7,73,18,78,47,66,
              75,58,76,70,68,69,67,60,77,13,8,5,2,3,1,4,6,74,22,45,
20,29,33,36,38,39,41,42,46,50,71,80,
              62,61,54,43,30,25,17,14,16,21,32,37,48,56

)


# check if it has duplicates and sort it to see if anthing is missing 
which(duplicated(cell_numbers))
sort(cell_numbers)
cell_types <- c(rep('Phloem',6),rep('Xylem',6),rep('Endodermis',12),rep('Cortex',10),rep('Epidermis',20),rep('Procambium',12),rep('Pericycle',14))

# length 
 length(cell_types)


data$CELL%in%cell_numbers
data$celltype <- NA
# Loop through the cell numbers
for (i in 1:length(cell_numbers)) {
    # Find rows where data$CELL matches the current cell number
    matching_rows <- data$CELL == cell_numbers[i]
    
    # Assign the corresponding cell type to the matched rows
    data$celltype[matching_rows] <- cell_types[i]
}
# merge to data 
data$AREA_cyto <- data$AREA_cell-data$AREA_nuc
# DATA normalized by the size of area
data$N_thres_Cyto <- data$N_thres_Total-data$N_thres_Nuc
data$N_Area_tot <- data$N_thres_Total/data$AREA_cell
data$N_Area_nuc <- data$N_thres_Nuc/data$AREA_nuc
data$N_Area_cyto <- data$N_thres_Cyto/data$AREA_cyto

# subset data 
c <- c("FILE",'CELL',"N_Area_tot","N_Area_nuc","N_Area_cyto","N_thres_Total","N_thres_Nuc","N_thres_Cyto",'celltype')
data
data <- data[c]

# create dataframe 
m <- cbind(data$celltype,data$N_thres_Total,data$N_thres_Nuc,data$N_thres_Cyto) 


m <- as.data.frame(m)
names(m) <- c('celltype','total','nuc','cyto')
# change the structure
s <- melt(m, id.vars = "celltype", 
          measure.vars = c("total", "nuc",'cyto'),
          variable.name = "compartment", 
          value.name = "Value")

str(s)
s$Value <- as.numeric(s$Value)
s$celltype <- as.factor(s$celltype)

s$celltype <- factor(s$celltype,levels=c("Epidermis", "Cortex", "Endodermis", "Pericycle", "Phloem", "Procambium", "Xylem"))
# or normalized by the cell quantity optional 
cell_count <- s %>%
    group_by(celltype, compartment) %>%
    summarize(Cell_Count = n())

s<- s %>%
      left_join(cell_count, by =c("celltype","compartment") )
s$norm <- (s$Value/s$Cell_Count)*10
# create data with normalized value 
s <- s[,c(1,5)]
names(s)
names(s)[2] <- 'asSOFL1'

## change the stucture , gene as vairable
s <- melt(s, id.vars = "celltype", variable.name = "gene", value.name = "expression")
# combine both together 
df <- rbind(df,s) # here we inlcude three PP2A, two asSOFL1

####################################################################### rep3 for asSOFL1
#setwd("/Volumes/files/EvoCell_lab/smFISH in cryosections/Final_data_each_figure/FigureS9_comparation among smFISH,scRNA-seq and bulk RNA-seq/A")
#setwd("Y:/EvoCell_lab/smFISH in cryosections/Fig2_NRT1.9_cell_type/results/corrected_results")
# load packages
library(ggplot2)
library(readxl)
library(dplyr)
# load data 
data <- read_excel("counts_per_celltype_asSOFL1.xlsx")
head(data)
# modify the structure of data 
names(data) <- c("replicates",'celltype',"transcripts")
data$replicates <- as.factor(data$replicates)
levels(data$replicates)
data <- data[data$replicates=="root5_con_sr_asSOFL1-2",] 
# normalize data by cell quantity
celltype_counts <- data %>%
    group_by(celltype, replicates) %>%
    summarize(total_cells = n())

data <- data %>%
    left_join(celltype_counts, by =c("celltype","replicates") )
data$norm <- (data$transcripts/data$total_cells)*10
# subset data 
data <- as.data.frame(data[,c(2,5)])
names(data)[2] <- 'asSOFL1'

# change the structure
data<- melt(data, id.vars = "celltype", variable.name = "gene", value.name = "expression")
# combine both together 
df <- rbind(df,data) # here we inlcude three PP2A, three asSOFL1

############################################################# for NRT1.9
data <- read_excel("cyto_nuclei_NRT1.9_update.xlsx") 

# only analysis 5-2
data <- data[data$FILE=="root13_CON_nrt1.9_PP2A_5-2_ch3_outline.txt",]
# for cell, we only keep the number 
data$CELL <- as.numeric(substring(data$CELL, 9, nchar(data$CELL)) )

cell_numbers <- c(35,34,28,29,31,43,20,
                  21,26,39,53,58,
                  56,38,19,8,10,14,42,61,67,65,
                  69,71,73,66,47,37,11,6,5,7,13,24,48,
                  68,72,76,74,75,70,60,45,36,22,77,78,3,1,2,4,9,79,30,52,63,
 44,46,51,18,15,27,40,49,55,
12,16,17,23,25,32,33,41,50,54,57,59,62,64,80)

cell_types <- c(rep('Phloem',7),rep('Xylem',5),rep('Endodermis',10),rep('Cortex',13),rep('Epidermis',21),rep('Procambium',9),rep('Pericycle',15))

# check if data$CELL is equal to cell_numbers
data$CELL%in%cell_numbers

data$celltype <- NA
# Loop through the cell numbers
for (i in 1:length(cell_numbers)) {
    # Find rows where data$CELL matches the current cell number
    matching_rows <- data$CELL == cell_numbers[i]
    
    # Assign the corresponding cell type to the matched rows
    data$celltype[matching_rows] <- cell_types[i]
}
data$celltype

data<-data[data$AREA_nuc>0,]
data <- as.data.frame(data)
data$AREA_cell <- as.numeric(data$AREA_cell)
data$AREA_nuc <- as.numeric(data$AREA_nuc)

# merge to data 
data$AREA_cyto <- data$AREA_cell-data$AREA_nuc
# DATA normalized by the size of area
data$N_thres_Cyto <- data$N_thres_Total-data$N_thres_Nuc
data$N_Area_tot <- data$N_thres_Total/data$AREA_cell
data$N_Area_nuc <- data$N_thres_Nuc/data$AREA_nuc
data$N_Area_cyto <- data$N_thres_Cyto/data$AREA_cyto

# subset data 
c <- c("FILE",'CELL',"N_Area_tot","N_Area_nuc","N_Area_cyto","N_thres_Total","N_thres_Nuc","N_thres_Cyto",'celltype')
data
data <- data[c]

# create dataframe 
m <- cbind(data$celltype,data$N_thres_Total,data$N_thres_Nuc,data$N_thres_Cyto) 


m <- as.data.frame(m)
names(m) <- c('celltype','total','nuc','cyto')
# change the structure
s <- melt(m, id.vars = "celltype", 
          measure.vars = c("total", "nuc",'cyto'),
          variable.name = "compartment", 
          value.name = "Value")

str(s)
s$Value <- as.numeric(s$Value)
s$celltype <- as.factor(s$celltype)

s$celltype <- factor(s$celltype,levels=c("Epidermis", "Cortex", "Endodermis", "Pericycle", "Phloem", "Procambium", "Xylem"))
# or normalized by the cell quantity optional 
cell_count <- s %>%
    group_by(celltype, compartment) %>%
    summarize(Cell_Count = n())

s<- s %>%
      left_join(cell_count, by =c("celltype","compartment") )
s$norm <- (s$Value/s$Cell_Count)*10
# create data with normalized value 
s <- s[,c(1,5)]
names(s)
names(s)[2] <- 'NRT1.9'

## change the stucture , gene as vairable
s <- melt(s, id.vars = "celltype", variable.name = "gene", value.name = "expression")
# combine both together 
df <- rbind(df,s) # here we inlcude three PP2A, three asSOFL1,one replicate for NRT1.9
########################################################################################## rep2,rep3 for NRT1.9
data <- read_excel("development_COUNTS_NRT1.9_celltype_update.xlsx") 
 data$FILE <- as.factor(data$FILE)
 levels(data$FILE) 
 ## I decided to mix two development stage 
 # modify the structure of data 
names(data) <- c("replicates",'celltype',"transcripts")
data$replicates <- as.factor(data$replicates)
levels(data$replicates)
 c <- c('root13_CON_SRnrt1.9_PP2A_4',"sr_Rroot4_TRANSVERSE_con2l")
 data <- data[data$replicates%in%c,] 

# normalize data by cell quantity
celltype_counts <- data %>%
    group_by(celltype, replicates) %>%
    summarize(total_cells = n())

data <- data %>%
    left_join(celltype_counts, by =c("celltype","replicates") )
data$norm <- (data$transcripts/data$total_cells)*10

# subset data 
data <- as.data.frame(data[,c(2,5)])
names(data)[2] <- 'NRT1.9'

# change the structure
data<- melt(data, id.vars = "celltype", variable.name = "gene", value.name = "expression")
# combine both together 
df <- rbind(df,data) # here we inlcude three PP2A, three asSOFL1,three replicate for NRT1.9
# check the structure 
str(df)
df$celltype <- as.factor(df$celltype)
levels(df$celltype)[9:11] <- 'Phloem'
levels(df$celltype)[9] <- 'Xylem'
levels(df$celltype)[8] <- 'Pericycle'
df$celltype <- factor(df$celltype,levels=c("Epidermis", "Cortex", "Endodermis", "Pericycle", "Phloem", "Procambium", "Xylem"))

df$gene <- factor(df$gene,levels=c("NRT1.9","PP2A","asSOFL1"))
# save data write_xlsx(df,'three_gene_expression_smFISHcelltype.xlsx')
# draw figures 
p <- ggplot(df, aes(x=celltype, y=expression)) + 
    geom_violin(aes(fill=celltype),trim = FALSE, scale = "width")+geom_boxplot(width=0.1)+ scale_fill_manual(values=c("#989B7A", "#E2E8AA", "#5D9C76", "#A0D5B2", "#684FA1","#D1EBE7",'#3A4A92'))+facet_grid(gene~., scales = 'free_y')
#scale_fill_manual(values=c("#A68BC2", "#85C680", "#D09B7E", "#79AED2", "#EFB266","#DCCA83",'#EE7677'))
#scale_fill_manual(values=c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2","#8491B4B2",'#91D1C2B2'))

p <- p + theme(panel.grid.major = element_blank(),axis.line = element_line(colour = "black"), 
               panel.background = element_blank(),
               axis.text = element_text(size = 20),
               axis.title = element_text(size =20),
               axis.title.x = element_blank(),  # Remove x-axis title
               axis.text.x = element_text(size = 20, angle = 35, hjust = 1, vjust = 1),  # Set the size and angle of x-axis text
               plot.title = element_text(colour = "black",size = 25, face = "bold",hjust = 0.5, vjust = 0.5),
               legend.text = element_text(size = 16 ),
               panel.grid = element_blank(),legend.title  = element_blank(),strip.background.y = element_blank(),strip.text.y = element_text(size=16,face='bold.italic'),panel.spacing.y = unit(0.5,'cm'),
               legend.position = " ")+
    labs(title = "cryo-smFISH",
         x = "Cell type", y = "Transcript count per cell")


         # stastics analysis
test <- list()
k<- c()
for (i in c('NRT1.9','PP2A','asSOFL1')) {
    print(i)
    a1 <- aov(df[df$gene == i, ]$expression ~ df[df$gene == i, ]$celltype)
    test[[i]] <- TukeyHSD(a1) 
    m <- as.data.frame( test[[i]]$`df[df$gene == i, ]$celltype`)
m$Significance <- ifelse(m$`p adj` > 0.05, "ns", 
                      ifelse(m$`p adj` <= 0.05 & m$`p adj` > 0.01, "*", 
                      ifelse(m$`p adj` <= 0.01 & m$`p adj` > 0.001, "**",
                      ifelse(m$`p adj` <= 0.001 & m$`p adj` > 0.0001, "***",
                       "****"))))
k[[i]]<- m
#save
        write.xlsx(k[[i]],paste0('/Volumes/Expansion/zhang_et_al_2024_finaldata/Results/FigS8/D/smFISH', i,'_statistics.xlsx'))
}

 ggsave("smFISH_celltype_three_genes.pdf",p,width = 8, height = 6.5)

