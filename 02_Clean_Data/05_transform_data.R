### LOAD IN PACKAGES ##########################################################
library(vroom)
library(phyloseq)
library(robCompositions)
library(vegan)
library(ggplot2)
library(gridExtra)

### SCRIPT SETUP ##############################################################
date <- format(Sys.Date(),"_%Y%m%d")
pwd <- "/projectnb/talbot-lab-data/Katies_data/Street_Trees_Dysbiosis/"
amplicon <- "16S" # options: 16S or ITS
yourname <- "atherton" # user's last name for file storage purposes
edit_metadata <- "N" # options: Y or N
check_batch_effect <- "N" # options: Y or N
batch_correct <- "N" # options: Y or N
rarefy <- "N" # options: Y or N
clr_transform <- "N" # options: Y or N
z_transform <- "N" # options: Y or N
aitchison_distance <- "N" # options: Y or N

setwd(pwd)
source("00_functions.R")

### LOAD IN DATA ##############################################################
setwd(paste0(pwd, "02_Clean_Data/04_Decontam_ASV_Tables/", amplicon))
ps_all <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                       "_allsampletypes_phyloseq_decontam_"), 
                       ".RDS")

### CHECK FOR BATCH EFFECT ####################################################
setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon))
if(check_batch_effect == "Y"){
  check_batch_effect(ps_all, amplicon, yourname, date)
}

### BATCH CORRECTION ##########################################################
if(batch_correct == "Y"){
  
}

### RAREFY DATA ###############################################################
if(rarefy == "Y"){
  
}

### CLR-TRANSFORM #############################################################
if(clr_transform == "Y"){
  
}

### Z-SCORE TRANSFORM #########################################################
if(z_transform == "Y"){
  
}

### CALCULATE AITCHISON DISTANCE MATRIX #######################################
if(aitchison_distance == "Y"){
  
}

