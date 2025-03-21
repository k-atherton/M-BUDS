### LOAD IN PACKAGES ##########################################################
library(vroom)
library(phyloseq)
library(robCompositions)
library(vegan)
library(ggplot2)
library(gridExtra)
library(sva)
library(dplyr)
library(tidyr)
library(compositions)

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
ps_leaf <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                        "_leaf_phyloseq_decontam_"), 
                        ".RDS")
ps_root <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                        "_root_phyloseq_decontam_"), 
                        ".RDS")
ps_msoil <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                         "_msoil_phyloseq_decontam_"), 
                         ".RDS")
ps_osoil <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                         "_osoil_phyloseq_decontam_"), 
                         ".RDS")

### CHECK FOR BATCH EFFECT ####################################################
setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon))
if(check_batch_effect == "Y"){
  check_batch_effect(ps_all, amplicon, yourname, date)
}

### BATCH CORRECTION ##########################################################
setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon))
if(batch_correct == "Y"){
  corrected <- batch_correct(ps_all, amplicon, yourname, date)
}

write.csv(corrected, paste0("Batch_Corrected_ASV_Tables/", yourname, "_", 
                            amplicon, "_batch_corrected_ASV_table", date, 
                            ".csv"))
### RAREFY DATA ###############################################################
setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon, 
             "/Rarefied_ASV_Tables"))
# rarefy separately for each sample type
if(rarefy == "Y"){
  rarefy_data(ps_leaf, "leaf", amplicon, yourname, date)
  rarefy_data(ps_root, "root", amplicon, yourname, date)
  rarefy_data(ps_msoil, "msoil", amplicon, yourname, date)
  rarefy_data(ps_osoil, "osoil", amplicon, yourname, date)
}

### CLR-TRANSFORM #############################################################
setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon, 
             "/CLR_Transformed_Tables"))
# clr-transform separately for each sample type and all together
if(clr_transform == "Y"){
  clr_data(ps_leaf, "leaf", amplicon, yourname, date)
  clr_data(ps_root, "root", amplicon, yourname, date)
  clr_data(ps_msoil, "msoil", amplicon, yourname, date)
  clr_data(ps_osoil, "osoil", amplicon, yourname, date)
  clr_data(ps_all, "allsampletypes", amplicon, yourname, date)
}

### Z-SCORE TRANSFORM #########################################################
setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon, 
             "/ZScore_Transformed_Tables"))
# z-score transform separately for each sample type and all together
if(z_transform == "Y"){
  zscore_data(ps_leaf, "leaf", amplicon, yourname, date)
  zscore_data(ps_root, "root", amplicon, yourname, date)
  zscore_data(ps_msoil, "msoil", amplicon, yourname, date)
  zscore_data(ps_osoil, "osoil", amplicon, yourname, date)
  zscore_data(ps_all, "allsampletypes", amplicon, yourname, date)
}

### CALCULATE AITCHISON DISTANCE MATRIX #######################################
setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon, 
             "/Aitchison_Distance_Tables"))
# calculate aitchison_distance sparately for each sample type and all together
# for the all together distance matrix, use the batch-corrected data
if(aitchison_distance == "Y"){
  aitchison_data(ps_leaf, "leaf", amplicon, yourname, date)
  aitchison_data(ps_root, "root", amplicon, yourname, date)
  aitchison_data(ps_msoil, "msoil", amplicon, yourname, date)
  aitchison_data(ps_osoil, "osoil", amplicon, yourname, date)
  aitchison_data(ps_all, "allsampletypes", amplicon, yourname, date)
}

