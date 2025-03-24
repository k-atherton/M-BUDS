### LOAD IN PACKAGES ##########################################################
print("LOADING IN PACKAGES:")
library(optparse)
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
print("SETTING UP SCRIPT:")
date <- format(Sys.Date(),"_%Y%m%d")
pwd <- "/projectnb/talbot-lab-data/Katies_data/Street_Trees_Dysbiosis/"

option_list = list(
  make_option(c("-a", "--amplicon"), type="character", default="16S", 
              help="amplicon dataset to filter; options: 16S or ITS [default= %default]", 
              metavar="character"),
  make_option(c("-n", "--name"), type="character", default="atherton", 
              help="last name for output file naming scheme [default= %default]", 
              metavar="character"),
  make_option(c("-e", "--edit"), type="character", default="N", 
              help="do you want to edit the metadata file? options: Y or N [default= %default]", 
              metavar="character"),
  make_option(c("-c", "--checkbatch"), type="character", default="N", 
              help="do you want to check for a batch effect? options: Y or N [default= %default]", 
              metavar="character"),
  make_option(c("-b", "--batchcorrect"), type="character", default="N", 
              help="do you want to run batch correction? options: Y or N [default= %default]", 
              metavar="character"),
  make_option(c("-r", "--rarefy"), type="character", default="N", 
              help="do you want to rarefy your data? options: Y or N [default= %default]", 
              metavar="character"),
  make_option(c("-l", "--clr"), type="character", default="N", 
              help="do you want to clr-transform your data? options: Y or N [default= %default]", 
              metavar="character"),
  make_option(c("-z", "--zscore"), type="character", default="N", 
              help="do you want to z-score transform your data? options: Y or N [default= %default]", 
              metavar="character"),
  make_option(c("-d", "--aitchisondistance"), type="character", default="N", 
              help="do you want to calculate the aitchison distance matrix for your samples? options: Y or N [default= %default]", 
              metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

amplicon <- opt$amplicon
yourname <- opt$name
edit_metadata <- opt$edit
check_batch_effect <- "N" # options: Y or N
batch_correct <- "N" # options: Y or N
rarefy <- "N" # options: Y or N
clr_transform <- "N" # options: Y or N
z_transform <- "N" # options: Y or N
aitchison_distance <- "N" # options: Y or N

setwd(pwd)
source("00_functions.R")

### LOAD IN DATA ##############################################################
print("LOADING IN DATA:")
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
  print("CHECKING FOR BATCH EFFECT:")
  check_batch_effect(ps_all, amplicon, yourname, date)
}

### BATCH CORRECTION ##########################################################
setwd(paste0(pwd, "02_Clean_Data/05_Transform_Data/", amplicon))
if(batch_correct == "Y"){
  print("RUNNING BATCH CORRECTION: ")
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
  print("RAREFYING DATA: ")
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
  print("CLR-TRANSFORMING DATA: ")
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
  print("Z-SCORE TRANSFORMING DATA:")
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
  print("CALCULATING AITCHISON DISTANCE MATRIX FOR DATA:")
  aitchison_data(ps_leaf, "leaf", amplicon, yourname, date)
  aitchison_data(ps_root, "root", amplicon, yourname, date)
  aitchison_data(ps_msoil, "msoil", amplicon, yourname, date)
  aitchison_data(ps_osoil, "osoil", amplicon, yourname, date)
  aitchison_data(ps_all, "allsampletypes", amplicon, yourname, date)
}

