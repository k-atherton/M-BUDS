### LOAD IN PACKAGES ##########################################################
library(phyloseq)

### SCRIPT SETUP ##############################################################
args <- commandArgs(trailingOnly = TRUE)
date <- format(Sys.Date(),"_%Y%m%d")
yourname <- args[1] # user's last name for file storage purposes
amplicon <- args[2] # options: 16S or ITS
edit_metadata <- args[3] # options: Y or N

### CHECK FOR AMPLICON TYPE ###################################################
if(amplicon %in% c("16s", "its", "16S", "ITS")){
  if(amplicon %in% c("16s", "16S")){
    amplicon <- "16S"
  } else{
    amplicon <- "ITS"
  }
  ### CHECK FOR EDITING METADATA ###############################################
  if(edit_metadata %in% c("Y", "N", "y", "n")){
    setwd("/projectnb/talbot-lab-data/Katies_data/Street_Trees_Dysbiosis/02_Clean_Data")
    source("00_functions.R")
    ### READ IN FORMATTED ASV TABLES ###########################################
    setwd(paste0("02_DADA2_ASV_Tables/",amplicon))
    leaf_raw <- read_in_file(getwd(), paste0(yourname, "_", amplicon, 
                               "_ASV_table_leaf_raw_"), ".csv")
    root_raw <- read_in_file(getwd(), paste0(yourname, "_", amplicon, 
                                             "_ASV_table_root_raw_"), ".csv")
    msoil_raw <- read_in_file(getwd(), paste0(yourname, "_", amplicon, 
                                             "_ASV_table_msoil_raw_"), ".csv")
    osoil_raw <- read_in_file(getwd(), paste0(yourname, "_", amplicon, 
                                             "_ASV_table_osoil_raw_"), ".csv")
    
    write.csv(leaf_raw_tax, paste0(getwd(), yourname, "_", amplicon, 
                                   "_taxonomy_leaf_raw_", date, ".csv"))
    write.csv(root_raw_tax, paste0(getwd(), yourname, "_", amplicon, 
                                   "_taxonomy_root_raw_", date, ".csv"))
    write.csv(msoil_raw_tax, paste0(getwd(), yourname, "_", amplicon, 
                                    "_taxonomy_msoil_raw_", date, ".csv"))
    write.csv(osoil_raw_tax, paste0(getwd(), yourname, "_", amplicon, 
                                    "_taxonomy_osoil_raw_", date, ".csv"))
  } else{
    print("Error: metadata trailing argument not recognized. \nOptions: Y (yes, edit the metadata file with raw DADA2 data) or N (no, do not edit metadata file)")
  }
} else{
  print("Error: amplicon trailing argument not recognized. \nOptions: 16S (bacterial) or ITS (fungal)")
}