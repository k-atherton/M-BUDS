### LOAD IN PACKAGES ##########################################################
library(vroom)
library(phyloseq)
library(gridExtra)
library(decontam)

### SCRIPT SETUP ##############################################################
date <- format(Sys.Date(),"_%Y%m%d")
pwd <- "/projectnb/talbot-lab-data/Katies_data/Street_Trees_Dysbiosis/"
amplicon <- "ITS" # options: 16S or ITS
yourname <- "atherton" # user's last name for file storage purposes
edit_metadata <- "N" # options: Y or N

setwd(pwd)
source("00_functions.R")

### READ IN PHYLOSEQ OBJECTS ##################################################
setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon))
ps_leaf <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                        "_leaf_phyloseq_filteredsamples_"), 
                        ".RDS")
ps_root <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                        "_root_phyloseq_filteredsamples_"), 
                        ".RDS")
ps_msoil <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                        "_msoil_phyloseq_filteredsamples_"), 
                        ".RDS")
ps_osoil <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                        "_osoil_phyloseq_filteredsamples_"), 
                        ".RDS")

### RUN DECONTAM ##############################################################
setwd(paste0(pwd, "02_Clean_Data/04_Decontam_ASV_Tables/", amplicon))
ps_leaf_decontam <- decontaminate_samples(ps_leaf, "leaf", yourname, amplicon, 
                                          date)
ps_root_decontam <- decontaminate_samples(ps_leaf, "root", yourname, amplicon, 
                                          date)
ps_msoil_decontam <- decontaminate_samples(ps_leaf, "msoil", yourname, amplicon, 
                                          date)
ps_osoil_decontam <- decontaminate_samples(ps_leaf, "osoil", yourname, amplicon, 
                                          date)

### MERGE ALL DECONTAMINATED SAMPLE TYPES INTO ONE PHYLOSEQ OBJECT ############
ps_decontam <- merge_phyloseq(ps_leaf_decontam, ps_root_decontam)
ps_decontam <- merge_phyloseq(ps_decontam, ps_msoil_decontam)
ps_decontam <- merge_phyloseq(ps_decontam, ps_osoil_decontam)

asv_table <- as.data.frame(as.matrix(ps_decontam@otu_table))
write.csv(asv_table, paste0(yourname, "_", amplicon, "_", 
                            "allsampletypes_ASV_table_decontam", date, ".csv"))

### WRITE ALL PHYLOSEQ OBJECTS AS RDS TO SCC ##################################
saveRDS(ps_decontam, paste0(yourname, "_", amplicon, "_",
                            "allsampletypes_phyloseq_decontam", date, ".RDS"))
saveRDS(ps_leaf_decontam, paste0(yourname, "_", amplicon, "_",
                            "leaf_phyloseq_decontam", date, ".RDS"))
saveRDS(ps_root_decontam, paste0(yourname, "_", amplicon, "_",
                            "root_phyloseq_decontam", date, ".RDS"))
saveRDS(ps_msoil_decontam, paste0(yourname, "_", amplicon, "_",
                            "msoil_phyloseq_decontam", date, ".RDS"))
saveRDS(ps_osoil_decontam, paste0(yourname, "_", amplicon, "_",
                            "osoil_phyloseq_decontam", date, ".RDS"))