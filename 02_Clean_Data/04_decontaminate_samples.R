### LOAD IN PACKAGES ##########################################################
library(vroom)
library(phyloseq)
library(gridExtra)
library(decontam)

### SCRIPT SETUP ##############################################################
date <- format(Sys.Date(),"_%Y%m%d")
pwd <- "/projectnb/talbot-lab-data/Katies_data/Street_Trees_Dysbiosis/"
amplicon <- "16S" # options: 16S or ITS
yourname <- "atherton" # user's last name for file storage purposes
edit_metadata <- "N" # options: Y or N

setwd(pwd)
source("00_functions.R")

### READ IN PHYLOSEQ OBJECTS ##################################################
setwd(paste0(pwd, "02_Clean_Data/03_Filter_Samples_ASV_Tables/", amplicon))
ps_leaf <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                        "_leaf_phyloseq_filteredsamples_withnc_"), 
                        ".RDS")
ps_root <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                        "_root_phyloseq_filteredsamples_withnc_"), 
                        ".RDS")
ps_msoil <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                        "_msoil_phyloseq_filteredsamples_withnc_"), 
                        ".RDS")
ps_osoil <- read_in_file(getwd(), paste0("atherton_", amplicon, 
                                        "_osoil_phyloseq_filteredsamples_withnc_"), 
                        ".RDS")

### RUN DECONTAM ##############################################################
setwd(paste0(pwd, "02_Clean_Data/04_Decontam_ASV_Tables/", amplicon))
ps_leaf_decontam <- decontaminate_samples(ps_leaf, "leaf", yourname, amplicon, 
                                          date)
ps_root_decontam <- decontaminate_samples(ps_root, "root", yourname, amplicon, 
                                          date)
ps_msoil_decontam <- decontaminate_samples(ps_msoil, "msoil", yourname, amplicon, 
                                          date)
ps_osoil_decontam <- decontaminate_samples(ps_osoil, "osoil", yourname, amplicon, 
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

### WRITE POST-DECONTAM SEQUENCE COUNT TO METADATA ############################
if(edit_metadata == "Y") {
  metadata <- as.data.frame(as.matrix(ps_decontam@sam_data))
  
  # add decontam sequence count column
  metadata$seq_count_decontam <- NA
  
  # sum the post-dada2 processing sequence counts in each sample
  decontam_seq_count <- colSums(ps_decontam@otu_table)
  
  # record sequence counts in metadata table
  for (i in 1:length(decontam_seq_count)) {
    sample_name <- names(decontam_seq_count)[i]
    metadata$seq_count_decontam[which(metadata$sample_name == sample_name)] <-
      decontam_seq_count[i]
  }
  
  # write this data to the file
  write.csv(metadata, paste0(getwd(), yourname, "_sample_metadata_", amplicon, 
                             date, ".csv"), row.names = FALSE)
}

### VISUALIZE HOW DECONTAM IMPACTED SEQUENCING DEPTH ##########################
setwd(paste0(pwd,"/02_Clean_Data/04_Decontam_ASV_Tables/",amplicon,"/Figures"))
plot_decontam_seq_depth(ps_leaf_decontam@sam_data, "leaf", date)
plot_decontam_seq_depth(ps_root_decontam@sam_data, "root", date)
plot_decontam_seq_depth(ps_msoil_decontam@sam_data, "msoil", date)
plot_decontam_seq_depth(ps_osoil_decontam@sam_data, "osoil", date)
