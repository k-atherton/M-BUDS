read_in_file <- function(path, prefix, extension){
  files <- list.files(path, pattern=prefix)
  files <- files[endsWith(files, extension)]
  if(length(files) > 1){
    dates <- gsub(prefix, "", files)
    dates <- gsub("\\..*", "", dates)
    dates <- as.numeric(dates)
    dates <- dates[!is.na(dates)]
    dates <- dates[which(dates == max(dates))]
    if(endsWith(path, "/")){
      file <- paste0(path,prefix,dates,extension) # keep the most recent version of the file to use
    } else{
      file <- paste0(path,"/",prefix,dates,extension) # keep the most recent version of the file to use
    }
  } else{
    if(endsWith(path, "/")){
      file <- paste0(path,files)
    } else{
      file <- paste0(path,"/",files)
    }
  }
  if(extension == ".RDS"){
    data <- readRDS(file)
  } else{
    data <- vroom(file)    
  }
  return(data)
}

read_in_dada2_asv_table <- function(path, prefix, extension){
  data <- read_in_file(path, prefix, extension)
  colnames(data)[1] <- "ASV" # rename the first column for merging purposes
  return(data)
}

read_in_dada2_taxonomy <- function(path, prefix, extension){
  data <- read_in_file(path, prefix, extension)
  data <- data[,c(1,2)] # keep only the ASV ID and taxonomy information
  return(data)
}

subset_phyloseq <- function(ps_otu, ps_meta, type, nc_list){
  metadata$seq_bin <- NA
  if(type == "Leaf"){
    ps_subset <- subset_samples(ps_meta, sample_type == "Leaf") 
  } else if(type == "Root"){
    ps_subset <- subset_samples(ps_meta, sample_type == "Root") 
  } else if(type == "MSoil"){
    ps_subset <- subset_samples(ps_meta, (sample_type == "Soil") & (soil_horizon == "M"))
  } else if(type == "OSoil"){
    ps_subset <- subset_samples(ps_meta, (sample_type == "Soil") & (soil_horizon == "O"))
  }
  # add negative controls
  nc <- colnames(ps_otu)[grep(paste(nc_list, collapse = "|"), colnames(ps_otu))]
  ps_nc <- prune_samples(nc, ps_meta)
  ps_subset_w_nc <- merge_phyloseq(ps_subset, ps_nc)
  return(ps_subset_w_nc)
}

bin_seq_depth <- function(metadata){
  for(i in 1:nrow(metadata)){
    depth <- metadata$seq_count_dada2[i]
    if(depth < 10000){
      metadata$seq_bin[i] <- "0-10k"
    } else if(depth < 20000){
      metadata$seq_bin[i] <- "10-20k"
    } else if(depth < 30000){
      metadata$seq_bin[i] <- "20-30k"
    } else if(depth < 40000){
      metadata$seq_bin[i] <- "30-40k"
    } else if(depth < 50000){
      metadata$seq_bin[i] <- "40-50k"
    } else if(depth < 60000){
      metadata$seq_bin[i] <- "50-60k"
    } else if(depth < 70000){
      metadata$seq_bin[i] <- "60-70k"
    } else if(depth < 80000){
      metadata$seq_bin[i] <- "70-80k"
    } else if(depth < 90000){
      metadata$seq_bin[i] <- "80-90k"
    } else if(depth < 100000){
      metadata$seq_bin[i] <- "90-100k"
    } else{
      metadata$seq_bin[i] <- "100k"
    }
  }
  return(metadata)
}

plot_prefilter_seq_depth <- function(metadata, sample_type, proposed_threshold,
                                     date){
  ggplot(metadata, 
         aes(as.numeric(seq_count_dada2), 
             fill = sequencing_batch)) + 
    geom_histogram(binwidth = 1000) + 
    theme_bw() + 
    ggtitle(paste0(sample_type, " sequencing depth histogram")) + 
    geom_vline(aes(xintercept=median(as.numeric(seq_count_dada2))), 
               color="blue", 
               linetype="dashed") + 
    geom_vline(aes(xintercept=proposed_threshold), 
               color="red", 
               linetype="dashed") +
    labs(x = "Read count after DADA2", y = "count", fill = "Sequencing batch")
  
  ggsave(paste0(sample_type, "/", yourname, "_", amplicon, "_", sample_type, 
                "_sequencing_depth_prefilter", date,".png"), width = 9, 
         height = 5, units = "in", dpi = 300)
}

id_outliers_evaluate_seq_depth <- function(data, metadata, sample_type, yourname, 
                                    amplicon, date, rep){
  # make a transposed verison of the leaf sequencing data
  data_t <- as.data.frame(t(data))

  # calculate Aitchison distance for leaf samples
  aitch_data <- aDist(data_t+1, y = NULL) #need to add +1 otherwise all values are NA due to 0s in df
  aitch_data <- as.matrix(aitch_data)

  # test for sequencing batch and depth effect
  print(paste0("Testing for sequencing batch effect in ", sample_type, 
               " samples:"))
  print(adonis2(aitch_data~as.factor(metadata$sequencing_batch)))
  print(paste0("Testing for sequencing depth effect in ", sample_type, 
               " samples:"))
  print(adonis2(aitch_data~as.numeric(metadata$seq_count_dada2)))

  all.MDS <- metaMDS(aitch_data,k=2,zerodist="add")
  coordinates <- data.frame(scores(all.MDS))
  coordinates <- cbind(coordinates, metadata)

  # plot samples by different variables
  ggplot(coordinates,
         aes(x = NMDS1,
             y = NMDS2)) +
    geom_text(label = rownames(coordinates),
              size = 2) +
    stat_ellipse() +
    theme_bw()
  ggsave(paste0(yourname, "_", amplicon, "_", sample_type, 
                "_NMDS_identify_outliers_", rep, date, ".png"), width = 7, 
         height = 5, units = "in", dpi = 300)

  by_batch <- ggplot(coordinates,
         aes(x = NMDS1,
             y = NMDS2,
             col = as.factor(sequencing_batch))) +
    geom_point() +
    stat_ellipse() +
    theme_bw()

  by_seqdepth <- ggplot(coordinates,
         aes(x = NMDS1,
             y = NMDS2,
             col = as.factor(seq_bin))) +
    geom_point() +
    stat_ellipse() +
    theme_bw()

  by_treeage <- ggplot(coordinates,
         aes(x = NMDS1,
             y = NMDS2,
             col = as.factor(tree_age))) +
    geom_point() +
    stat_ellipse() +
    theme_bw()

  by_treepittype <- ggplot(coordinates,
         aes(x = NMDS1,
             y = NMDS2,
             col = as.factor(tree_pit_type))) +
    geom_point() +
    stat_ellipse() +
    theme_bw()

  by_species <- ggplot(coordinates,
         aes(x = NMDS1,
             y = NMDS2,
             col = as.factor(tree_species))) +
    geom_point() +
    stat_ellipse() +
    theme_bw()

  by_site <- ggplot(coordinates,
         aes(x = NMDS1,
             y = NMDS2,
             col = as.factor(site_name))) +
    geom_point() +
    stat_ellipse() +
    theme_bw()
  
  multipanel <- grid.arrange(by_batch, by_seqdepth, by_treeage, by_treepittype, 
                             by_species, by_site, nrow = 2, ncol = 3)
  
  ggsave(paste0(yourname, "_", amplicon, "_", sample_type, 
                "_NMDS_predrop_datastructure_", rep, date, ".png"), multipanel, 
         width = 21, height = 10, units = "in", dpi = 300)
}

test_drop_threshold <- function(data, metadata, sample_type, yourname, 
                                 amplicon, date, threshold){
  # drop samples < threshold
  metadata_drop <- metadata[metadata$seq_count_dada2 > threshold,]
  data_drop <- data[,colnames(data) %in% rownames(metadata_drop)]
  
  # make a transposed verison of the sequencing data
  data_t <- as.data.frame(t(data_drop))
  
  # calculate Aitchison distance for samples
  aitch <- aDist(data_t+1, y = NULL) #need to add +1 otherwise all values are NA due to 0s in df
  aitch <- as.matrix(aitch)
  
  # test for sequencing batch and depth effect
  print("Testing sequencing batch effect:")
  print(adonis2(aitch ~ as.factor(metadata_drop$sequencing_batch)))
  
  print("Testing sequencing depth effect:")
  adonis2(aitch ~ as.numeric(metadata_drop$seq_count_dada2))
  
  # calculate NMDS
  all.MDS <- metaMDS(aitch, k=2, zerodist="add")
  coordinates<-data.frame(scores(all.MDS))
  coordinates <- cbind(coordinates, metadata_drop)
  
  # plot samples by different variables
  by_batch <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2,
                                      col = as.factor(sequencing_batch))) +
    geom_point() + stat_ellipse() + theme_bw()
  
  by_seqdepth <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2,
                                         col = as.factor(seq_bin))) +
    geom_point() + stat_ellipse() + theme_bw()
  
  by_treeage <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2,
                                        col = as.factor(tree_age))) +
    geom_point() + stat_ellipse() + theme_bw()
  
  by_treepittype <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2,
                                            col = as.factor(tree_pit_type))) +
    geom_point() + stat_ellipse() + theme_bw()
  
  by_species <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                        col = as.factor(tree_species))) +
    geom_point() + stat_ellipse() + theme_bw()
  
  by_site <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                     col = as.factor(site_name))) +
    geom_point() + stat_ellipse() + theme_bw()
  
  multipanel <- grid.arrange(by_batch, by_seqdepth, by_treeage, by_treepittype, 
                             by_species, by_site, nrow = 2, ncol = 3)
  
  ggsave(paste0(yourname, "_", amplicon, "_", sample_type, 
                "_NMDS_drop", threshold, date, ".png"), multipanel, 
         width = 21, height = 10, units = "in", dpi = 300)
  
}

decontaminate_samples <- function(ps, sample_type, yourname, amplicon, date){
  # get the batch IDs
  metadata <- as.data.frame(as.matrix(ps@sam_data))
  batches <- unique(metadata$sequencing_batch[which(metadata$sample_type != "Negative Control")])
  
  pre_decontam_pngs <- list()
  contaminants_pngs <- list()
  
  j <- 1
  
  for(i in 1:length(batches)){
    batch <- batches[i]
    print(batch)
    # subset the batch
    batch_data <- prune_samples(sample_data(ps)$sequencing_batch == batch, ps)
    meta_batch <- as.data.frame(as.matrix(batch_data@sam_data))
    
    # order for rank
    meta_batch <- meta_batch[order(meta_batch$seq_count_dada2),]
    meta_batch$index <- seq(nrow(meta_batch))
    
    # plot the pre-decontam data
    pre_decontam_pngs[[i]] <- ggplot(meta_batch, aes(x = index, 
                                                     y = as.numeric(seq_count_dada2), 
                                                     color = is_control)) +
                             geom_point() + theme_bw() + 
                             ggtitle(paste0(sample_type, " ", batch))
    
    # run decontam
    if(length(unique(meta_batch$is_control)) > 1){
      decontam_batch <- isContaminant(batch_data, neg = "is_control", 
                                      method = "prevalence")
      
      # visualize contaminants
      pa <- transform_sample_counts(batch_data, function(abund) 1*(abund>0))
      neg <- prune_samples(sample_data(batch_data)$is_control == TRUE, 
                           batch_data)
      pos <- prune_samples(sample_data(batch_data)$is_control == FALSE, 
                           batch_data)
      df_batch <- data.frame(pos=taxa_sums(pos), neg=taxa_sums(neg), 
                             contaminant=decontam_batch$contaminant)
      contaminants_pngs[[j]] <- ggplot(df_batch, aes(x= neg, y = pos, 
                                                     color = contaminant)) +
        geom_point() + xlab("Prevalence (Negative Controls)") + 
        ylab("Prevalence (True Samples)") + theme_bw() +
        ggtitle(paste0(sample_type, " ", batch))
      
      j <- j + 1
    } else{
      decontam_batch <- isContaminant(batch_data, conc = "dna_conc", 
                                      method = "frequency")
      
      # visualize contaminants
      freq_contaminants <- plot_frequency(batch_data, 
                                          taxa_names(batch_data)[sample(which(decontam_batch$contaminant),
                                                                        sum(decontam_batch$contaminant, 
                                                                            na.rm = TRUE))], 
                                          conc="dna_conc") + 
                                          xlab("DNA Concentration")
      
      png(filename=paste0("Figures/", sample_type, "/", 
                          yourname, "_", amplicon, "_", sample_type, "_",
                          batch, "_decontam_frequencymethod", date, 
                          ".png"))
      plot(freq_contaminants)
      dev.off()
    }
    # identify contaminants
    taxonomy <- as.data.frame(as.matrix(batch_data@tax_table))
    contaminants <- taxonomy[which(rownames(taxonomy) %in% 
                                     rownames(decontam_batch)[which(decontam_batch$contaminant == TRUE)]),]
    
    # write information to table
    write.csv(contaminants, paste0("Contaminant_Taxonomy/", sample_type, "/", 
                                   yourname, "_", amplicon, "_", sample_type, 
                                   "_", batch, "_contaminants", date, 
                                   ".csv"))
    
    # remove contaminants
    decontam_data <- prune_taxa(!decontam_batch$contaminant, batch_data)
    
    # merge back into one dataset
    if(i > 1){
      final_data <- merge_phyloseq(final_data, decontam_data)
    } else{
      final_data <- decontam_data
    }
  }
  # save multipanel figures
  multipanel_pre_decontam <- do.call(grid.arrange, pre_decontam_pngs)
  multipanel_contaminants <- do.call(grid.arrange, contaminants_pngs)
  
  ggsave(paste0("Figures/", sample_type, "/",
                yourname, "_", amplicon, "_", sample_type, 
                "_predecontam", date, ".png"), multipanel_pre_decontam, 
         width = 21, height = 10, units = "in", dpi = 300)
  
  ggsave(paste0("Figures/", sample_type, "/",
                yourname, "_", amplicon, "_", sample_type, 
                "_decontam_prevalencemethod", date, ".png"), 
         multipanel_contaminants, width = 21, height = 10, units = "in", 
         dpi = 300)
  
  # remove negative controls
  if(sample_type == "leaf"){
    final_data <- subset_samples(final_data, sample_type == "Leaf")
  } else if(sample_type == "root"){
    final_data <- subset_samples(final_data, sample_type == "Root")
  } else {
    final_data <- subset_samples(final_data, sample_type == "Soil")
  }
  
  if(amplicon == "16S"){
    # remove chloroplasts, mitochondria, and anything not a bacteria
    final_data <- subset_taxa(final_data, Kingdom != "Unassigned")
    final_data <- subset_taxa(final_data, Kingdom != "d__Eukaryota")
    final_data <- subset_taxa(final_data, Family != " f__Mitochondria")
    final_data <- subset_taxa(final_data, Order != " o__Chloroplast")
  } else{
    # remove anything not a fungi
    final_data <- subset_taxa(final_data, Kingdom == "k__Fungi")
  }
  
  # write ASV table to CSV
  asv_table <- as.data.frame(as.matrix(final_data@otu_table))
  write.csv(asv_table, paste0(yourname, "_", amplicon, "_", sample_type, 
                              "ASV_table_decontam", date, ".csv"))
  return(final_data)
}

plot_decontam_seq_depth <- function(metadata, sample_type, date){
  ggplot(metadata, 
         aes(as.numeric(seq_count_decontam), 
             fill = sequencing_batch)) + 
    geom_histogram(binwidth = 1000) + 
    theme_bw() + 
    ggtitle(paste0(sample_type, " sequencing depth histogram")) + 
    geom_vline(aes(xintercept=median(as.numeric(seq_count_dada2))), 
               color="blue", 
               linetype="dashed") + 
    geom_vline(aes(xintercept=median(as.numeric(seq_count_decontam), 
                                     na.rm = TRUE)), 
               color="red", 
               linetype="dashed") +
    labs(x = "Read count after decontam", y = "count", 
         fill = "Sequencing batch")
  
  ggsave(paste0(sample_type, "/", yourname, "_", amplicon, "_", sample_type, 
                "_sequencing_depth_postdecontam", date,".png"), width = 9, 
         height = 5, units = "in", dpi = 300)
}

check_batch_effect <- function(ps, amplicon, yourname, date){
  # extract asv table
  sequence_data <- as.data.frame(as.matrix(ps@otu_table))
  
  # extract metadata
  metadata <- as.data.frame(as.matrix(ps@sam_data))
  metadata$sample_type[which(metadata$sample_type == "Soil")] <- paste0(metadata$soil_horizon[which(metadata$sample_type == "Soil")], 
                                                               " Soil")
  
  # make a transposed version of asv table
  sequence_data_t <- as.data.frame(t(sequence_data))
  
  # calculate Aitchison distance for o samples
  aitch <- aDist(sequence_data_t+1, y = NULL) #need to add +1 otherwise all values are NA due to 0s in df
  aitch <- as.matrix(aitch)
  
  # test for sequencing batch and depth effect
  print(adonis2(aitch ~ as.factor(metadata$sequencing_batch)))
  print(adonis2(aitch ~ as.factor(metadata$tree_type)))
  print(adonis2(aitch ~ as.factor(metadata$tree_age)))
  print(adonis2(aitch ~ as.factor(metadata$tree_pit_type)))
  print(adonis2(aitch ~ as.factor(metadata$tree_species)))
  print(adonis2(aitch ~ as.factor(metadata$site_name)))
  print(adonis2(aitch ~ as.factor(metadata$sample_type)))
  
  
  all.MDS <- metaMDS(aitch, k=2, zerodist="add")
  coordinates <- data.frame(scores(all.MDS))
  coordinates <- cbind(coordinates, metadata)
  
  # plot samples by different variables
  seq_batch <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                          col = as.factor(sequencing_batch))) + 
    geom_point() + stat_ellipse() + theme_bw()
  tree_type <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                       col = as.factor(tree_type))) + 
    geom_point() + stat_ellipse() + theme_bw()
  tree_age <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                      col = as.factor(tree_age))) + 
    geom_point() + stat_ellipse() + theme_bw()
  tree_pit_type <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                           col = as.factor(tree_pit_type))) + 
    geom_point() + stat_ellipse() + theme_bw()
  tree_species <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                          col = as.factor(tree_species))) + 
    geom_point() + stat_ellipse() + theme_bw()
  site_name <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                       col = as.factor(site_name))) + 
    geom_point() + stat_ellipse() + theme_bw()
  sample_type <- ggplot(coordinates, aes(x = NMDS1, y = NMDS2, 
                                         col = as.factor(sample_type))) + 
    geom_point() + stat_ellipse() + theme_bw()
  
  multipanel <- grid.arrange(seq_batch, tree_type, tree_age, tree_pit_type,
                             tree_species, site_name, sample_type, nrow = 2, 
                             ncol = 4)
  
  ggsave(paste0(yourname, "_", amplicon, "_check_batch_effect", date, ".png"), 
         multipanel, width = 28, height = 10, units = "in", dpi = 300)
}

batch_correct <- function(ps, amplicon, yourname, date){
  # get asv table from phyloseq object
  asv_table <- as.data.frame(as.matrix(ps@otu_table))
  
  # remove any taxa with 0 counts and any samples under 1000 counts as 
  # ComBat_seq will not work if you have samples with 10k+ counts AND really 
  # low counts in the same dataset
  asv_table <- asv_table[which(rowSums(asv_table) > 0), 
                         which(colSums(asv_table) > 1000)]
  asv_table <- asv_table[,order(colnames(asv_table))]
  asv_table <- as.matrix(asv_table)
  
  # asv_table <- as.matrix(asv_table)
  # asv_table <- asv_table %>% mutate_if(is.character, as.numeric)
  
  # get metadata from phyloseq object
  metadata <- as.data.frame(as.matrix(ps@sam_data))
  
  # keep metadata samples which are still in the asv table
  metadata <- metadata[which(rownames(metadata) %in% colnames(asv_table)),]
  
  # order the metadata by sample name
  metadata <- metadata[order(row.names(metadata)),]
  
  # batch is defined as the sample type, as the samples were extracted and 
  # amplified with different methods, but the sequencing batch did not have
  # a batch effect
  batch <- paste0(metadata$sample_type, " ", metadata$soil_horizon)
  # for leaves, remove the soil horizon NA
  batch <- gsub(" NA", "", batch) 
  # for roots, just keep roots as they were all amplified the same way
  batch <- gsub("Root.*", "Root", batch) 
  
  # the covariates we want to maintain the signature of are tree age and tree 
  # pit type
  age <- metadata$tree_age
  type <- metadata$tree_pit_type
  covar_mod <- data.frame(age, type)
  
  # batch correct with ComBat_seq
  batch_corrected <- ComBat_seq(asv_table, batch, group=NULL, covar_mod)
  
  # visualize how ComBat_seq changed the sample sums
  pre_combatseq <- as.data.frame(colSums(asv_table))
  post_combatseq <- as.data.frame(colSums(batch_corrected))
  
  pre_figure <- ggplot(pre_combatseq, aes(x = `colSums(asv_table)`)) + 
    geom_histogram() + xlab("Pre-Batch Correction Sample Count") + 
    ylab("Count") + theme_bw() 
  post_figure <- ggplot(post_combatseq, aes(x = `colSums(batch_corrected)`)) +
    geom_histogram() + xlab("Post-Batch Correction Sample Count") + 
    ylab("Count") + theme_bw()
  
  multipanel <- grid.arrange(pre_figure, post_figure, ncol = 2, nrow = 1)
  
  ggsave(paste0(yourname, "_", amplicon, "_batchcorrection_change", date, 
                ".png"), multipanel, width = 14, height = 5, units = "in", 
         dpi = 300)
  
  return(batch_corrected)
}

rarefy_data <- function(ps, sample_type, amplicon, yourname, date){
  # extract asv table from pre-rarefied data
  pre_rare_asv <- as.data.frame(as.matrix(ps@otu_table))
  
  # rarefy data
  set.seed(1)
  if(sample_type == "leaf"){
    ps_rare <- rarefy_even_depth(ps, rngseed=1, sample.size=200, replace=F)
  }
  else{
    ps_rare <- rarefy_even_depth(ps, rngseed=1, 
                                 sample.size=min(sample_sums(ps)), replace=F)
  }
  
  # extract the asv table from the rarefied data
  rare_asv <- as.data.frame(as.matrix(ps_rare@otu_table))
  
  # write to files
  write.csv(rare_asv, paste0(yourname, "_", amplicon, "_", sample_type, 
                             "_rarefied_ASV_table", date, ".csv"))
 
  # format data for plotting
  rare_asv$asv <- rownames(rare_asv)
  pre_rare_asv$asv <- rownames(pre_rare_asv)
  
  data_long <- pivot_longer(rare_asv, !asv, names_to = "sample", 
                            values_to = "rarefied_count")
  data_long_pre <- pivot_longer(pre_rare_asv, !asv, names_to = "sample", 
                                values_to = "pre_rarefied_count")

  data_compare <- merge(data_long, data_long_pre, by = c("asv", "sample"))
  
  tax <- as.data.frame(as.matrix(ps@tax_table))
  tax$asv <- rownames(tax)
  
  data_compare <- merge(data_compare, tax, by = "asv")
  
  # plot change
  by_phylum <- ggplot(data_compare, aes(x = pre_rarefied_count, 
                                        y = rarefied_count, col = Phylum)) + 
    geom_point() + theme_bw()
  
  by_sample <- ggplot(data_compare, aes(x = pre_rarefied_count, 
                                        y = rarefied_count, col = sample)) + 
    geom_point() + theme_bw()
  multipanel <- grid.arrange(by_phylum, by_sample, ncol = 1, nrow = 2)
  ggsave(paste0("Figures/", yourname, "_", amplicon, "_", sample_type, 
                "_pre_post_rarefaction_counts", date, ".png"), multipanel, 
         height = 10, width = 7, units = "in", dpi = 300)
  
  # plot rarefaction curves
  if(sample_type == "leaf"){
    png(filename=paste0(yourname, "_", amplicon, "_", sample_type, 
                        "rarefaction_curve", date, ".png"))
    rarecurve(t(as.data.frame(as.matrix(ps@otu_table))), step = 20, 
              sample = 200, col = "blue", label = FALSE)
    dev.off()
  } else{
    png(filename=paste0("Figures/", yourname, "_", amplicon, "_", sample_type, 
                        "rarefaction_curve", date, ".png"))
    rarecurve(t(as.data.frame(as.matrix(ps@otu_table))), step = 20, 
              sample = min(colSums(as.data.frame(as.matrix(ps@otu_table)))), 
              col = "blue", label = FALSE)
    dev.off()
  }
}

clr_data <- function(ps, sample_type, amplicon, yourname, date){
  # extract the ASV table
  asv_table <- as.data.frame(as.matrix(ps@otu_table))
  
  # do the CLR transformation
  clr_transformed <- as.data.frame(clr(asv_table+1))
  
  # write to file
  write.csv(clr_transformed, paste0(yourname, "_", amplicon, "_", sample_type, 
                             "_CLRtransformed_ASV_table", date, ".csv"))
}

zscore_data <- function(ps, sample_type, amplicon, yourname, date){
  # extract the ASV table
  asv_table <- as.data.frame(as.matrix(ps@otu_table))
  
  # do the Z-Score transformation
  z_transformed <- as.data.frame(scale(asv_table, center = TRUE, scale = TRUE))
  
  # write to file
  write.csv(z_transformed, paste0(yourname, "_", amplicon, "_", sample_type, 
                                    "_ZScoretransformed_ASV_table", date, ".csv"))
}

aitchison_data <- function(ps, sample_type, amplicon, yourname, date){
  # extract the ASV table
  asv_table <- as.data.frame(as.matrix(ps@otu_table))
  
  # calculate Aitchison distance matrix
  aitchison_distance_matrix <- aDist(t(asv_table+1), y = NULL)
  aitchison_distance_matrix <- as.matrix(aitchison_distance_matrix)
  
  # write to file
  write.csv(aitchison_distance_matrix, paste0(yourname, "_", amplicon, "_", 
                                              sample_type, 
                                              "_Aitchison_distance_matrix", 
                                              date, ".csv"))
}
