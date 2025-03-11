# load in packages
library(vroom)

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
  data <- vroom(file)
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
  if(type == "Leaf"){
    ps_subset <- subset_samples(ps_meta, sample_type == "Leaf") 
  } else if(type == "Root"){
    ps_subset <- subset_samples(ps_meta, sample_type == "Root") 
  } else if(type == "M Soil"){
    ps_subset <- subset_samples(ps_meta, (sample_type == "Soil") & (soil_horizon == "M"))
  } else if(type == "O Soil"){
    ps_subset <- subset_samples(ps_meta, (sample_type == "Soil") & (soil_horizon == "O"))
  }
  # add negative controls
  nc <- colnames(ps_otu)[grep(paste(nc_list, collapse = "|"), colnames(ps_otu))]
  ps_nc <- subset_samples(ps_meta, sample_names(ps_meta) %in% nc)
  ps_subset_w_nc <- merge_phyloseq(ps_subset, ps_nc)
  return(ps_subset_w_nc)
}