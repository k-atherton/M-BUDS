# Kathryn Atherton Dysbiosis in Urban Street Trees (Chapter 2)
The purpose of this project is to understand how urbanization impacts the tree-associated microbiome. 

## Requirements
To run this analysis you will need the following R packages (version I used for this analysis):
- dplyr (1.1.4)
- ggplot2 (3.5.1)
- gridExtra (2.3)
- phyloseq (1.44.0)
- plyr (1.8.9)
- robCompositions (2.4.1)
- vegan (2.6-4)
- vroom (1.6.5)

Additionally, the initial sequence data cleaning step is run with the BU16S pipeline, created by Michael Silverstein and published in [Nature Ecology and Evolution](https://www.nature.com/articles/s41559-024-02440-6). More information about the pipeline can be found [here](https://github.com/Boston-University-Microbiome-Initiative/BU16s). 

## About the directories and scripts
- `00_functions.R`: contains the custom functions used repeatedly throughout the analyses.

- `01_Collect_Data`: this is where I store my raw data/metadata files. These are not currently shared on GitHub, and no scripts are stored here.
- `02_Clean_Data`: this is where I store my cleaned data/metadata files and the files that I use to make them.
    - `01_DADA2`: this is where I store my scripts to run DADA2 through the BU16S pipeline to perform QC and taxonomic annotation on my sequence data.
      - There is one directory per 16S/ITS sequencing run here. Each directory contains the qsub script required to run the BU16S pipeline. The pipeline outputs are also stored here, but they are not currently shared on GitHub.  
    - `02_format_raw_dada2_asv_tables.R`: this script takes the BU16S taxonomy and ASV table outputs and formats them for downstream cleaning steps. It merges all sequencing runs into one ASV table and then writes the post-DADA2 sequence count information into the metadata table. It also merges the taxonomy assignments into one table and formats the taxonomic information into separate columns for each taxonomic rank. Finally, it separates the ASV tables into one table for each sample type (leaf, root, upper soil horizon, lower soil horizon) and saves a csv and phyloseq object for each sample type for downstream analyis.
        - This file takes in the following arguments, to be changed in the script: a name (used for file storage purposes), an amplicon ID (either 16S or ITS, for determining which dataset to format), and an edit metadata argument (Y or N, for determining whether the sequence count information should be written into the metadata table). 
    - `03_filter_samples.R`: this script takes in the phyloseq objects created by script 02 and creates visualizations of the data read distribution and structure and then facilitates removing outliers and choosing a threshold for dropping low-read samples. It records which samples were dropped to the metadata file. 
        - This file takes in the following arguments, to be changed in the script: a name (used for file storage purposes), an amplicon ID (either 16S or ITS, for determining which dataset to format), and an edit metadata argument (Y or N, for determining whether the sequence count information should be written into the metadata table). 
