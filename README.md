# M-BUDS: MicroBiome Urban Dysbiosis in Street Trees 
Kathryn Atherton Dysbiosis in Urban Street Trees (Chapter 2)
The purpose of this project is to understand how urbanization impacts the tree-associated microbiome. 

## Requirements
To run this analysis you will need the following R packages (version I used for this analysis):
- compositions (2.0-8)
- decontam (1.20.0)
- dplyr (1.1.4)
- emmeans (1.11.1)
- fields (16.3.1)
- geosphere (1.5-18)
- ggplot2 (3.5.1)
- ggrepel (0.9.3)
- ggsignif (0.6.4)
- ggtext (0.1.2)
- gridExtra (2.3)
- multcompView (0.1-10)
- MuMIn (1.48.4)
- nlme (3.1-166)
- optparse (1.7.5)
- phyloseq (1.44.0)
- plyr (1.8.9)
- readxl (1.4.2)
- reshape2 (1.4.4)
- robCompositions (2.4.1)
- sva (3.48.0)
- tidyr (1.3.1)
- vegan (2.6-4)
- vroom (1.6.5)

Additionally, the initial sequence data cleaning step is run with the BU16S pipeline, created by Michael Silverstein and published in [Nature Ecology and Evolution](https://www.nature.com/articles/s41559-024-02440-6). More information about the pipeline can be found [here](https://github.com/Boston-University-Microbiome-Initiative/BU16s). 

## About the directories and scripts
- `00_functions.R`: contains the custom functions used repeatedly throughout the analyses.

- `01_Collect_Data`: this is where I store my raw data/metadata files. These are not currently shared on GitHub, and no scripts are stored here.
- `02_Clean_Data`: this is where I store my cleaned data/metadata files and the files that I use to make them.
    - `01_DADA2`: this is where I store my scripts to run DADA2 through the BU16S pipeline to perform QC and taxonomic annotation on my sequence data.
      - There is one directory per 16S/ITS sequencing run here. Each directory contains the qsub script required to run the BU16S pipeline. The pipeline outputs are also stored here, but they are not currently shared on GitHub.
    - `00_Scripts`: This holds all of the Data Cleaning scripts, as described below. Each process is both an R script and a bash (.sh) script that you can use to run the R script on a computing cluster. 
        - `02_format_raw_dada2_asv_tables.R`: this script takes the BU16S taxonomy and ASV table outputs and formats them for downstream cleaning steps. It merges all sequencing runs into one ASV table and then writes the post-DADA2 sequence count information into the metadata table. It also merges the taxonomy assignments into one table and formats the taxonomic information into separate columns for each taxonomic rank. Finally, it separates the ASV tables into one table for each sample type (leaf, root, upper soil horizon, lower soil horizon) and saves a csv and phyloseq object for each sample type for downstream analyis.
            - This file takes in the following arguments: an amplicon ID (either 16S or ITS, for determining which dataset to format), a name (used for file storage purposes), and an edit metadata argument (Y or N, for determining whether the sequence count information should be written into the metadata table). 
        - `03_filter_samples.R`: this script takes in the phyloseq objects created by script 02 and creates visualizations of the data read distribution and structure and then facilitates removing outliers and choosing a threshold for dropping low-read samples. It records which samples were dropped to the metadata file. 
            - This file takes in the following arguments, to be changed in the script: an amplicon ID (either 16S or ITS, for determining which dataset to format), a name (used for file storage purposes), and an edit metadata argument (Y or N, for determining whether the sequence count information should be written into the metadata table).
        - `04_decontaminate_samples.R`: this script uses the decontam R package to remove contaminants from the samples, as determined by ASV prevalence in negative controls.
            - This file takes in the following arguments, to be changed in the script: an amplicon ID (either 16S or ITS, for determining which dataset to format), a name (used for file storage purposes), and an edit metadata argument (Y or N, for determining whether the sequence count information should be written into the metadata table).
        - `05_transform_data.R`: this script is used to manipulate data a number of ways: 1) check if batch correction is needed, 2) run batch correction if needed and wanted, and 3) perform the rarefaction, CLR, Z-Score, and/or Aitchison Distance transformationscalculations.
            - This file takes in the following arguments, to be changed in the script: an amplicon ID (either 16S or ITS, for determining which dataset to format), a name (used for file storage purposes), an edit metadata argument (Y or N, for determining whether the sequence count information should be written into the metadata table), a check for batch effect argument (Y or N, for determining whether to check for a batch effect), a run batch correction argument (Y or N, for determining whether to run batch correction, if the check batch correction function determined that batch correction was recommended), a rarefy argument (Y or N, for determining whether to rarefy the data), a CLR-transform argument (Y or N, for determining whether to perform a Center Log Ratio transformation on the data), a Z-Score transform argument (Y or N, for determining whether to Z-Score transform the data), and an Aitchison Distance argument (Y or N, for determining whether to calculate the Aitchison Distance matrix for the data). 

- `03_Data_Analysis`: this is where I keep my files and scripts for calculating new metrics about data
    - `01_Calculate_Diversity_and_Guild_Abundances`: This is where I keep scripts for calculating 16S and ITS diversity metrics and functional guild abundances.
        - `atherton_16S_calculate_diversity_and_guild_abundance.Rmd`: This file takes in the bacterial metadata files, calculates the sample distance from Boston, as defined by the Park Street MBTA station; the percent relative abundance of functional guilds, as defined by an in-house database; microbial diversity metrics (Shannon, Simpson, Inverse Simpson, and Fisher); and matches the bacterial ASVs to pathogen status, as defined by blasting their representative sequences from DADA2 to the Multiple Bacterial Pathogens Database (MBPD).
        -  `atherton_ITS_calculate_diversity_and_guild_abundance.Rmd`: This file takes in the fungal metadata files, calculates the sample distance from Boston, as defined by the Park Street MBTA station; the percent relative abundance of functional guilds, as defined by the FungalTraits' and microbial diversity metrics (Shannon, Simpson, Inverse Simpson, and Fisher).
     - `02_Identify_Pathogens`: This is where I keep scripts for BLASTing the DADA2-derived representative sequences from 16S data against the MBPD database.
         - `01_blast_16s_pathogens.sh`: This script concatenates the DADA2-derived representative ASV sequences from multiple sequencing runs into one file and removes duplicated sequences. It then unzips the MBPD database and BLASTs against it. It then formats the output into a table that can be read into R.

- `04_Make_Figures`: This contains the folder for all of my figures included in the final paper. There is one R script per figure panel. The statistics for the analysis paired with the figure is also included in each figure script. 
