### LOAD IN LIBRARIES #########################################################
library(dplyr)
library(ggplot2)
library(nlme)
library(emmeans)
library(ggtext)
library(MuMIn)
library(reshape2)

### SET UP SCRIPT #############################################################
shapes <- c(16, 17, 15, 8)
linear_pal <- c("#669BBC", "#003049")
categorical_pal <- c("#FF8080","#FF0000","#C00000","#800000","#000000")

mean_se <- function(x, mult = 1) {  
  x <- na.omit(x)
  se <- mult * sqrt(var(x) / length(x))
  mean <- mean(x)
  data.frame(y = mean, ymin = mean - se, ymax = mean + se)
}
standard_error <- function(x) sd(x) / sqrt(length(x)) 

# read in metadata
setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/01_Collect_Data/01_Sample_Metadata/")
metadata_16s <- read.csv("atherton_sample_metadata_16S_rareto200_20250121.csv")
metadata_16s <- distinct(metadata_16s)
metadata_its <- read.csv("atherton_sample_metadata_ITS_20250125.csv")
metadata_its <- distinct(metadata_its)

# read in sample names of cleaned data
setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/02_Clean_Data/05_Transform_Data/ITS/Aitchison_Distance_Tables/")
all_aitch_its <- read.csv("atherton_ITS_allsampletypes_aitchisondistance_20240919.csv", 
                          row.names = 1)
setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/02_Clean_Data/05_Transform_Data/16S/Aitchison_Distance_Tables/")
bg_aitch_16s <- read.csv("atherton_16S_allsampletypes_aitchisondistance_20240925.csv", 
                         row.names = 1)
leaf_aitch_16s <- read.csv("atherton_16S_leaf_aitchisondistance_20241112.csv", 
                           row.names = 1)

### FORMAT METADATA ###########################################################
# Combine ITS soil horizon with sample type
metadata_its$sample_type<- paste0(metadata_its$soil_horizon, " ", 
                                  metadata_its$sample_type)
metadata_its$sample_type <- gsub("NA ", "", 
                                 metadata_its$sample_type)

# Replace ITS soil horizon category with depth measurement
metadata_its$sample_type <- gsub("O Soil", "Soil, 0-15 cm", 
                                 metadata_its$sample_type)
metadata_its$sample_type <- gsub("M Soil", "Soil, 15-30 cm", 
                                 metadata_its$sample_type)

metadata_its$sample_type <- gsub("O Root", "Roots, 0-15 cm", 
                                 metadata_its$sample_type)
metadata_its$sample_type <- gsub("M Root", "Roots, 15-30 cm", 
                                 metadata_its$sample_type)

# Replace ITS street tree separate sample types with just Street Tree
metadata_its$tree_pit_type <- gsub("Pit", "Street Tree",
                                   metadata_its$tree_pit_type)
metadata_its$tree_pit_type <- gsub("Grass", "Street Tree",
                                   metadata_its$tree_pit_type)

# Filter the sample names for just the samples from the cleaned ITS dataset
metadata_its <- metadata_its[which(metadata_its$sample_name %in% 
                                     colnames(all_aitch_its)),]

# Make tree pit type a factor with the levels order
metadata_its$tree_pit_type <- factor(metadata_its$tree_pit_type,
                                     levels = c("Street Tree", "Urban Edge",
                                                "Urban Interior", "Rural Edge",
                                                "Rural Interior"))

# Format the ITS environmental data
metadata_its$soil_percent_moisture <- gsub("%", "", 
                                           metadata_its$soil_percent_moisture)
metadata_its$soil_percent_organic_matter <- gsub("%", "",
                                                 metadata_its$soil_percent_organic_matter)
metadata_its$soil_percent_moisture <- as.numeric(metadata_its$soil_percent_moisture)
metadata_its$soil_percent_organic_matter <- as.numeric(metadata_its$soil_percent_organic_matter)
metadata_its$tree_distance_from_edge <- as.numeric(metadata_its$tree_distance_from_edge)

# Combine 16S soil horizon with sample type
metadata_16s$sample_type<- paste0(metadata_16s$soil_horizon, " ",
                                  metadata_16s$sample_type)
metadata_16s$sample_type <- gsub("NA ", "", 
                                 metadata_16s$sample_type)

# Replace 16S soil horizon category with depth measurement
metadata_16s$sample_type <- gsub("O Soil", "Soil, 0-15 cm",
                                 metadata_16s$sample_type)
metadata_16s$sample_type <- gsub("M Soil", "Soil, 15-30 cm",
                                 metadata_16s$sample_type)

metadata_16s$sample_type <- gsub("O Root", "Roots, 0-15 cm",
                                 metadata_16s$sample_type)
metadata_16s$sample_type <- gsub("M Root", "Roots, 15-30 cm",
                                 metadata_16s$sample_type)

# Replace 16S street tree separate sample types with just Street Tree
metadata_16s$tree_pit_type <- gsub("Pit", "Street Tree",
                                   metadata_16s$tree_pit_type)
metadata_16s$tree_pit_type <- gsub("Grass", "Street Tree",
                                   metadata_16s$tree_pit_type)

# Filter the sample names for just the samples from the cleaned ITS dataset
metadata_16s <- metadata_16s[which(metadata_16s$sample_name %in% 
                                     c(colnames(bg_aitch_16s), 
                                       colnames(leaf_aitch_16s))),]
metadata_16s <- metadata_16s[which(metadata_16s$sequences_dropped == "No"),]

# Make tree pit type a factor with the levels order
metadata_16s$tree_pit_type <- factor(metadata_16s$tree_pit_type,
                                     levels = c("Street Tree", "Urban Edge",
                                                "Urban Interior", 
                                                "Rural Edge",
                                                "Rural Interior"))

# Format the 16S environmental data
metadata_16s$soil_percent_moisture <- gsub("%", "",
                                           metadata_16s$soil_percent_moisture)
metadata_16s$soil_percent_organic_matter <- gsub("%", "",
                                                 metadata_16s$soil_percent_organic_matter)
metadata_16s$soil_percent_moisture <- as.numeric(metadata_16s$soil_percent_moisture)
metadata_16s$soil_percent_organic_matter <- as.numeric(metadata_16s$soil_percent_organic_matter)
metadata_16s$tree_distance_from_edge <- as.numeric(metadata_16s$tree_distance_from_edge)

# Normalize the data to be used in this figure
metadata_its$perc_dung_saprotroph <- log(metadata_its$perc_dung_saprotroph+1)
metadata_its$perc_wood_saprotroph <- log(metadata_its$perc_wood_saprotroph+1)
metadata_its$perc_litter_saprotroph <- log(metadata_its$perc_litter_saprotroph+1)
metadata_its$perc_soil_saprotroph <- log(metadata_its$perc_soil_saprotroph+1)

metadata_16s$perc_chitinolytic <- log(metadata_16s$perc_chitinolytic+1)
metadata_16s$perc_partial_nitrification <- log(metadata_16s$perc_partial_nitrification+1)
metadata_16s$perc_dissim_nitrate_reduction <- log(metadata_16s$perc_dissim_nitrate_reduction+1)
metadata_16s$perc_denitrification <- log(metadata_16s$perc_denitrification+1)
metadata_16s$perc_methanotroph <- log(metadata_16s$perc_methanotroph+1)

# Make a list of fungal guild names
fungal_trait_names <- c("perc_algal_parasite", "perc_animal_endosymbiont", 
                        "perc_animal_parasite", 
                        "perc_amf", "perc_arthropod_associated", 
                        "perc_dung_saprotroph", "perc_ecm", "perc_epiphyte",
                        "perc_foliar_endophyte", "perc_lichen_parasite",
                        "perc_lichenized", "perc_litter_saprotroph",
                        "perc_moss_symbiont", "perc_mycoparasite",
                        "perc_nectar_tap_saprotroph", 
                        "perc_plant_pathogen", 
                        "perc_pollen_saprotroph", "perc_root_endophyte",
                        "perc_soil_saprotroph", "perc_sooty_mold",
                        "perc_wood_saprotroph", 
                        "perc_animal_pathogen", 
                        "perc_human_pathogen", "perc_opportunistic_pathogen",
                        "perc_root_associated_pathogen", 
                        "perc_leaf_fruit_seed_pathogen", "perc_wood_pathogen",
                        "shannon", "shannon_nopath", "perc_all_saprotroph", 
                        "perc_all_pathotroph", "perc_all_symbiont")

# Make a list of bacterial guild names
bacterial_trait_names <- c("perc_c_fixation", "perc_cellulolytic", 
                           "perc_chitinolytic", "perc_lignolytic",
                           "perc_carbon_monoxide_oxidation", "perc_n_fixation",
                           "perc_copiotroph", "perc_denitrification",
                           "perc_dissim_nitrate_reduction", 
                           "perc_hydrocarbon_degradation",
                           "perc_sulfonate_desulfurization",
                           "perc_other_p_cycling",
                           "perc_oxidize_reduced_sulfur", "perc_methanotroph",
                           "perc_oligotroph", "perc_partial_nitrification", 
                           "perc_archaea", "perc_plant_pathogen", 
                           "perc_animal_pathogen", "perc_zoonotic_pathogen", 
                           "shannon", "shannon.x", "shannon_nopath")

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

### FIGURE S2A #################################################################
# Select soil, top layer data from 16S
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 52)
print("Number of trees used in the root analysis:")
length(unique(metadata_soil_16s$tree_id))

# Format the x-axis tick marks
metadata_soil_16s$site_edge <- paste0(metadata_soil_16s$tree_urban_type, " ", 
                                      metadata_soil_16s$tree_edge_type)
metadata_soil_16s$site_edge <- gsub("Urban Street Tree Urban Street Tree", "Urban Street Tree",
                                    metadata_soil_16s$site_edge)
metadata_soil_16s$site_edge <- factor(metadata_soil_16s$site_edge, 
                                      levels = c("Urban Street Tree", 
                                                 "Urban Forest Edge",
                                                 "Urban Forest Interior",
                                                 "Rural Forest Edge", 
                                                 "Rural Forest Interior"))

# Take the functional guild abundances part of the metadata from roots
vdata <- metadata_soil_16s[which(colnames(metadata_soil_16s) %in% 
                                          c("sample_name", "tree_id", 
                                            "plot_name", "tree_edge_type", 
                                            "tree_urban_type", 
                                            "soil_percent_organic_matter", 
                                            "site_edge"))]

# Rename columns
colnames(vdata) <- c("sample_name", "tree_id", "plot_name", 
                            "tree_edge_type", "tree_urban_type", 
                            "soil_percent_organic_matter", "site_edge")

### STATISTICS SOIL ORGANIC MATTER ###
model_som <- lme(formula(soil_percent_organic_matter ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_16s)
summary(model_som)
summary(emmeans(model_som, pairwise ~ site_edge))

# Make the text for adding Tukey groups and record the location of the group id on the y-axis
vdata$tukey <- NA
vdata$location <- NA

for(i in 1:nrow(vdata)){
  site_edge <- as.character(vdata$site_edge[i])
  
    if(site_edge == "Urban Street Tree"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$soil_percent_organic_matter[which(vdata$site_edge == site_edge)]) + 0.02
    } else if(site_edge == "Urban Forest Edge"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$soil_percent_organic_matter[which(vdata$site_edge == site_edge)]) + 0.02
    } else if(site_edge == "Urban Forest Interior"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$soil_percent_organic_matter[which(vdata$site_edge == site_edge)]) + 0.02
    } else if(site_edge == "Rural Forest Edge"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$soil_percent_organic_matter[which(vdata$site_edge == site_edge)]) + 0.02
    } else if(site_edge == "Rural Forest Interior"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$soil_percent_organic_matter[which(vdata$site_edge == site_edge)]) + 0.02
    }
}

# Format the x-axis tick marks for plotting
vdata$site_edge <- gsub(" ", "\n", as.character(vdata$site_edge))
vdata$site_edge <- factor(vdata$site_edge, 
                                      levels = c("Urban\nStreet\nTree", 
                                                 "Urban\nForest\nEdge",
                                                 "Urban\nForest\nInterior",
                                                 "Rural\nForest\nEdge",
                                                 "Rural\nForest\nInterior"))


# Visualize
figs2a <- ggplot() +
  stat_summary(data = vdata, aes(x = site_edge, y = soil_percent_organic_matter, 
                                 group = 1),
               fun = median,
               geom = "line",
               linetype = "dashed",
               linewidth = 0.8) +
  geom_boxplot(data = vdata, outlier.shape = NA,
               aes(x=site_edge, y=soil_percent_organic_matter, fill=site_edge), 
               alpha = 0.5) + 
  geom_jitter(data=vdata, aes(x = site_edge, y = soil_percent_organic_matter, 
                                          color = site_edge, alpha = 0.75, 
                                          stroke = NA)) +
  theme_classic() + labs(y="Soil organic matter content (%)", x="", shape = "") +
  # this adds the Tukey groups labels defined above
  geom_text(data=vdata, family = "Arial",
            aes(x = site_edge, y = location, label = tukey), size = 3, 
            parse = TRUE) +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_markdown(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), 
        legend.position = "none", 
        legend.background = element_rect(fill = NA)) +
  scale_color_manual(values = categorical_pal, guide = "none") + 
  scale_fill_manual(values = categorical_pal, guide = "none") +
  ggtitle("Soil organic matter content")

figs2a
ggsave(paste0("atherton_figures2a_organicmattersoil_vs_treetype_boxplot_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs2a, dpi = 300, width = 6, height = 3.5, units = "in")

rm(list = c("metadata_soil_16s", "model_som", "vdata"))

### FIGURE S2B #################################################################
# Select the soil, top layer data from 16S
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 52)
print("Number of trees used in the root analysis:")
length(unique(metadata_soil_16s$tree_id))

# Format the x-axis tick marks
metadata_soil_16s$site_edge <- paste0(metadata_soil_16s$tree_urban_type, " ", 
                                      metadata_soil_16s$tree_edge_type)
metadata_soil_16s$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                                    "Urban Street Tree",
                                    metadata_soil_16s$site_edge)
metadata_soil_16s$site_edge <- factor(metadata_soil_16s$site_edge, 
                                      levels = c("Urban Street Tree", 
                                                 "Urban Forest Edge", 
                                                 "Urban Forest Interior", 
                                                 "Rural Forest Edge", 
                                                 "Rural Forest Interior"))

# Keep the inorganic N deposition data
data <- cbind(metadata_soil_16s$rindy_2019_nh4n_deposition.kgn_ha_yr, 
              metadata_soil_16s$rindy_2019_no3n_deposition.kgn_ha_yr)
colnames(data) <- c("NH<sub>4</sub><sup>+</sup> annual throughfall", 
                    "NO<sub>3</sub><sup>-</sup> annual throughfall")

# Take the mean and standard error of the data
mean <- aggregate(data~site_edge, data = metadata_soil_16s, FUN= "mean",
                  na.action = na.omit)
se <- aggregate(data~site_edge, data = metadata_soil_16s, FUN= standard_error,
                na.action = na.omit)

# Reformat the data for plotting
melt.mean <- melt(mean)
melt.se <- melt(se)

# Combine the mean and standard error into one data frame
vdata <- cbind(melt.mean, melt.se[,ncol(melt.se)])
colnames(vdata) <- c("site_edge","Nutrient","Abundance","SE")

# Reformat the x-axis tick marks for plotting
vdata$site_edge <- gsub(" ", "\n", vdata$site_edge)
vdata$site_edge <- factor(vdata$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

# Reorder the nutrients for plotting
vdata$Nutrient <- factor(vdata$Nutrient, 
                         levels = c("NH<sub>4</sub><sup>+</sup> annual throughfall", 
                                    "NO<sub>3</sub><sup>-</sup> annual throughfall"))

# Plot figure s2b
figs2b <- ggplot() +
  # This adds the dotted line tracking average abundance across groups
  geom_line(data = vdata, 
            aes(x=site_edge, y=Abundance, shape = Nutrient, group = Nutrient), 
            linetype = 2) +
  # This adds the average abundance point
  geom_point(data = vdata, 
             aes(x=site_edge, y=Abundance, color=site_edge, shape = Nutrient, 
                 group = Nutrient), size=5) + 
  # This adds the error bars around the average abundance point
  geom_errorbar(data = vdata, 
                aes(x=site_edge, ymin=Abundance-SE, ymax=Abundance+SE, 
                    color=site_edge), width = 0.2) +
  theme_classic() + labs (y="Atmospheric inorganic N deposition (kg/ha)", x="") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_markdown(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_markdown(size = 10), 
        legend.position = c(0.85, 1), 
        legend.background = element_rect(fill = NA)) +
  scale_color_manual(values = categorical_pal) + 
  scale_fill_manual(values = categorical_pal) +
  guides(color = "none", shape = guide_legend(title = ""))+
  ggtitle("Annual atmospheric N throughfall")

figs2b
ggsave(paste0("atherton_figures2b_nthroughfall_vs_treetype_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figs2b, dpi = 300, width = 6, height = 3.5, units = "in")

rm(list = c("data", "mean", "melt.mean", "melt.se", "metadata_soil_16s", "se", 
            "vdata"))

### FIGURE S2C #################################################################
# Select the soil, top layer data from ITS
metadata_soil_its <- metadata_its[which(metadata_its$sample_type == "Soil, 0-15 cm"),]

# Select the soil, top layer data from 16S
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 0-15 cm"),]

# Select the leaf data from 16S
metadata_leaf_16s <- metadata_16s[which(metadata_16s$sample_type == "Leaf"),]

# Keep the functional group data from the metadata
func_soil_its <- metadata_soil_its[,colnames(metadata_soil_its) %in%
                                             c(fungal_trait_names)]
func_soil_16s <- metadata_soil_16s[,colnames(metadata_soil_16s) %in% 
                                                c(bacterial_trait_names)]

func_leaf_16s <- metadata_leaf_16s[,colnames(metadata_leaf_16s) %in% 
                                                c(bacterial_trait_names)]

# Format the soil ITS column names for plotting
func_soil_its <- as.data.frame(func_soil_its[,which(colnames(func_soil_its) %in%
                                                      c("perc_dung_saprotroph", 
                                                        "perc_wood_saprotroph",
                                                        "perc_litter_saprotroph", 
                                                        "perc_soil_saprotroph"))])

colnames(func_soil_its) <- c("Fungal dung decomposers", 
                             "Fungal litter decomposers",
                             "Fungal soil decomposers", 
                             "Fungal wood decomposers")
colnames(func_soil_its) <- paste0(colnames(func_soil_its), " in soil")

# Format the soil 16S column names for plotting
func_soil_16s <- as.data.frame(func_soil_16s[,which(colnames(func_soil_16s) %in% 
                                                      c("perc_chitinolytic", 
                                                        "perc_n_fixation",
                                                        "perc_dissim_nitrate_reduction", 
                                                        "perc_denitrification"))])
colnames(func_soil_16s) <- c("Chitinolytic bacteria", 
                             "Nitrifying bacteria",
                             "Denitrifying bacteria",
                             "Dissimilatory nitrate reducing bacteria")
colnames(func_soil_16s) <- paste0(colnames(func_soil_16s), " in soil")

# Format the leaf 16S column names for plotting
func_leaf_16s <- as.data.frame(func_leaf_16s[,which(colnames(func_leaf_16s) %in% 
                                                      c("perc_methanotroph"))])
colnames(func_leaf_16s) <- c("Methanotrophic bacteria")
colnames(func_leaf_16s) <- paste0(colnames(func_leaf_16s), " in leaves")

# Keep the environmental data
metadata_soil_heatmap_its <- metadata_soil_its[,c(17,18,20,23,30:36,38,49,53:55,59:61,63:65,67:76)]
metadata_soil_heatmap_16s <- metadata_soil_16s[,c(17,18,20,23,30:36,38,49,53:55,59:61,63:65,67:76)]
metadata_leaf_heatmap_16s <- metadata_leaf_16s[,c(17,18,20,23,30:36,38,49,53:55,59:61,63:65,67:76)]

# Rename the columns
colnames(metadata_soil_heatmap_its) <- c("Tree distance from edge", 
                                         "Tree distance from Boston", 
                                         "Tree DBH", "Tree growth rate",
                                         "Soil temperature", 
                                         "Soil bulk density", 
                                         "Soil pH", "Soil moisture", 
                                         "Soil organic matter", 
                                         "Soil available NH<sub>4</sub><sup>+</sup>",
                                         "Soil available NO<sub>3</sub><sup>-</sup>", 
                                         "Normalized root biomass", 
                                         "Soil soluble salt",
                                         "NO concentration",
                                         "NO<sub>2</sub> concentration",
                                         "NO<sub>x</sub> concentration",
                                         "NH<sub>4</sub><sup>+</sup> annual throughfall",
                                         "NO<sub>3</sub><sup>-</sup> annual throughfall",
                                         "Total inorganic N annual throughfall",
                                         "O<sub>3</sub> concentration",
                                         "Foliar N %",
                                         "Foliar C:N",
                                         "Ca<sup>2+</sup> annual wet deposition",
                                         "Cl<sup>-</sup> annual wet deposition",
                                         "H<sup>+</sup> annual wet deposition",
                                         "K<sup>+</sup> annual wet deposition",
                                         "Mg<sup>2+</sup> annual wet deposition",
                                         "Na<sup>+</sup> annual wet deposition",
                                         "NH<sub>4</sub><sup>+</sup> annual wet deposition",
                                         "NO<sub>3</sub><sup>-</sup> annual wet deposition",
                                         "SO<sub>4</sub><sup>2-</sup> annual wet deposition",
                                         "Total inorganic N annual wet deposition")

# Rename the columns
colnames(metadata_soil_heatmap_16s) <- c("Tree distance from edge", 
                                         "Tree distance from Boston", 
                                         "Tree DBH", "Tree growth rate",
                                         "Soil temperature", 
                                         "Soil bulk density", 
                                         "Soil pH", "Soil moisture", 
                                         "Soil organic matter", 
                                         "Soil available NH<sub>4</sub><sup>+</sup>",
                                         "Soil available NO<sub>3</sub><sup>-</sup>", 
                                         "Normalized root biomass", 
                                         "Soil soluble salt",
                                         "NO concentration",
                                         "NO<sub>2</sub> concentration",
                                         "NO<sub>x</sub> concentration",
                                         "NH<sub>4</sub><sup>+</sup> annual throughfall",
                                         "NO<sub>3</sub><sup>-</sup> annual throughfall",
                                         "Total inorganic N annual throughfall",
                                         "O<sub>3</sub> concentration",
                                         "Foliar N %",
                                         "Foliar C:N",
                                         "Ca<sup>2+</sup> annual wet deposition",
                                         "Cl<sup>-</sup> annual wet deposition",
                                         "H<sup>+</sup> annual wet deposition",
                                         "K<sup>+</sup> annual wet deposition",
                                         "Mg<sup>2+</sup> annual wet deposition",
                                         "Na<sup>+</sup> annual wet deposition",
                                         "NH<sub>4</sub><sup>+</sup> annual wet deposition",
                                         "NO<sub>3</sub><sup>-</sup> annual wet deposition",
                                         "SO<sub>4</sub><sup>2-</sup> annual wet deposition",
                                         "Total inorganic N annual wet deposition")

# Rename columns
colnames(metadata_leaf_heatmap_16s) <- c("Tree distance from edge", 
                                         "Tree distance from Boston", 
                                         "Tree DBH", "Tree growth rate",
                                         "Soil temperature", 
                                         "Soil bulk density", 
                                         "Soil pH", "Soil moisture", 
                                         "Soil organic matter", 
                                         "Soil available NH<sub>4</sub><sup>+</sup>",
                                         "Soil available NO<sub>3</sub><sup>-</sup>", 
                                         "Normalized root biomass", 
                                         "Soil soluble salt",
                                         "NO concentration",
                                         "NO<sub>2</sub> concentration",
                                         "NO<sub>x</sub> concentration",
                                         "NH<sub>4</sub><sup>+</sup> annual throughfall",
                                         "NO<sub>3</sub><sup>-</sup> annual throughfall",
                                         "Total inorganic N annual throughfall",
                                         "O<sub>3</sub> concentration",
                                         "Foliar N %",
                                         "Foliar C:N",
                                         "Ca<sup>2+</sup> annual wet deposition",
                                         "Cl<sup>-</sup> annual wet deposition",
                                         "H<sup>+</sup> annual wet deposition",
                                         "K<sup>+</sup> annual wet deposition",
                                         "Mg<sup>2+</sup> annual wet deposition",
                                         "Na<sup>+</sup> annual wet deposition",
                                         "NH<sub>4</sub><sup>+</sup> annual wet deposition",
                                         "NO<sub>3</sub><sup>-</sup> annual wet deposition",
                                         "SO<sub>4</sub><sup>2-</sup> annual wet deposition",
                                         "Total inorganic N annual wet deposition")

# Make coefficient matrix
# Soil ITS
Result_soil_its <- NULL
for (i in 1:ncol(func_soil_its)){
  estimate_soil_its <- NULL
  for (j in 1:ncol(metadata_soil_heatmap_its)){
    cor_soil_its <- cor.test(func_soil_its[,i],metadata_soil_heatmap_its[,j], method="p")
    estimate_soil_its <- c(estimate_soil_its, cor_soil_its$estimate)
  }
  Result_soil_its <- rbind(Result_soil_its, estimate_soil_its)
}
colnames(Result_soil_its) <- colnames(metadata_soil_heatmap_its)
rownames(Result_soil_its) <- colnames(func_soil_its)

# Soil 16S
Result_soil_16s <- NULL
for (i in 1:ncol(func_soil_16s)){
  estimate_soil_16s <- NULL
  for (j in 1:ncol(metadata_soil_heatmap_16s)){
    cor_soil_16s <- cor.test(func_soil_16s[,i],metadata_soil_heatmap_16s[,j], method="p")
    estimate_soil_16s <- c(estimate_soil_16s, cor_soil_16s$estimate)
  }
  Result_soil_16s <- rbind(Result_soil_16s, estimate_soil_16s)
}
colnames(Result_soil_16s) <- colnames(metadata_soil_heatmap_16s)
rownames(Result_soil_16s) <- colnames(func_soil_16s)

# Leaf ITS
Result_leaf_16s <- NULL
estimate_leaf_16s <- NULL
for (j in 1:ncol(metadata_leaf_heatmap_16s)){
  cor_leaf_16s <- cor.test(func_leaf_16s[,1],metadata_leaf_heatmap_16s[,j], method="p")
  estimate_leaf_16s <- c(estimate_leaf_16s, cor_leaf_16s$estimate)
}
Result_leaf_16s <- rbind(Result_leaf_16s, estimate_leaf_16s)
colnames(Result_leaf_16s) <- colnames(metadata_leaf_heatmap_16s)
rownames(Result_leaf_16s) <- "Methanotrophic bacteria in leaves" 

# Make Pval matrix
# Soil ITS
Result.Pval_soil_its <- NULL
for (i in 1:ncol(func_soil_its)){
  Pval_soil_its <- NULL
  for (j in 1:ncol(metadata_soil_heatmap_its)){
    cor_soil_its <- cor.test(func_soil_its[,i],metadata_soil_heatmap_its[,j], method="p")
    Pval_soil_its <- c(Pval_soil_its, cor_soil_its$p.value)
  }
  Result.Pval_soil_its <- rbind (Result.Pval_soil_its, Pval_soil_its)
}
colnames(Result.Pval_soil_its) <- colnames(metadata_soil_heatmap_its)
rownames(Result.Pval_soil_its) <- colnames(func_soil_its)

# Soil 16S
Result.Pval_soil_16s <- NULL
for (i in 1:ncol(func_soil_16s)){
  Pval_soil_16s <- NULL
  for (j in 1:ncol(metadata_soil_heatmap_16s)){
    cor_soil_16s <- cor.test(func_soil_16s[,i],metadata_soil_heatmap_16s[,j], method="p")
    Pval_soil_16s <- c(Pval_soil_16s, cor_soil_16s$p.value)
  }
  Result.Pval_soil_16s <- rbind (Result.Pval_soil_16s, Pval_soil_16s)
}
colnames(Result.Pval_soil_16s) <- colnames(metadata_soil_heatmap_16s)
rownames(Result.Pval_soil_16s) <- colnames(func_soil_16s)

# Leaf 16S
Result.Pval_leaf_16s <- NULL
Pval_leaf_16s <- NULL
for (j in 1:ncol(metadata_leaf_heatmap_16s)){
  cor_leaf_16s <- cor.test(func_leaf_16s[,1],metadata_leaf_heatmap_16s[,j], method="p")
  Pval_leaf_16s <- c(Pval_leaf_16s, cor_leaf_16s$p.value)
}
Result.Pval_leaf_16s <- rbind (Result.Pval_leaf_16s, Pval_leaf_16s)
colnames(Result.Pval_leaf_16s) <- colnames(metadata_leaf_heatmap_16s)
rownames(Result.Pval_leaf_16s) <- "Methanotrophic bacteria in leaves"

# Make Asterisk matrix
# Soil ITS
Result.Pval.combined_soil_its <- Result.Pval_soil_its
Result.asterisk_soil_its <- Result.Pval.combined_soil_its
for (k in 1:ncol(Result.Pval.combined_soil_its)){
  for (i in 1:nrow(Result.Pval.combined_soil_its)){
    if(!is.na(Result.Pval.combined_soil_its[i,k])){
      if (Result.Pval.combined_soil_its[i,k] < 0.001/(9*32)){
        Result.asterisk_soil_its[i,k] <- "***"
      } else if (Result.Pval.combined_soil_its[i,k] < 0.01/(9*32)){
        Result.asterisk_soil_its[i,k] <- "**"
      } else if (Result.Pval.combined_soil_its[i,k] < 0.05/(9*32)){
        Result.asterisk_soil_its[i,k] <- "*"
      } else {
        Result.asterisk_soil_its[i,k]<- " "
      }
    } 
  }
}

# Soil 16S
Result.Pval.combined_soil_16s <- Result.Pval_soil_16s
Result.asterisk_soil_16s <- Result.Pval.combined_soil_16s
for (k in 1:ncol(Result.Pval.combined_soil_16s)){
  for (i in 1:nrow(Result.Pval.combined_soil_16s)){
    if(!is.na(Result.Pval.combined_soil_16s[i,k])){
      if (Result.Pval.combined_soil_16s[i,k] < 0.001/(9*32)){
        Result.asterisk_soil_16s[i,k] <- "***"
      } else if (Result.Pval.combined_soil_16s[i,k] < 0.01/(9*32)){
        Result.asterisk_soil_16s[i,k] <- "**"
      } else if (Result.Pval.combined_soil_16s[i,k] < 0.05/(9*32)){
        Result.asterisk_soil_16s[i,k] <- "*"
      } else {
        Result.asterisk_soil_16s[i,k]<- " "
      }
    } 
  }
}

# Leaf 16S
Result.Pval.combined_leaf_16s <- Result.Pval_leaf_16s
Result.asterisk_leaf_16s <- Result.Pval.combined_leaf_16s
for (k in 1:ncol(Result.Pval.combined_leaf_16s)){
  for (i in 1:nrow(Result.Pval.combined_leaf_16s)){
    if(!is.na(Result.Pval.combined_leaf_16s[i,k])){
      if (Result.Pval.combined_leaf_16s[i,k] < 0.001/(9*32)){
        Result.asterisk_leaf_16s[i,k] <- "***"
      } else if (Result.Pval.combined_leaf_16s[i,k] < 0.01/(9*32)){
        Result.asterisk_leaf_16s[i,k] <- "**"
      } else if (Result.Pval.combined_leaf_16s[i,k] < 0.05/(9*32)){
        Result.asterisk_leaf_16s[i,k] <- "*"
      } else {
        Result.asterisk_leaf_16s[i,k]<- " "
      }
    } 
  }
}

# Resutls with p-values and asterisks
Result.combined_soil_its <- Result_soil_its
melt_soil_its <- melt(Result.combined_soil_its)
melt.asterisk_soil_its <- melt(Result.asterisk_soil_its)

Result.combined_soil_16s <- Result_soil_16s
melt_soil_16s <- melt(Result.combined_soil_16s)
melt.asterisk_soil_16s <- melt(Result.asterisk_soil_16s)

Result.combined_leaf_16s <- Result_leaf_16s
melt_leaf_16s <- melt(Result.combined_leaf_16s)
melt.asterisk_leaf_16s <- melt(Result.asterisk_leaf_16s)

# Combine all together
melt_all <- rbind(melt_soil_16s, melt_soil_its, melt_leaf_16s)
melt_asterisk_all <- rbind(melt.asterisk_soil_16s, melt.asterisk_soil_its, 
                           melt.asterisk_leaf_16s)

# Visualize
vdata_its <- cbind(melt_all, melt_asterisk_all$value)
colnames(vdata_its) <- c("microbe", "environment", "coef", "Pval")

# Order the environmental factors
order_of_env_factors <- vdata_its[which(vdata_its$microbe == "Methanotrophic bacteria in leaves"),]

vdata_its$environment <- factor(vdata_its$environment, 
                                levels = order_of_env_factors$environment[order(order_of_env_factors$coef)])

# Order the microbes
vdata_its$microbe <- factor(vdata_its$microbe, levels = c("Chitinolytic bacteria in soil",
                                                          "Dissimilatory nitrate reducing bacteria in soil",
                                                          "Denitrifying bacteria in soil",
                                                          "Fungal dung decomposers in soil",
                                                          "Fungal wood decomposers in soil",
                                                          "Fungal litter decomposers in soil",
                                                          "Nitrifying bacteria in soil",
                                                          "Fungal soil decomposers in soil",
                                                          "Methanotrophic bacteria in leaves"))

# Plot figure s2c
figs2c <- ggplot(vdata_its) +
  geom_tile(aes(y=microbe, x=environment, fill=coef)) +
  geom_text(aes(y=microbe, x=environment, label = Pval), size = 3) +
  labs(x="",y="") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.3)) +
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.3)) + 
  scale_fill_gradient2(low = "#003049", mid = "#FDF0D5", high = "#C00000", 
                       name = "Correlation\ncoefficient") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_markdown(size = 8, color = "black"), 
        axis.title.y = element_markdown(size = 8, color = "black"), 
        axis.text.y = element_markdown(size = 8, color = "black"),
        axis.text.x = element_markdown(size = 8, color = "black"),
        legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8), 
        plot.title = element_text(size = 8)) + 
  ggtitle("") + coord_fixed()

figs2c
ggsave(paste0("atherton_figures2c_functionalguilds_vs_environment_heatmap_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), plot = figs2c, 
       device = cairo_pdf, width = 5, height = 10, units = "in", dpi = 300)
