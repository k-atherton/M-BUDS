### LOAD IN LIBRARIES #########################################################
library(dplyr)
library(ggplot2)
library(nlme)
library(emmeans)
library(ggtext)
library(MuMIn)
library(reshape2)
library(ggsignif)
library(ggrepel)

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

### FIGURE 3A #################################################################
# Make new x-axis variable: combination of urban and edge types
metadata_16s$site_edge <- paste0(metadata_16s$tree_urban_type, " ", 
                                 metadata_16s$tree_edge_type)
metadata_16s$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                               "Urban Street Tree", metadata_16s$site_edge)

# Make variable a factor and order it for the x-axis
metadata_16s$site_edge <- factor(metadata_16s$site_edge, 
                                 levels = c("Urban Street Tree", 
                                            "Urban Forest Edge", 
                                            "Urban Forest Interior", 
                                            "Rural Forest Edge", 
                                            "Rural Forest Interior"))

# Select the data for only the soil top layer
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 88)
print("Number of trees used in the soil analysis:")
length(unique(metadata_soil_16s$tree_id))

# Create statistical data frames
Results_anova <- NULL
Results_multcomp <- NULL

# Keep the data that will be plotted in this figure
data <- as.data.frame(metadata_soil_16s$soil_no3)

# Rename the column as it will appear in the figure
colnames(data) <- c("Available soil NO<sub>3</sub><sup>-</sup>")

# Take the mean and standard error of the data for each tree type
mean <- aggregate(data$`Available soil NO<sub>3</sub><sup>-</sup>`~metadata_soil_16s$site_edge, 
                  FUN= "mean", na.action = na.omit)
se <- aggregate(data$`Available soil NO<sub>3</sub><sup>-</sup>`~metadata_soil_16s$site_edge, 
                FUN= standard_error, na.action = na.omit)

# Reformat the data for the figure plotting
melt.mean <- melt(mean)
melt.se <- melt(se)

# Combine the mean and standard error for plotting
vdata <- cbind(melt.mean, melt.se[,ncol(melt.se)])

# Rename the columns
colnames(vdata) <- c("site_edge","Nutirent","Abundance","SE")

### STATISTICS SOIL NO3- ###
model_no3 <- lme(formula(soil_no3 ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_16s)
summary(model_no3)
summary(emmeans(model_no3, pairwise ~ site_edge))

# Create the Tukey groups
vdata$tukey <- c(
  "a", 
  "ab", 
  "b", 
  "ab", 
  "b"
)

# Reformat the x-axis tickmark labels
metadata_soil_16s$site_edge <- gsub(" ", "\n", metadata_soil_16s$site_edge)
metadata_soil_16s$site_edge <- factor(metadata_soil_16s$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

vdata$site_edge <- gsub(" ", "\n", vdata$site_edge)
vdata$site_edge <- factor(vdata$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

# Visualize

fig3a <- ggplot() +
  # this creates the dotted line between points to track the average across tree types
  geom_line(data = vdata, 
            aes(x=site_edge, y=Abundance, group=Nutirent), 
            linetype = 2, color = "black") +
  # this creates the points that show the average abundance for each guild and tree type
  geom_point(data = vdata, 
             aes(x=site_edge, y=Abundance, color=site_edge, shape=Nutirent), 
             size=5) + 
  # this creates the standard error bars around the average abundance points
  geom_errorbar(data = vdata, 
                aes(x=site_edge, ymin=Abundance-SE, ymax=Abundance+SE, 
                    color=site_edge), width = 0.2) +
  # this adds the Tukey groups labels defined above
  geom_text(data = vdata, direction = "y",
            aes(x=site_edge, y=Abundance + SE + 0.02, label = tukey), 
            size = 3, parse = TRUE) +
  # this adds the individual data points to the figure
  geom_point(data = metadata_soil_16s,
             aes(x = site_edge, y = soil_no3, color = site_edge),
             position = position_jitter(width = 0.2), size = 2, alpha = 0.2, 
             shape = 16, stroke = NA) +
  theme_classic() + labs (y="Soil inorganic NO<sub>3</sub><sup>-</sup> content (ug/g)", 
                          x="") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_markdown(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_markdown(size = 10, color = "black"), 
        legend.text = element_markdown(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_markdown(size = 10), legend.position = "none", 
        legend.background = element_rect(fill = NA)) +
  scale_color_manual(values = categorical_pal) + 
  guides(color = "none", shape = guide_legend(title = ""))+
  ggtitle("Available soil NO<sub>3</sub><sup>-</sup>")

fig3a
ggsave(paste0("atherton_figure3a_no3soil_vs_treetype_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig3a, dpi = 300, width = 6, height = 3.5, units = "in")

rm(list = c("data", "mean", "melt.mean", "melt.se", "metadata_soil_16s",
            "model_no3", "se", "vdata"))

### FIGURE 3B ##################################################################
# Select the soil, top layer data from 16S
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 88)
print("Number of trees used in the soil analysis:")
length(unique(metadata_soil_16s$tree_id))

# Create the x-axis tick marks
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

# Transform the data for the figure
data <- cbind(metadata_soil_16s$perc_chitinolytic, 
              metadata_soil_16s$perc_dissim_nitrate_reduction,
              metadata_soil_16s$perc_partial_nitrification, 
              metadata_soil_16s$perc_denitrification)

# Rename the columns
colnames(data) <- c("Chitinolytic bacteria", 
                    "Dissimilatory nitrate reducing bacteria",
                    "Nitrifying bacteria", "Denitrifying bacteria")

# Take the mean and standard error of the data
mean <- aggregate(data~site_edge, data = metadata_soil_16s, FUN= "mean", 
                  na.action = na.omit)
se <- aggregate(data~site_edge, data = metadata_soil_16s, FUN= standard_error, 
                na.action = na.omit)

# Format the data frame
melt.mean <- melt(mean)
melt.se <- melt(se)

# Combine the mean and standard error data frames
vdata <- cbind(melt.mean, melt.se[,ncol(melt.se)])

# Rename the columns
colnames(vdata) <- c("site_edge","Microbe","Abundance","SE")

# Order the microbial groups for the plot
vdata$Microbe <- factor(vdata$Microbe, levels = c("Chitinolytic bacteria", 
                                                  "Nitrifying bacteria",
                                                  "Dissimilatory nitrate reducing bacteria",
                                                  "Denitrifying bacteria"))

model_chitin <- lme(formula(perc_chitinolytic ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_16s)
summary(model_chitin)
summary(emmeans(model_chitin, pairwise ~ site_edge))

model_dnr <- lme(formula(perc_dissim_nitrate_reduction ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_16s)
summary(model_dnr)
summary(emmeans(model_dnr, pairwise ~ site_edge))

model_nit <- lme(formula(perc_partial_nitrification ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_16s)
summary(model_nit)
summary(emmeans(model_nit, pairwise ~ site_edge))

model_denit <- lme(formula(perc_denitrification ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_16s)
summary(model_denit)
summary(emmeans(model_denit, pairwise ~ site_edge))

# Create the Tukey groups text
vdata$tukey <- c(
  "a[\"\\u25A0\"]", # chitinolytic street tree
  "b[\"\\u25A0\"]", # chitinolytic urban forest edge
  "b[\"\\u25A0\"]", # chitinolytic urban forest interior
  "b[\"\\u25A0\"]", # chitinolytic rural forest edge
  "b[\"\\u25A0\"]", # chitinolytic rural forest interior
  "a[\"\\u25B2\"]", # dissim street tree
  "b[\"\\u25B2\"]", # dissim urban forest edge
  "b[\"\\u25B2\"]", # dissim urban forest interior
  "b[\"\\u25B2\"]", # dissim rural forest edge
  "b[\"\\u25B2\"]", # dissim rural forest interior
  "a[\"\\u25CF\"]", # nitrifying street tree
  "b[\"\\u25CF\"]", # nitrifying urban forest edge
  "b[\"\\u25CF\"]", # nitrifying urban forest interior
  "b[\"\\u25CF\"]", # nitrifying rural forest edge
  "b[\"\\u25CF\"]", # nitrifying rural forest interior
  "a[\"\\u25C6\"]", # denitrifying street tree
  "b[\"\\u25C6\"]", # denitrifying urban forest edge
  "b[\"\\u25C6\"]", # denitrifying urban forest interior
  "b[\"\\u25C6\"]", # denitrifying rural forest edge
  "b[\"\\u25C6\"]" # denitrifying rural forest interior
)

# Reformat the x-axis tick marks for the plot
vdata$site_edge <- gsub(" ", "\n", vdata$site_edge)
vdata$site_edge <- factor(vdata$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

metadata_soil_16s$site_edge <- gsub(" ", "\n", metadata_soil_16s$site_edge)
metadata_soil_16s$site_edge <- factor(metadata_soil_16s$site_edge, 
                                      levels = c("Urban\nStreet\nTree",
                                                 "Urban\nForest\nEdge",
                                                 "Urban\nForest\nInterior",
                                                 "Rural\nForest\nEdge",
                                                 "Rural\nForest\nInterior"))

# Plot figure 3b
fig3b <- ggplot(vdata) +
  # this creates the dotted line between points to track the average across tree types
  geom_line(data = vdata, 
            aes(x=site_edge, y=Abundance, group=Microbe), 
            linetype = 2, color = "black") +
  # this creates the points that show the average abundance for each guild and tree type
  geom_point(data = vdata, 
             aes(x=site_edge, y=Abundance, color=site_edge, shape=Microbe), 
             size=5) + 
  # this creates the standard error bars around the average abundance points
  geom_errorbar(data = vdata, 
                aes(x=site_edge, ymin=Abundance-SE, ymax=Abundance+SE, 
                    color=site_edge), width = 0.2) +
  # this adds the Tukey groups labels defined above
  geom_text_repel(data = vdata, direction = "y",
                  aes(x=site_edge, y=Abundance, label = tukey), 
                  hjust = 1, segment.color = NA,
                  size = 3, parse = TRUE, nudge_x = -0.1) +
  # this adds the individual data points to the figure
  geom_point(data = metadata_soil_16s,
             aes(x = site_edge, y = perc_chitinolytic, color = site_edge),
             position = position_jitter(width = 0.2), size = 3, alpha = 0.2,
             shape = 15, stroke = NA) +
  geom_point(data = metadata_soil_16s,
             aes(x = site_edge, y = perc_partial_nitrification, 
                 color = site_edge), position = position_jitter(width = 0.2), 
             size = 3, alpha = 0.2, shape = 1) +
  geom_point(data = metadata_soil_16s,
             aes(x = site_edge, y = perc_dissim_nitrate_reduction, 
                 color = site_edge), position = position_jitter(width = 0.2), 
             size = 3, alpha = 0.2, shape = 17, stroke = NA) +
  geom_point(data = metadata_soil_16s,
             aes(x = site_edge, y = perc_denitrification, color = site_edge),
             position = position_jitter(width = 0.2), size = 3, alpha = 0.2,
             shape = 5) +
  theme_classic() + labs (y="log(Functional group relative abundance (%))", x="") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_markdown(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), legend.position = c(0.6, 0.8), 
        legend.background = element_rect(fill = NA)) +
  scale_color_manual(values = categorical_pal) + 
  scale_shape_manual(values = c(15,16,17,18)) + 
  guides(color = "none", shape = guide_legend(title = ""))+
  ggtitle("Nitrogen cycling bacteria in soil")

fig3b
ggsave(paste0("atherton_figure3b_ncyclingsoil_vs_treetype_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig3b, dpi = 300, width = 5.75, height = 3.5, units = "in")

rm(list = c("data", "mean", "melt.mean", "melt.se", "metadata_soil_16s",
            "model_chitin", "model_denit", "model_dnr", "model_nit", "se", 
            "vdata"))

### FIGURE 3C #################################################################
# Select the soil, top layer from ITS
metadata_soil_its <- metadata_its[which(metadata_its$sample_type == "Soil, 0-15 cm"),]

# For figure caption: number of trees' data used (n = 87)
print("Number of trees used in the soil analysis:")
length(unique(metadata_soil_its$tree_id))

# Create the x-axis tick marks
metadata_soil_its$site_edge <- paste0(metadata_soil_its$tree_urban_type, " ", 
                                      metadata_soil_its$tree_edge_type)
metadata_soil_its$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                                    "Urban Street Tree",
                                    metadata_soil_its$site_edge)
metadata_soil_its$site_edge <- factor(metadata_soil_its$site_edge, 
                                      levels = c("Urban Street Tree", 
                                                 "Urban Forest Edge", 
                                                 "Urban Forest Interior", 
                                                 "Rural Forest Edge", 
                                                 "Rural Forest Interior"))

# Select the abundances of the saprotrophs
data <- cbind(metadata_soil_its$perc_dung_saprotroph, 
              metadata_soil_its$perc_wood_saprotroph, 
              metadata_soil_its$perc_litter_saprotroph, 
              metadata_soil_its$perc_soil_saprotroph)

# Rename the columns
colnames(data) <- c("Dung","Wood","Litter","Soil")

# Calculate the mean abundance
mean <- aggregate(data~site_edge, data = metadata_soil_its, FUN= "mean", 
                  na.action = na.omit)

# Reformat the dataframe
vdata <- melt(mean)

# Define the color pallette
barchart_pal <- c("#B3CDDE", "#669BBC", "#003049", "#001825")

# Order the decomposers in decreasing abundance
vdata <- vdata[order(vdata$variable, decreasing = TRUE),]

### STATISTICS TOTAL DECOMPOSER ABUNDANCE ###
metadata_soil_its$sap_tot <- metadata_soil_its$perc_dung_saprotroph + 
  metadata_soil_its$perc_wood_saprotroph + 
  metadata_soil_its$perc_litter_saprotroph + 
  metadata_soil_its$perc_soil_saprotroph

model_total <- lme(formula(sap_tot ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_its)
summary(model_total)
summary(emmeans(model_total, pairwise ~ site_edge))

### STATISTICS DUNG DECOMPOSERS ###
model_dung <- lme(formula(perc_dung_saprotroph ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_its)
summary(model_dung)
summary(emmeans(model_dung, pairwise ~ site_edge))

### STATISTICS WOOD DECOMPOSERS ###
model_wood <- lme(formula(perc_wood_saprotroph ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_its)
summary(model_wood)
summary(emmeans(model_wood, pairwise ~ site_edge))

### STATISTICS LITTER DECOMPOSERS ###
model_litter <- lme(formula(perc_litter_saprotroph ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_its)
summary(model_litter)
summary(emmeans(model_litter, pairwise ~ site_edge))

### STATISTICS SOIL DECOMPOSERS ###
model_soil <- lme(formula(perc_soil_saprotroph ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_soil_its)
summary(model_soil)
summary(emmeans(model_soil, pairwise ~ site_edge))

# Format the data and find the point where the Tukey groups should be placed
vdata <- vdata %>%
  group_by(site_edge) %>%
  mutate(cumulative = cumsum(value),            # Cumulative sum
         midpoint = cumulative - value / 2)     # Midpoint of each segment
vdata$midpoint[17:20] <- vdata$cumulative[17:20] + 0.4
vdata_cumsum <- vdata %>%
  group_by(site_edge) %>%
  summarise(cumulative_sum = sum(value))

# Make the Tukey groups text labels
vdata$tukey <- c("b", "a", "a", "a", "a", rep(c("a", "b", "b", "b", "b"), 3))
vdata$color <- c(rep("white",10), rep("black",10))

# Format the x-axis labels
vdata$site_edge <- gsub(" ", "\n", vdata$site_edge)
vdata$site_edge <- factor(vdata$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

# Plot Figure 3c
fig3c <- ggplot(vdata,  mapping=aes(x=site_edge, y=value, fill=variable))+
  # this makes the stacked barchart
  geom_bar(aes(), stat="identity", position="stack",color="black")+
  # this makes the Tukey group labels
  geom_text(aes(y = midpoint, label = tukey), color = vdata$color, family = "Arial") +
  # this makes the between-group whole bar asterisks
  theme_classic()+ geom_signif(
    comparisons = list(c("Urban\nStreet\nTree", "Urban\nForest\nEdge"),
                       c("Urban\nStreet\nTree", "Urban\nForest\nInterior"),
                       c("Urban\nStreet\nTree", "Rural\nForest\nEdge"),
                       c("Urban\nStreet\nTree", "Rural\nForest\nInterior")),  # Specify pairwise comparisons
    y_position = c(max(vdata_cumsum$cumulative_sum) + 0.1,
                   max(vdata_cumsum$cumulative_sum) + 0.7,
                   max(vdata_cumsum$cumulative_sum) + 1.3,
                   max(vdata_cumsum$cumulative_sum) + 1.9),  # Position above the bars
    annotations = c("**","***", "**", "***"), textsize = 2.5) +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_markdown(size = 8), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), 
        legend.position = c(0.7,1), 
        legend.direction="horizontal") + 
  scale_fill_manual(values = barchart_pal) + 
  labs (x="",y="log(Fungal decomposer relative abundance (%))", fill="") + 
  ggtitle("Fungal decomposers in soil")

fig3c
ggsave(paste0("atherton_figure3c_decomposerssoil_vs_treetype_", 
              format(Sys.Date(), "%Y%m%d"),".pdf"), device = cairo_pdf, 
       plot = fig3c, dpi = 300, width = 5.75, height = 3.5, units = "in")

rm(list = c("data", "mean", "metadata_soil_its", "vdata", "vdata_cumsum"))

### FIGURE 3D #################################################################
# Select leaf data from 16S
metadata_leaf_16s <- metadata_16s[which(metadata_16s$sample_type == "Leaf"),]

# For figure caption: number of trees' data used (n = 62)
print("Number of trees used in the soil analysis:")
length(unique(metadata_leaf_16s$tree_id))

# Format the x-axis tick marks
metadata_leaf_16s$site_edge <- paste0(metadata_leaf_16s$tree_urban_type, " ", 
                                      metadata_leaf_16s$tree_edge_type)
metadata_leaf_16s$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                                    "Urban Street Tree",
                                    metadata_leaf_16s$site_edge)
metadata_leaf_16s$site_edge <- factor(metadata_leaf_16s$site_edge, 
                                      levels = c("Urban Street Tree", 
                                                 "Urban Forest Edge", 
                                                 "Urban Forest Interior",
                                                 "Rural Forest Edge", 
                                                 "Rural Forest Interior"))

# Take the mean and standard error of the data
mean_leaf <- aggregate(perc_methanotroph~site_edge, data = metadata_leaf_16s,
                       FUN= "mean", na.action = na.omit)
se_leaf <- aggregate(perc_methanotroph~site_edge, data = metadata_leaf_16s,
                     FUN= standard_error, na.action = na.omit)

# Make plotting dataframe
vdata <- cbind(mean_leaf, se_leaf[,ncol(se_leaf)],"Leaf")

# Rename columns
colnames(vdata) <- c("site_edge","Abundance","SE","Microbe")

### STATISTICS METHANOTROPHS
model_meth <- lme(formula(perc_methanotroph ~ site_edge), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_leaf_16s)
summary(model_meth)
summary(emmeans(model_meth, pairwise ~ site_edge))

# Create location for Tukey group labels on plot
vdata$max <- c(0.1, 2.6, 2.5, 2, 2.9)

# Define Tukey group labels
vdata$tukey <- c("b", "ab", "a", "ab", "a")

# Format the x-axis tick marks for the plot
vdata$site_edge <- gsub(" ", "\n", vdata$site_edge)
vdata$site_edge <- factor(vdata$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

metadata_leaf_16s$site_edge <- gsub(" ", "\n", metadata_leaf_16s$site_edge)
metadata_leaf_16s$site_edge <- factor(metadata_leaf_16s$site_edge, 
                                      levels = c("Urban\nStreet\nTree", 
                                                 "Urban\nForest\nEdge",
                                                 "Urban\nForest\nInterior",
                                                 "Rural\nForest\nEdge",
                                                 "Rural\nForest\nInterior"))

# Plot figure 3d
fig3d <- ggplot() +
  # This makes the line that follows the averages across the points
  geom_line(data = vdata, 
            aes(x=site_edge, y=Abundance, group=Microbe), 
            linetype = 2, color = "black") +
  # This makes the average abundance point
  geom_point(data = vdata,
             aes(x=site_edge, y=Abundance, color=site_edge, shape=Microbe), 
             size=5) +
  # This makes the individual data points
  geom_point(data = metadata_leaf_16s,
             aes(x = site_edge, y = perc_methanotroph, color = site_edge),
             position = position_jitter(width = 0.2), size = 2, alpha = 0.2, 
             stroke = NA) +
  # This makes the error bars
  geom_errorbar(data = vdata,
                aes(x=site_edge, ymin=Abundance-SE, ymax=Abundance+SE,
                    color=site_edge), width = 0.2) +
  # This makes the Tukey group labels
  geom_text(data = vdata, aes(x = site_edge, y = max, label = tukey), size = 3,
            family = "Arial") +
  theme_classic() + labs(y="log(Methanotrophic bacteria relative abundance (%))", 
                         x="") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_markdown(size = 8), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10),
        legend.position = "none") +
  scale_color_manual(values = categorical_pal) + 
  scale_fill_manual(values = categorical_pal) +
  guides(color = "none", shape = guide_legend(title = ""))+
  ggtitle("Methanotrophic bacterial abundance in leaves")

fig3d
ggsave(paste0("atherton_figure3d_methanotrophicleaves_vs_treetype_", 
              format(Sys.Date(), "%Y%m%d"),".pdf"), device = cairo_pdf, 
       plot = fig3d, dpi = 300, width = 5.75, height = 3.5, units = "in")
