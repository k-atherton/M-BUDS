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
metadata_its <- read.csv("atherton_sample_metadata_ITS_20250125.csv")
metadata_its <- distinct(metadata_its)

# read in sample names of cleaned data
setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/02_Clean_Data/05_Transform_Data/ITS/Aitchison_Distance_Tables/")
all_aitch_its <- read.csv("atherton_ITS_allsampletypes_aitchisondistance_20240919.csv", 
                          row.names = 1)

# read in sequence data
setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/02_Clean_Data/05_Transform_Data/ITS/Rarefied_ASV_Tables/")
o_rare_its <- read.csv("atherton_ITS_osoil_ASV_table_rare_taxonomy_guild_ASVID_20250130.csv",
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

# Normalize the data to be used in this figure
metadata_its$perc_wood_saprotroph <- log(metadata_its$perc_wood_saprotroph+1)

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

### FIGURE S1C #################################################################
# Subset the wood saprotrophs from the reads in the rarefied O Soil ITS data
o_sap_paths <- subset(o_rare_its, primary_lifestyle %in% c("wood_saprotroph"))

# Subset the Plant pathogens from the wood saprotrophs
o_sap_wood_paths <- subset(o_sap_paths, Plant_pathogenic_capacity_template %in% 
                             c("leaf/fruit/seed_pathogen", "root_pathogen", 
                               "wood_pathogen"))

# Subset the human pathogens from the wood saprotrophs
o_sap_human_paths <- subset(o_sap_paths, Animal_biotrophic_capacity_template 
                            %in% c("opportunistic_human_parasite"))

# Find the total sum of the plant pathogen / wood saprotrophs
plant_path_abund <- colSums(o_sap_wood_paths[12:253])
plant_path_abund <- data.frame(plant_path_abund)

# Find the total sum of the human pathogen / wood saprotrophs
human_path_abund <- colSums(o_sap_human_paths[12:253])
human_path_abund <- data.frame(human_path_abund)

# Make the human and plant pathogen data one dataframe
microbe.data<-merge(human_path_abund, plant_path_abund, by = 0)

# Format the metadata
metadata_its$sample_name <- gsub(" ","", metadata_its$sample_name)

# Pull the relevant metadata for the samples for plotting
sample_type<-metadata_its$sample_type[which(metadata_its$sample_name %in% 
                                              microbe.data$Row.names)]
soil_horizon<-metadata_its$soil_horizon[which(metadata_its$sample_name %in% 
                                                microbe.data$Row.names)]
DFB<-metadata_its$tree_dist_to_boston_km[which(metadata_its$sample_name %in% 
                                                 microbe.data$Row.names)]
DFE<-metadata_its$tree_distance_from_edge[which(metadata_its$sample_name %in% 
                                                  microbe.data$Row.names)]
plot_name<-metadata_its$plot_name[which(metadata_its$sample_name %in% 
                                          microbe.data$Row.names)]
site_name<-metadata_its$site_name[which(metadata_its$sample_name %in% 
                                          microbe.data$Row.names)]
tree_id<-metadata_its$tree_id[which(metadata_its$sample_name %in% 
                                      microbe.data$Row.names)]
tree_edge_type<-metadata_its$tree_edge_type[which(metadata_its$sample_name %in%
                                                    microbe.data$Row.names)]
tree_urban_type<-metadata_its$tree_urban_type[which(metadata_its$sample_name %in% 
                                                      microbe.data$Row.names)]
new_tree_type<-paste(tree_urban_type, tree_edge_type, sep = " ")
new_tree_type<-gsub("Urban Street Tree Urban Street Tree", "Urban Street Tree", 
                    new_tree_type)
wood_paths<-microbe.data$plant_path_abund[which(microbe.data$Row.names %in% 
                                                  metadata_its$sample_name)]
human_paths<-microbe.data$human_path_abund[which(microbe.data$Row.names %in% 
                                                   metadata_its$sample_name)]

wood_paths_df <- as.data.frame(cbind(wood_paths, new_tree_type))
human_paths_df <- as.data.frame(cbind(human_paths, new_tree_type))
wood_paths_df$wood_paths <- as.numeric(wood_paths)
human_paths_df$human_paths <- as.numeric(human_paths)
colnames(wood_paths_df) <- c("rel_abund", "tree_type")
colnames(human_paths_df) <- c("rel_abund", "tree_type")
wood_paths_df$pathogen_type <- "Potential plant pathogen"
human_paths_df$pathogen_type <- "Potential human pathogen"
vdata <- rbind(wood_paths_df, human_paths_df)

### STATISTICS HUMAN PATHOGENS ###
model_human <- lme(formula(log(human_paths+1) ~ new_tree_type),
              random = ~1|plot_name/tree_id, na.action = na.exclude,
              method = "ML")
summary(model_human)
summary(emmeans(model_human, pairwise ~ new_tree_type))
anova(model_human)

### STATISTICS WOOD PATHOGENS ###
model_plant <- lme(formula(log(wood_paths+1) ~ new_tree_type), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML")
summary(model_plant)
summary(emmeans(model_plant, pairwise ~ new_tree_type))
anova(model_plant)

# Format the x-axis tickmarks
vdata$site_edge <- gsub(" ", "\n", vdata$tree_type)
vdata$site_edge <- factor(vdata$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

# Order the functional guilds as a factor for plotting
vdata$pathogen_type <- factor(vdata$pathogen_type, 
                              levels = c("Potential plant pathogen",
                                         "Potential human pathogen"))

  # Summarize the data for adding the error bars, mean, and standard error
  summary_data <- vdata %>%
  group_by(site_edge, pathogen_type) %>%
  summarise(
    mean = mean(rel_abund, na.rm = TRUE),
    se = sd(rel_abund, na.rm = TRUE) / sqrt(n()),
    sd = sd(rel_abund, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot Figure 1b
figS1c <- ggplot(vdata, aes(x = site_edge, color = pathogen_type)) +
  # this adds the jittered points
  # geom_jitter(aes(y = rel_abund, group = pathogen_type, stroke = NA), 
  #             size = 0.75,  alpha = 0.2, 
  #             position = position_jitterdodge(jitter.width = 0.5, 
  #                                             dodge.width = 0.75)) +
  # this adds the lines that follow the mean across tree types
  stat_summary(data = vdata,
               aes(y = rel_abund, group = pathogen_type, color = pathogen_type),
               fun = mean, geom = "line",
               position = position_dodge(width = 0.75), linetype = "dashed",
               linewidth = 0.8, alpha = 0.5) +
  # this adds the error bars
  geom_errorbar(data = summary_data,
                aes(x = site_edge, ymin = mean - sd, ymax = mean + sd,
                    group = pathogen_type, color = pathogen_type),
                position = position_dodge(width = 0.75), width = 0.1,
                alpha = 0.5) +
  # this adds the box with the mean and standard error
  geom_crossbar(data = summary_data,
                aes(x = site_edge, y = mean, ymin = mean - se, ymax = mean + se,
                    fill = pathogen_type),
                position = position_dodge(width = 0.75), width = 0.5,
                alpha = 0.5, fatten = 1) +
  theme_classic() +
  labs(y = "log(Potential pathogens\nrelative abundance (%))", x = "", 
       color = "", fill = "") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color = "black"), 
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 0),
        legend.key.size = unit(0.15, "in"),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(size = 10),
        legend.position = c(0.8, 0.9)) +
  scale_color_manual(values = categorical_pal[c(1, 3, 5)]) +
  scale_fill_manual(values = categorical_pal[c(1, 3, 5)]) +
  ggtitle("Wood decomposer potential human and plant pathogen relative abundance in soils")


figS1c
ggsave(paste0("atherton_figureS1c_pathogenicwooddecomposersoils_vs_treetype_boxplot_",
              format(Sys.Date(), "%Y%m%d"),".pdf"), 
       plot = figS1c, dpi = 300, width = 5.75, height = 3.5, units = "in", 
       device = cairo_pdf)
