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

# read in sample names of cleaned data
setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/02_Clean_Data/05_Transform_Data/16S/Aitchison_Distance_Tables/")
bg_aitch_16s <- read.csv("atherton_16S_allsampletypes_aitchisondistance_20240925.csv", 
                         row.names = 1)
leaf_aitch_16s <- read.csv("atherton_16S_leaf_aitchisondistance_20241112.csv", 
                           row.names = 1)

### FORMAT METADATA ###########################################################
# Combine 16S soil horizon with sample type
metadata_16s$sample_type<- paste0(metadata_16s$soil_horizon, " ",
                                  metadata_16s$sample_type)
metadata_16s$sample_type <- gsub("NA ", "", 
                                 metadata_16s$sample_type)

# Replace 16S soil horizon category with depth measurement
metadata_16s$sample_type <- gsub("O Soil", "Soil, 0-15 cm",
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

# Normalize the data to be used in this figure
metadata_16s$perc_chitinolytic <- log(metadata_16s$perc_chitinolytic+1)
metadata_16s$perc_partial_nitrification <- log(metadata_16s$perc_partial_nitrification+1)
metadata_16s$perc_dissim_nitrate_reduction <- log(metadata_16s$perc_dissim_nitrate_reduction+1)
metadata_16s$perc_denitrification <- log(metadata_16s$perc_denitrification+1)

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

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

# Take the functional guild abundances part of the metadata from soil
metadata_chitin <- metadata_soil_16s[which(colnames(metadata_soil_16s) %in% 
                                             c("sample_name", "tree_id", 
                                               "plot_name", "tree_edge_type", 
                                               "tree_urban_type", 
                                               "perc_chitinolytic", 
                                               "site_edge"))]
metadata_nit <- metadata_soil_16s[which(colnames(metadata_soil_16s) %in%
                                          c("sample_name", "tree_id", 
                                            "plot_name", "tree_edge_type", 
                                            "tree_urban_type", 
                                            "perc_partial_nitrification", 
                                            "site_edge"))]
metadata_dnr<- metadata_soil_16s[which(colnames(metadata_soil_16s) %in%
                                         c("sample_name", "tree_id", 
                                           "plot_name", "tree_edge_type", 
                                           "tree_urban_type", 
                                           "perc_dissim_nitrate_reduction", 
                                           "site_edge"))]
metadata_denit <- metadata_soil_16s[which(colnames(metadata_soil_16s) %in%
                                            c("sample_name", "tree_id", 
                                              "plot_name", "tree_edge_type", 
                                              "tree_urban_type", 
                                              "perc_denitrification", 
                                              "site_edge"))]

# Rename the columns
colnames(metadata_chitin) <- c("sample_name", "tree_id", "plot_name", 
                               "tree_edge_type", "tree_urban_type", 
                               "rel_abund", "site_edge")
colnames(metadata_nit) <- c("sample_name", "tree_id", "plot_name",
                            "tree_edge_type", "tree_urban_type", 
                            "rel_abund", "site_edge")
colnames(metadata_dnr) <- c("sample_name", "tree_id", "plot_name", 
                            "tree_edge_type", "tree_urban_type", 
                            "rel_abund", "site_edge")
colnames(metadata_denit) <- c("sample_name", "tree_id", "plot_name", 
                              "tree_edge_type", "tree_urban_type", 
                              "rel_abund", "site_edge")

# Add the name of the functional guild
metadata_chitin$guild <- "Chitinolytic bacteria"
metadata_nit$guild <- "Nitrifying bacteria"
metadata_dnr$guild <- "Dissimilatory nitrate reducing bacteria"
metadata_denit$guild <- "Denitrifying bacteria"

# Combine the soil and leaf datasets for plotting
vdata <- rbind(metadata_chitin, metadata_nit, metadata_dnr, metadata_denit)

# Order the microbial groups for the plot
vdata$guild <- factor(vdata$guild, levels = c("Chitinolytic bacteria", 
                                              "Nitrifying bacteria",
                                              "Dissimilatory nitrate reducing bacteria",
                                              "Denitrifying bacteria"))

### STATISTICS CHITINOLYTIC BACTERIA ###
model_chitin <- lme(formula(perc_chitinolytic ~ site_edge), 
                    random = ~1|plot_name/tree_id, na.action = na.exclude, 
                    method = "ML", data = metadata_soil_16s)
summary(model_chitin)
summary(emmeans(model_chitin, pairwise ~ site_edge))

### STATISTICS DISSIMILATORY NITRATE REDUCING BACTERIA ###
model_dnr <- lme(formula(perc_dissim_nitrate_reduction ~ site_edge), 
                 random = ~1|plot_name/tree_id, na.action = na.exclude, 
                 method = "ML", data = metadata_soil_16s)
summary(model_dnr)
summary(emmeans(model_dnr, pairwise ~ site_edge))

### STATISTICS NITRIFYING BACTERIA ###
model_nit <- lme(formula(perc_partial_nitrification ~ site_edge), 
                 random = ~1|plot_name/tree_id, na.action = na.exclude, 
                 method = "ML", data = metadata_soil_16s)
summary(model_nit)
summary(emmeans(model_nit, pairwise ~ site_edge))

### STATISTICS DENITRIFYING BACTERIA ###
model_denit <- lme(formula(perc_denitrification ~ site_edge), 
                   random = ~1|plot_name/tree_id, na.action = na.exclude, 
                   method = "ML", data = metadata_soil_16s)
summary(model_denit)
summary(emmeans(model_denit, pairwise ~ site_edge))

# Make the text for adding Tukey groups and record the location of the group id on the y-axis
vdata$tukey <- NA
vdata$location <- NA

for(i in 1:nrow(vdata)){
  guild <- as.character(vdata$guild[i])
  site_edge <- as.character(vdata$site_edge[i])
  
  if(guild == "Chitinolytic bacteria"){
    if(site_edge == "Urban Street Tree"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    }
  } else if(guild == "Nitrifying bacteria"){
    if(site_edge == "Urban Street Tree"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    }
  } else if(guild == "Dissimilatory nitrate reducing bacteria"){
    if(site_edge == "Urban Street Tree"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    }
  } else if(guild == "Denitrifying bacteria"){
    if(site_edge == "Urban Street Tree"){
      vdata$tukey[i] <- "a"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Urban Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Edge"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    } else if(site_edge == "Rural Forest Interior"){
      vdata$tukey[i] <- "b"
      vdata$location[i] <- max(vdata$rel_abund[which(vdata$site_edge == site_edge & vdata$guild == guild)]) + 0.2
    }
  }
}

# Reformat the x-axis tick marks for the plot
vdata$site_edge <- gsub(" ", "\n", vdata$site_edge)
vdata$site_edge <- factor(vdata$site_edge, 
                          levels = c("Urban\nStreet\nTree", 
                                     "Urban\nForest\nEdge",
                                     "Urban\nForest\nInterior",
                                     "Rural\nForest\nEdge",
                                     "Rural\nForest\nInterior"))

# Summarize the data for adding the error bars, mean, and standard error
summary_data <- vdata %>%
  group_by(site_edge, guild) %>%
  summarise(
    mean = mean(rel_abund, na.rm = TRUE),
    se = sd(rel_abund, na.rm = TRUE) / sqrt(n()),
    sd = sd(rel_abund, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot Figure 3b
fig3b <- ggplot(vdata, aes(x = site_edge, color = guild)) +
  # this adds the jittered points
  geom_jitter(aes(y = rel_abund, group = guild, stroke = NA), 
              size = 0.75,  alpha = 0.2, 
              position = position_jitterdodge(jitter.width = 0.5, 
                                              dodge.width = 0.75)) +
  # this adds the lines that follow the mean across tree types
  stat_summary(data = vdata,
               aes(y = rel_abund, group = guild, color = guild),
               fun = mean, geom = "line",
               position = position_dodge(width = 0.75), linetype = "dashed",
               linewidth = 0.8, alpha = 0.5) +
  # this adds the error bars
  geom_errorbar(data = summary_data,
                aes(x = site_edge, ymin = mean - sd, ymax = mean + sd,
                    group = guild, color = guild),
                position = position_dodge(width = 0.75), width = 0.1,
                alpha = 0.5) +
  # this adds the box with the mean and standard error
  geom_crossbar(data = summary_data,
                aes(x = site_edge, y = mean, ymin = mean - se, ymax = mean + se,
                    fill = guild),
                position = position_dodge(width = 0.75), width = 0.5,
                alpha = 0.5, fatten = 1) +
  # this adds the tukey groups
  geom_text(data = vdata, family = "Arial",
            aes(x = site_edge, y = location, color = guild, group = guild,
                label = tukey),
            position = position_dodge(width = 0.75), size = 3, parse = TRUE,
            show.legend = FALSE) +
  theme_classic() +
  labs(y = "log(Functional group\nrelative abundance (%))", x = "", color = "", fill = "") +
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
        legend.position = c(0.7, 1.0)) +
  scale_color_manual(values = categorical_pal[c(1,2,3,4)]) +
  scale_fill_manual(values = categorical_pal[c(1,2,3,4)]) +
  ggtitle("Nitrogen cycling bacteria in soil")

fig3b
ggsave(paste0("atherton_figure3b_ncyclingsoil_vs_treetype_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig3b, dpi = 300, width = 5.75, height = 3.5, units = "in")
