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

# Format the 16S environmental data
metadata_16s$soil_percent_organic_matter <- gsub("%", "",
                                                 metadata_16s$soil_percent_organic_matter)
metadata_16s$soil_percent_organic_matter <- as.numeric(metadata_16s$soil_percent_organic_matter)

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

### FIGURE S2A #################################################################
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
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == 
                                          "Soil, 0-15 cm"),]

metadata_soil_16s <- metadata_soil_16s[,which(colnames(metadata_soil_16s) %in%
                                                c("tree_id", "plot_name", 
                                                  "soil_percent_organic_matter", 
                                                  "site_edge"))]

# For figure caption: number of trees' data used (n = 88)
print("Number of trees used in the soil analysis:")
length(unique(metadata_soil_16s$tree_id))

### STATISTICS SOIL NO3- ###
model_no3 <- lme(formula(soil_percent_organic_matter ~ site_edge), 
                 random = ~1|plot_name/tree_id, na.action = na.exclude, 
                 method = "ML", data = metadata_soil_16s)
summary(model_no3)
summary(emmeans(model_no3, pairwise ~ site_edge))
anova(model_no3)

# Create the Tukey groups
for(i in 1:nrow(metadata_soil_16s)){
  site_edge <- as.character(metadata_soil_16s$site_edge[i])
  if(site_edge == "Urban Street Tree"){
    metadata_soil_16s$tukey[i] <- "a"
    metadata_soil_16s$location[i] <- max(metadata_soil_16s$soil_percent_organic_matter[which(metadata_soil_16s$site_edge == site_edge)]) + 0.05
  } else if(site_edge == "Urban Forest Edge"){
    metadata_soil_16s$tukey[i] <- "ab"
    metadata_soil_16s$location[i] <- max(metadata_soil_16s$soil_percent_organic_matter[which(metadata_soil_16s$site_edge == site_edge)]) + 0.05
  } else if(site_edge == "Urban Forest Interior"){
    metadata_soil_16s$tukey[i] <- "b"
    metadata_soil_16s$location[i] <- max(metadata_soil_16s$soil_percent_organic_matter[which(metadata_soil_16s$site_edge == site_edge)]) + 0.05
  } else if(site_edge == "Rural Forest Edge"){
    metadata_soil_16s$tukey[i] <- "ab"
    metadata_soil_16s$location[i] <- max(metadata_soil_16s$soil_percent_organic_matter[which(metadata_soil_16s$site_edge == site_edge)]) + 0.05
  } else if(site_edge == "Rural Forest Interior"){
    metadata_soil_16s$tukey[i] <- "b"
    metadata_soil_16s$location[i] <- max(metadata_soil_16s$soil_percent_organic_matter[which(metadata_soil_16s$site_edge == site_edge)]) + 0.05
  }
}

# Reformat the x-axis tickmark labels
metadata_soil_16s$site_edge <- gsub(" ", "\n", metadata_soil_16s$site_edge)
metadata_soil_16s$site_edge <- factor(metadata_soil_16s$site_edge, 
                                      levels = c("Urban\nStreet\nTree", 
                                                 "Urban\nForest\nEdge",
                                                 "Urban\nForest\nInterior",
                                                 "Rural\nForest\nEdge",
                                                 "Rural\nForest\nInterior"))


# Summarize the data for adding the error bars, mean, and standard error
summary_data <- metadata_soil_16s %>%
  group_by(site_edge) %>%
  summarise(
    mean = mean(soil_percent_organic_matter, na.rm = TRUE),
    se = sd(soil_percent_organic_matter, na.rm = TRUE) / sqrt(n()),
    sd = sd(soil_percent_organic_matter, na.rm = TRUE),
    .groups = 'drop'
  )

# Visualize
figS2a <- ggplot(metadata_soil_16s, aes(x = site_edge, y = soil_percent_organic_matter, 
                                       group = 1)) +
  # this adds the jittered points
  geom_jitter(aes(stroke = NA), 
              size = 0.75,  alpha = 0.2, color = "#FF0000") +
  # this adds the lines that follow the mean across tree types
  stat_summary(fun = mean, geom = "line",
               position = position_dodge(width = 0.75), linetype = "dashed",
               linewidth = 0.8, alpha = 0.5, color = "#FF0000") +
  # this adds the error bars
  geom_errorbar(data = summary_data,
                aes(x = site_edge, y = mean, ymin = mean - sd, ymax = mean + sd,
                    group = site_edge),
                position = position_dodge(width = 0.75), width = 0.1,
                alpha = 0.5, color = "#FF0000") +
  # # this adds the box with the mean and standard error
  geom_crossbar(data = summary_data,
                aes(x = site_edge, y = mean, ymin = mean - se, 
                    ymax = mean + se),
                position = position_dodge(width = 0.75), width = 0.5,
                alpha = 0.5, fatten = 1, , color = "#FF0000", 
                fill = "#FF0000") +
  # # this adds the tukey groups
  geom_text(family = "Arial",
            aes(y = location, label = tukey),
            position = position_dodge(width = 0.75), size = 3, parse = TRUE,
            show.legend = FALSE) +
  theme_classic() + 
  labs (y="Soil organic matter content (%)", x="") +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_markdown(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_markdown(size = 10, color = "black"), 
        legend.text = element_markdown(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_markdown(size = 10), 
        legend.position = "none", 
        legend.background = element_rect(fill = NA)) +
  guides(color = "none", shape = guide_legend(title = ""))+
  ggtitle("Soil organic matter content")

figS2a
ggsave(paste0("atherton_figureS2a_organicmattersoil_vs_treetype_boxplot_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = figS2a, dpi = 300, width = 6, height = 3.5, units = "in")
