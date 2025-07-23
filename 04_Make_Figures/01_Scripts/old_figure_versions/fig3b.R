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
metadata_soil_16s <- metadata_16s[which(metadata_16s$sample_type == 
                                          "Soil, 0-15 cm"),]

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
vdata$Microbe <- factor(vdata$Microbe, 
                        levels = c("Chitinolytic bacteria", 
                                   "Nitrifying bacteria",
                                   "Dissimilatory nitrate reducing bacteria",
                                   "Denitrifying bacteria"))

### STATISTICS: CHITINOLYTIC BACTERIA ###
model_chitin <- lme(formula(perc_chitinolytic ~ site_edge), 
                    random = ~1|plot_name/tree_id, na.action = na.exclude, 
                    method = "ML", data = metadata_soil_16s)
summary(model_chitin)
summary(emmeans(model_chitin, pairwise ~ site_edge))

### STATISTICS: DISSIMILATORY NITRATE REDUCING BACTERIA ###
model_dnr <- lme(formula(perc_dissim_nitrate_reduction ~ site_edge), 
                 random = ~1|plot_name/tree_id, na.action = na.exclude, 
                 method = "ML", data = metadata_soil_16s)
summary(model_dnr)
summary(emmeans(model_dnr, pairwise ~ site_edge))

### STATISTICS: NITRIFYING BACTERIA ###
model_nit <- lme(formula(perc_partial_nitrification ~ site_edge), 
                 random = ~1|plot_name/tree_id, na.action = na.exclude, 
                 method = "ML", data = metadata_soil_16s)
summary(model_nit)
summary(emmeans(model_nit, pairwise ~ site_edge))

### STATISTICS: DENITRIFYING BACTERIA ###
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

# Plot figure 3b
fig3b <- ggplot(vdata) +
  # this creates the dotted line between points to track the average across tree types
  geom_line(data = vdata, 
            aes(x=site_edge, y=Abundance, group=Microbe), 
            linetype = 2, color = "black") +
  # this creates the points that show the average abundance for each guild and tree type
  geom_point(data = vdata, 
             aes(x=site_edge, y=Abundance, color=site_edge, shape=Microbe), 
             size=3) + 
  # this creates the standard error bars around the average abundance points
  geom_errorbar(data = vdata, 
                aes(x=site_edge, ymin=Abundance-SE, ymax=Abundance+SE, 
                    color=site_edge), width = 0.2, size = 0.3) +
  # this adds the Tukey groups labels defined above
  geom_text_repel(data = vdata, direction = "y",
                  aes(x=site_edge, y=Abundance, label = tukey), 
                  hjust = 1, segment.color = NA,
                  size = 3, parse = TRUE, nudge_x = -0.1) +
  theme_classic() + labs (y="log(Functional group\nrelative abundance (%))", x="") +
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
