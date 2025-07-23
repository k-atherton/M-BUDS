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
metadata_its$perc_epiphyte <- log(metadata_its$perc_epiphyte+1)

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

### FIGURE 1A #################################################################
# Select soil, top layer data from ITS
metadata_soil_its <- metadata_its[which(metadata_its$sample_type %in% 
                                          c("Soil, 0-15 cm")),]

# For figure caption: number of trees' data used (n = 87)
print("Number of trees used in the soil analysis:")
length(unique(metadata_soil_its$tree_id))

# Select leaf data from ITS
metadata_leaf_its <- metadata_its[which(metadata_its$sample_type %in% 
                                          c("Leaf")),]

# For figure caption: number of trees' data used (n = 83)
print("Number of trees used in the leaf analysis:")
length(unique(metadata_leaf_its$tree_id))

# Take the fungal guild abundances part of the metadata from soil and leaves
fungal_trait_soil <- metadata_soil_its[, colnames(metadata_soil_its) %in% 
                                         "perc_ecm"]
fungal_trait_leaf <- metadata_leaf_its[, colnames(metadata_leaf_its) %in% 
                                         "perc_epiphyte"]

# Calculate the mean of the guild abundances for each tree type
mean_soil <- aggregate(fungal_trait_soil,
                       by=list(metadata_soil_its$tree_urban_type,
                               metadata_soil_its$tree_edge_type), FUN= "mean")
colnames(mean_soil)[3] <- "ecm_mean"
mean_leaf <- aggregate(fungal_trait_leaf,
                       by=list(metadata_leaf_its$tree_urban_type,
                               metadata_leaf_its$tree_edge_type), FUN= "mean")
colnames(mean_leaf)[3] <- "epiphyte_mean"

# Calculate the standard error of the guild abundances for each tree type
se_soil <- aggregate(fungal_trait_soil,
                     by=list(metadata_soil_its$tree_urban_type,
                             metadata_soil_its$tree_edge_type), 
                     FUN= standard_error)
colnames(se_soil)[3] <- "ecm_se"
se_leaf <- aggregate(fungal_trait_leaf,
                     by=list(metadata_leaf_its$tree_urban_type,
                             metadata_leaf_its$tree_edge_type), 
                     FUN= standard_error)
colnames(se_leaf)[3] <- "epiphyte_se"

# Combine the data for just the guilds analyzed in this figure (ECM and epiphytes)
vdata <- merge(mean_soil, mean_leaf, by = intersect(names(mean_soil), 
                                                    names(mean_leaf)))
vdata <- merge(vdata, se_soil, by = intersect(names(vdata), 
                                                    names(se_soil)))
vdata <- merge(vdata, se_leaf, by = intersect(names(vdata), 
                                              names(se_leaf)))

# Name the columns 
colnames(vdata)[1:2] <- c("Site", "Edge")

# Make a long version of the data for visualization
vdata_long <- as.data.frame(matrix(data = NA, nrow = 10, ncol = 5))

# Name the columns
colnames(vdata_long) <- c("Site", "Edge", "Guild", "Abundance", "SE")

# Since we have 2 guilds in the figure, repeat the site and edge information twice
vdata_long$Site <- rep(vdata$Site, 2)
vdata_long$Edge <- rep(vdata$Edge, 2)

# Since there are 5 tree types, repeat the guilds five times
vdata_long$Guild <- c(rep("ECM in soil", 5), 
                      rep("Epiphytes in leaves", 5))

# Input the abundance data for each guild, multiply the epiphyte abundance by 20
# for visualization -- this will be adjusted for in the secondary y-axis
vdata_long$Abundance <- c(vdata$ecm_mean, 20*vdata$epiphyte_mean)

# Input the standard error data for each guild
vdata_long$SE <- c(vdata$ecm_se, 20*vdata$epiphyte_se)

# Create the x-axis categories
vdata_long$site_edge <- paste0(vdata_long$Site, " ", vdata_long$Edge)
vdata_long$site_edge <- gsub("Urban Street Tree Urban Street Tree", 
                             "Urban Street Tree",
                             vdata_long$site_edge)

# Make the tree type a factor with the order for the x-axis
vdata_long$site_edge <- factor(vdata_long$site_edge, 
                               levels = c("Urban Street Tree", 
                                          "Urban Forest Edge",
                                          "Urban Forest Interior",
                                          "Rural Forest Edge", 
                                          "Rural Forest Interior"))

# Make guild a factor for plotting the legend in order
vdata_long$Guild <- factor(vdata_long$Guild, levels = c("ECM in soil", 
                                                        "Epiphytes in leaves"))

### STATISTICS ECM ###
model_ecm <- lme(formula(perc_ecm ~ tree_pit_type), 
                 random = ~1|plot_name/tree_id, na.action = na.exclude, 
                 method = "ML", data = metadata_soil_its)
summary(model_ecm)
summary(emmeans(model_ecm, pairwise ~ tree_pit_type))

### STATISTICS EPIPHYTES ###
model_epiphyte <- lme(formula(perc_epiphyte ~ tree_pit_type), 
                      random = ~1|plot_name/tree_id, na.action = na.exclude, 
                      method = "ML", data = metadata_leaf_its)
summary(model_epiphyte)
summary(emmeans(model_epiphyte, pairwise ~ tree_pit_type))


# Make the text for adding Tukey groups
vdata_long$tukey <- c(
  "ab[\"\\u25CF\"]", # ecm rural edge
  "a[\"\\u25CF\"]", # ecm rural interior
  "a[\"\\u25CF\"]", # ecm urban edge
  "a[\"\\u25CF\"]", # ecm urban interior
  "b[\"\\u25CF\"]", # ecm street tree
  "a[\"\\u25B2\"]", # epiphyte rural edge
  "a[\"\\u25B2\"]", # epiphyte rural interior
  "a[\"\\u25B2\"]", # epiphyte urban edge
  "a[\"\\u25B2\"]", # epiphyte urban interior
  "b[\"\\u25B2\"]" # epiphyte street tree
)

# Format labels for plotting
vdata_long$site_edge <- gsub(" ", "\n", vdata_long$site_edge)
vdata_long$site_edge <- factor(vdata_long$site_edge, 
                               levels = c("Urban\nStreet\nTree", 
                                          "Urban\nForest\nEdge",
                                          "Urban\nForest\nInterior",
                                          "Rural\nForest\nEdge",
                                          "Rural\nForest\nInterior"))

# Plot Figure 1a
fig1a <- ggplot(vdata_long) +
  # this creates the dotted line between points to track the average across tree types
  geom_line(aes(x=site_edge, y=Abundance, group = Guild),
            linetype = 2) +
  # this creates the points that show the average abundance for each guild and tree type
  geom_point(aes(x=site_edge, y=Abundance, color=site_edge, shape = Guild, 
                 group = Guild), size=3) + 
  # this creates the standard error bars around the average abundance points
  geom_errorbar(aes(x=site_edge, ymin=Abundance-SE, ymax=Abundance+SE, 
                    color=site_edge), width = 0.2, size = 0.3) +
  # this adds the Tukey groups labels defined above
  geom_text(family = "Arial",
            aes(x = site_edge, y = Abundance + SE + 3, label = tukey), size = 3, 
            parse = TRUE) +
  theme_classic() + labs (y="ECM relative abundance (%)", x="", shape = "") + 
  # this creates the secondary axis, dividing by 20 to adjust for the multiplication by 20 above
  scale_y_continuous(sec.axis = sec_axis(~./20, 
                                         name = "log(Epiphyte relative abundance (%))")) +
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1),
        legend.background = element_rect(fill='transparent'),
        plot.title = element_text(size = 10), legend.position = c(0.5, 0.2)) +
  scale_color_manual(values = categorical_pal, guide = "none") + 
  ggtitle("Tree-associated fungal symbionts")

fig1a
ggsave(paste0("atherton_figure1a_ecmsoil_epiphytesleaves_vs_treetype_",
              format(Sys.Date(), "%Y%m%d"),".pdf"), 
       plot = fig1a, dpi = 300, width = 6, height = 3.5, units = "in", 
       device = cairo_pdf)
