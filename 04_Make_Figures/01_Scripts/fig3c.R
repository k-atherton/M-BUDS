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
metadata_its$perc_dung_saprotroph <- log(metadata_its$perc_dung_saprotroph+1)
metadata_its$perc_wood_saprotroph <- log(metadata_its$perc_wood_saprotroph+1)
metadata_its$perc_litter_saprotroph <- log(metadata_its$perc_litter_saprotroph+1)
metadata_its$perc_soil_saprotroph <- log(metadata_its$perc_soil_saprotroph+1)

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

### FIGURE 3C #################################################################
# Select the soil, top layer from ITS
metadata_soil_its <- metadata_its[which(metadata_its$sample_type == 
                                          "Soil, 0-15 cm"),]

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
        legend.text = element_markdown(size = 10), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), 
        legend.position = c(0.7,0.7), 
        legend.direction="horizontal") + 
  scale_fill_manual(values = barchart_pal) + 
  labs (x="",y="log(Fungal decomposer\nrelative abundance (%))", fill="") + 
  ggtitle("Fungal decomposers in soil")

fig3c
ggsave(paste0("atherton_figure3c_decomposerssoil_vs_treetype_", 
              format(Sys.Date(), "%Y%m%d"),".pdf"), device = cairo_pdf, 
       plot = fig3c, dpi = 300, width = 5.75, height = 3.5, units = "in")
