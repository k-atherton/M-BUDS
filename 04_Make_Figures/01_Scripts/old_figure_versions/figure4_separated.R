### LOAD IN LIBRARIES #########################################################
library(dplyr)
library(tidyr)
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

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

### FIGURE 4A #################################################################
# Select the leaf data from ITS
metadata_leaf_its <- metadata_its[which(metadata_its$sample_type == "Leaf"),]

# Select the variables for the figure
metadata_leaf_its_div <- metadata_leaf_its[,colnames(metadata_leaf_its) %in% 
                                             c("sample_name", 
                                               "tree_dist_to_boston_km", 
                                               "shannon", "shannon_nopath")]
# Format the data for plotting
diversity_long <- metadata_leaf_its_div %>% 
  pivot_longer(!c(sample_name, tree_dist_to_boston_km), 
               names_to = "diversity_type", values_to = "value")
diversity_long$diversity_type <- gsub("shannon_nopath", "non-pathogenic fungi",
                                      diversity_long$diversity_type)
diversity_long$diversity_type <- gsub("shannon", "all fungi",
                                      diversity_long$diversity_type)

### STATISTICS ALL FUNGI ###
model_all <- lme(formula(shannon ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_leaf_its)
anova(model_all)
r.squaredGLMM(model_all)

### STATISTICS NON-PATHOGENIC FUNGI ###
model_nopath <- lme(formula(shannon_nopath ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_leaf_its)
anova(model_nopath)
r.squaredGLMM(model_nopath)

# Plot figure 4a
fig4a <- ggplot(diversity_long, aes(x = log(tree_dist_to_boston_km +1), 
                                    y = value, color = diversity_type, 
                                    fill = diversity_type)) + 
  geom_point(shape = 18, alpha = 0.5, stroke = NA, size = 2.5) +
  geom_smooth(method = "lm", alpha = 0.12, 
              aes(ymax = after_stat(y + se * sqrt(length(y))),
                  ymin = after_stat(y - se * sqrt(length(y))))) + 
  scale_color_manual(values = linear_pal) + 
  scale_fill_manual(values = linear_pal) + 
  theme_classic() + 
  labs(x = "log(Distance from Boston (km))", y = "Shannon diversity", legend = NA) + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), legend.position = "none") + 
  ggtitle("Fungal diversity in leaves")

fig4a
ggsave(paste0("atherton_figure4a_ITSdiversityleaves_vs_dfb_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig4a, dpi = 300, width = 5.75, height = 3.5, units = "in")

rm(list = c("diversity_long" ,"metadata_leaf_its", "metadata_leaf_its_div", 
            "model_all", "model_nopath"))

### FIGURE 4B ##################################################################
# Select the leaf data from 16S
metadata_leaf_16s <- metadata_16s[which(metadata_16s$sample_type == "Leaf"),]

# Select the variables for the figure
metadata_leaf_16s_div <- metadata_leaf_16s[,colnames(metadata_leaf_16s) %in% 
                                             c("sample_name", 
                                               "tree_dist_to_boston_km", 
                                               "shannon.x", "shannon_nopath")]
# Format the data for plotting
diversity_long <- metadata_leaf_16s_div %>% 
  pivot_longer(!c(sample_name, tree_dist_to_boston_km), 
               names_to = "diversity_type", values_to = "value")
diversity_long$diversity_type <- gsub("shannon_nopath", "non-pathogenic fungi",
                                      diversity_long$diversity_type)
diversity_long$diversity_type <- gsub("shannon", "all fungi",
                                      diversity_long$diversity_type)

### STATISTICS ALL BACTERIA ###
model_all <- lme(formula(shannon.x ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_leaf_16s)
anova(model_all)
r.squaredGLMM(model_all)

### STATISTICS NON-PATHOGENIC BACTERIA ###
model_nopath <- lme(formula(shannon_nopath ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_leaf_16s)
anova(model_nopath)
r.squaredGLMM(model_nopath)

# Plot figure 4b
fig4b <- ggplot(diversity_long, aes(x = log(tree_dist_to_boston_km +1), 
                                    y = value, color = diversity_type, 
                                    fill = diversity_type)) + 
  geom_point(shape = 18, alpha = 0.5, stroke = NA, size = 2.5) +
  geom_smooth(method = "lm", alpha = 0.12, 
              aes(ymax = after_stat(y + se * sqrt(length(y))),
                  ymin = after_stat(y - se * sqrt(length(y))))) +  
  scale_color_manual(values = linear_pal) + 
  scale_fill_manual(values = linear_pal) + 
  theme_classic() + 
  labs(x = "log(Distance from Boston (km))", y = "Shannon diversity", legend = NA) + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), 
        legend.position = "none") + 
  ggtitle("Bacterial diversity in leaves")

fig4b
ggsave(paste0("atherton_figure4b_16Sdiversityleaves_vs_dfb_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig4b, dpi = 300, width = 5.75, height = 3.5, units = "in")

rm(list = c("diversity_long" ,"metadata_leaf_16s", "metadata_leaf_16s_div", 
            "model_all", "model_nopath"))

### FIGURE 4C #################################################################
# Select the belowground data from ITS
metadata_osoil_its <- metadata_its[which(metadata_its$sample_type == "Soil, 0-15 cm"),]
metadata_msoil_its <- metadata_its[which(metadata_its$sample_type == "Soil, 15-30 cm"),]
metadata_root_its <- metadata_its[which(metadata_its$sample_type %in% 
                                          c("Roots, 0-15 cm", 
                                            "Roots, 15-30 cm")),]

# Select the variables for the figure
metadata_osoil_its_div <- metadata_osoil_its[,colnames(metadata_osoil_its) %in% 
                                             c("sample_name", "sample_type", 
                                               "tree_dist_to_boston_km", 
                                               "shannon", "shannon_nopath")]
metadata_msoil_its_div <- metadata_msoil_its[,colnames(metadata_msoil_its) %in% 
                                               c("sample_name", "sample_type", 
                                                 "tree_dist_to_boston_km", 
                                                 "shannon", "shannon_nopath")]
metadata_root_its_div <- metadata_root_its[,colnames(metadata_root_its) %in% 
                                               c("sample_name", "sample_type", 
                                                 "tree_dist_to_boston_km", 
                                                 "shannon", "shannon_nopath")]

# Format the data for plotting
diversity_long_osoil <- metadata_osoil_its_div %>% 
  pivot_longer(!c(sample_name, sample_type, tree_dist_to_boston_km), 
               names_to = "diversity_type", values_to = "value")
diversity_long_osoil$diversity_type <- gsub("shannon_nopath", "non-pathogenic fungi",
                                            diversity_long_osoil$diversity_type)
diversity_long_osoil$diversity_type <- gsub("shannon", "all fungi",
                                            diversity_long_osoil$diversity_type)

diversity_long_msoil <- metadata_msoil_its_div %>% 
  pivot_longer(!c(sample_name, sample_type, tree_dist_to_boston_km), 
               names_to = "diversity_type", values_to = "value")
diversity_long_msoil$diversity_type <- gsub("shannon_nopath", "non-pathogenic fungi",
                                            diversity_long_msoil$diversity_type)
diversity_long_msoil$diversity_type <- gsub("shannon", "all fungi",
                                            diversity_long_msoil$diversity_type)

diversity_long_root <- metadata_root_its_div %>% 
  pivot_longer(!c(sample_name, sample_type, tree_dist_to_boston_km), 
               names_to = "diversity_type", values_to = "value")
diversity_long_root$diversity_type <- gsub("shannon_nopath", "non-pathogenic fungi",
                                           diversity_long_root$diversity_type)
diversity_long_root$diversity_type <- gsub("shannon", "all fungi",
                                           diversity_long_root$diversity_type)

### STATISTICS ALL FUNGI ###
model_all_osoil <- lme(formula(shannon ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
              random = ~1|plot_name/tree_id, na.action = na.exclude, 
              method = "ML", data = metadata_osoil_its)
anova(model_all_osoil)
r.squaredGLMM(model_all_osoil)

model_all_msoil <- lme(formula(shannon ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
                       random = ~1|plot_name/tree_id, na.action = na.exclude, 
                       method = "ML", data = metadata_msoil_its)
anova(model_all_msoil)
r.squaredGLMM(model_all_msoil)

model_all_root <- lme(formula(shannon ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
                       random = ~1|plot_name/tree_id, na.action = na.exclude, 
                       method = "ML", data = metadata_root_its)
anova(model_all_root)
r.squaredGLMM(model_all_root)

### STATISTICS NON-PATHOGENIC FUNGI ###
model_nopath_osoil <- lme(formula(shannon_nopath ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
                       random = ~1|plot_name/tree_id, na.action = na.exclude, 
                       method = "ML", data = metadata_osoil_its)
anova(model_nopath_osoil)
r.squaredGLMM(model_nopath_osoil)

model_nopath_msoil <- lme(formula(shannon_nopath ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
                       random = ~1|plot_name/tree_id, na.action = na.exclude, 
                       method = "ML", data = metadata_msoil_its)
anova(model_nopath_msoil)
r.squaredGLMM(model_nopath_msoil)

model_nopath_root <- lme(formula(shannon_nopath ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
                      random = ~1|plot_name/tree_id, na.action = na.exclude, 
                      method = "ML", data = metadata_root_its)
anova(model_nopath_root)
r.squaredGLMM(model_nopath_root)

# Merge all sample types into one
diversity_long_osoil$color <- paste0(diversity_long_osoil$diversity_type, " ", 
                                     diversity_long_osoil$sample_type)
diversity_long_msoil$color <- paste0(diversity_long_msoil$diversity_type, " ", 
                                     diversity_long_msoil$sample_type)
diversity_long_root$color <- paste0(diversity_long_root$diversity_type, " ", 
                                     diversity_long_root$sample_type)

diversity_long <- rbind(diversity_long_osoil, diversity_long_msoil, 
                        diversity_long_root)

new_pal <- c(categorical_pal[c(1,2,3,4)], linear_pal, "#808080", "#000000")

# Plot figure 4c
fig4c_osoil <- ggplot(diversity_long_osoil, aes(x = log(tree_dist_to_boston_km +1), 
                                    y = value, color = diversity_type,
                                    fill = diversity_type)) + 
  geom_point(aes(shape = sample_type), alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", alpha = 0.12, 
              aes(ymax = after_stat(y + se * sqrt(length(y))), 
                  ymin = after_stat(y - se * sqrt(length(y))))) + 
  scale_color_manual(values = linear_pal) + 
  scale_fill_manual(values = linear_pal) + 
  theme_classic() + 
  labs(x = "log(Distance from Boston (km))", y = "Shannon diversity", 
       legend = NA) + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), 
        legend.position = "none") + 
  ggtitle("Fungal diversity belowground") +
  scale_shape_manual(values = shapes)

fig4c_osoil
ggsave(paste0("atherton_figure4c_ITSdiversityosoil_vs_dfb_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig4c_osoil, dpi = 300, width = 5.75, height = 3.5, units = "in")

fig4c_msoil <- ggplot(diversity_long_msoil, aes(x = log(tree_dist_to_boston_km +1), 
                                                y = value, color = diversity_type,
                                                fill = diversity_type)) + 
  geom_point(aes(shape = sample_type), alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", alpha = 0.12, 
              aes(ymax = after_stat(y + se * sqrt(length(y))), 
                  ymin = after_stat(y - se * sqrt(length(y))))) + 
  scale_color_manual(values = linear_pal) + 
  scale_fill_manual(values = linear_pal) + 
  theme_classic() + 
  labs(x = "log(Distance from Boston (km))", y = "Shannon diversity", 
       legend = NA) + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), 
        legend.position = "none") + 
  ggtitle("Fungal diversity belowground") +
  scale_shape_manual(values = shapes)

fig4c_msoil
ggsave(paste0("atherton_figure4c_ITSdiversitymsoil_vs_dfb_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig4c_msoil, dpi = 300, width = 5.75, height = 3.5, units = "in")

fig4c_root <- ggplot(diversity_long_root, aes(x = log(tree_dist_to_boston_km +1), 
                                                y = value, color = diversity_type,
                                                fill = diversity_type)) + 
  geom_point(aes(shape = sample_type), alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", alpha = 0.12, 
              aes(ymax = after_stat(y + se * sqrt(length(y))), 
                  ymin = after_stat(y - se * sqrt(length(y))))) + 
  scale_color_manual(values = linear_pal) + 
  scale_fill_manual(values = linear_pal) + 
  theme_classic() + 
  labs(x = "log(Distance from Boston (km))", y = "Shannon diversity", 
       legend = NA) + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), 
        legend.position = "none") + 
  ggtitle("Fungal diversity belowground") +
  scale_shape_manual(values = shapes)

fig4c_root
ggsave(paste0("atherton_figure4c_ITSdiversityroot_vs_dfb_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig4c_root, dpi = 300, width = 5.75, height = 3.5, units = "in")


fig4c_wlegend <- ggplot(diversity_long, aes(x = log(tree_dist_to_boston_km +1), 
                                            y = value, color = diversity_type, 
                                            fill = diversity_type)) + 
  geom_point(aes(shape = sample_type), alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", alpha = 0.12, 
              aes(ymax = after_stat(y + se * sqrt(length(y))), 
                  ymin = after_stat(y - se * sqrt(length(y))))) + 
  scale_color_manual(values = linear_pal) + 
  scale_fill_manual(values = linear_pal) + 
  theme_classic() + 
  labs(x = "log(Distance from Boston (km))", y = "Shannon diversity", shape = NA) + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_blank(), 
        plot.title = element_text(size = 10), 
        legend.position = "bottom", legend.box = "vertical", 
        legend.spacing = unit(0, "pt")) + 
  ggtitle("Fungal diversity belowground") +
  scale_shape_manual(values = shapes) + 
  guides(color = guide_legend(order = 1),
         fill = guide_legend(order = 1),
         shape = guide_legend(order = 2))

fig4c_wlegend
ggsave(paste0("atherton_figure4c_ITSdiversitybelowground_vs_dfb_wlegend_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig4c_wlegend, dpi = 300, width = 5.75, height = 5, 
       units = "in")

rm(list = c("metadata_soil_its", "metadata_soil_its_div", "model_all", 
            "model_nopath"))

### FIGURE 4D #################################################################
# Select the belowground data from 16s
metadata_osoil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 0-15 cm"),]
metadata_msoil_16s <- metadata_16s[which(metadata_16s$sample_type == "Soil, 15-30 cm"),]
metadata_root_16s <- metadata_16s[which(metadata_16s$sample_type %in% 
                                          c("Roots, 0-15 cm", 
                                            "Roots, 15-30 cm")),]

# Select the variables for the figure
metadata_osoil_16s_div <- metadata_osoil_16s[,colnames(metadata_osoil_16s) %in% 
                                               c("sample_name", "sample_type", 
                                                 "tree_dist_to_boston_km", 
                                                 "shannon.x", "shannon_nopath")]
metadata_msoil_16s_div <- metadata_msoil_16s[,colnames(metadata_msoil_16s) %in% 
                                               c("sample_name", "sample_type", 
                                                 "tree_dist_to_boston_km", 
                                                 "shannon.x", "shannon_nopath")]
metadata_root_16s_div <- metadata_root_16s[,colnames(metadata_root_16s) %in% 
                                             c("sample_name", "sample_type", 
                                               "tree_dist_to_boston_km", 
                                               "shannon.x", "shannon_nopath")]

# Format the data for plotting
diversity_long_osoil <- metadata_osoil_16s_div %>% 
  pivot_longer(!c(sample_name, sample_type, tree_dist_to_boston_km), 
               names_to = "diversity_type", values_to = "value")
diversity_long_osoil$diversity_type <- gsub("shannon_nopath", "non-pathogenic fungi",
                                            diversity_long_osoil$diversity_type)
diversity_long_osoil$diversity_type <- gsub("shannon.x", "all fungi",
                                            diversity_long_osoil$diversity_type)

diversity_long_msoil <- metadata_msoil_16s_div %>% 
  pivot_longer(!c(sample_name, sample_type, tree_dist_to_boston_km), 
               names_to = "diversity_type", values_to = "value")
diversity_long_msoil$diversity_type <- gsub("shannon_nopath", "non-pathogenic fungi",
                                            diversity_long_msoil$diversity_type)
diversity_long_msoil$diversity_type <- gsub("shannon.x", "all fungi",
                                            diversity_long_msoil$diversity_type)

diversity_long_root <- metadata_root_16s_div %>% 
  pivot_longer(!c(sample_name, sample_type, tree_dist_to_boston_km), 
               names_to = "diversity_type", values_to = "value")
diversity_long_root$diversity_type <- gsub("shannon_nopath", "non-pathogenic fungi",
                                           diversity_long_root$diversity_type)
diversity_long_root$diversity_type <- gsub("shannon.x", "all fungi",
                                           diversity_long_root$diversity_type)

### STATISTICS ALL FUNGI ###
model_all_osoil <- lme(formula(shannon.x ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
                       random = ~1|plot_name/tree_id, na.action = na.exclude, 
                       method = "ML", data = metadata_osoil_16s)
anova(model_all_osoil)
r.squaredGLMM(model_all_osoil)

model_all_msoil <- lme(formula(shannon.x ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
                       random = ~1|plot_name/tree_id, na.action = na.exclude, 
                       method = "ML", data = metadata_msoil_16s)
anova(model_all_msoil)
r.squaredGLMM(model_all_msoil)

model_all_root <- lme(formula(shannon.x ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
                      random = ~1|plot_name/tree_id, na.action = na.exclude, 
                      method = "ML", data = metadata_root_16s)
anova(model_all_root)
r.squaredGLMM(model_all_root)

### STATISTICS NON-PATHOGENIC FUNGI ###
model_nopath_osoil <- lme(formula(shannon_nopath ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
                          random = ~1|plot_name/tree_id, na.action = na.exclude, 
                          method = "ML", data = metadata_osoil_16s)
anova(model_nopath_osoil)
r.squaredGLMM(model_nopath_osoil)

model_nopath_msoil <- lme(formula(shannon_nopath ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
                          random = ~1|plot_name/tree_id, na.action = na.exclude, 
                          method = "ML", data = metadata_msoil_16s)
anova(model_nopath_msoil)
r.squaredGLMM(model_nopath_msoil)

model_nopath_root <- lme(formula(shannon_nopath ~ log(tree_dist_to_boston_km+1) * as.numeric(tree_distance_from_edge)), 
                         random = ~1|plot_name/tree_id, na.action = na.exclude, 
                         method = "ML", data = metadata_root_16s)
anova(model_nopath_root)
r.squaredGLMM(model_nopath_root)

# Merge all sample types into one
diversity_long_osoil$color <- paste0(diversity_long_osoil$diversity_type, " ", 
                                     diversity_long_osoil$sample_type)
diversity_long_msoil$color <- paste0(diversity_long_msoil$diversity_type, " ", 
                                     diversity_long_msoil$sample_type)
diversity_long_root$color <- paste0(diversity_long_root$diversity_type, " ", 
                                    diversity_long_root$sample_type)

diversity_long <- rbind(diversity_long_osoil, diversity_long_msoil, 
                        diversity_long_root)

new_pal <- c(categorical_pal[c(1,2,3,4)], linear_pal, "#808080", "#000000")

# Plot figure 4c
fig4d_osoil <- ggplot(diversity_long_osoil, aes(x = log(tree_dist_to_boston_km +1), 
                                                y = value, color = diversity_type,
                                                fill = diversity_type)) + 
  geom_point(aes(shape = sample_type), alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", alpha = 0.12, 
              aes(ymax = after_stat(y + se * sqrt(length(y))), 
                  ymin = after_stat(y - se * sqrt(length(y))))) + 
  scale_color_manual(values = linear_pal) + 
  scale_fill_manual(values = linear_pal) + 
  theme_classic() + 
  labs(x = "log(Distance from Boston (km))", y = "Shannon diversity", 
       legend = NA) + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), 
        legend.position = "none") + 
  ggtitle("Bacterial diversity belowground") +
  scale_shape_manual(values = shapes)

fig4d_osoil
ggsave(paste0("atherton_figure4c_16sdiversityosoil_vs_dfb_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig4d_osoil, dpi = 300, width = 5.75, height = 3.5, units = "in")

fig4d_msoil <- ggplot(diversity_long_msoil, aes(x = log(tree_dist_to_boston_km +1), 
                                                y = value, color = diversity_type,
                                                fill = diversity_type)) + 
  geom_point(aes(shape = sample_type), alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", alpha = 0.12, 
              aes(ymax = after_stat(y + se * sqrt(length(y))), 
                  ymin = after_stat(y - se * sqrt(length(y))))) + 
  scale_color_manual(values = linear_pal) + 
  scale_fill_manual(values = linear_pal) + 
  theme_classic() + 
  labs(x = "log(Distance from Boston (km))", y = "Shannon diversity", 
       legend = NA) + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), 
        legend.position = "none") + 
  ggtitle("Bacterial diversity belowground") +
  scale_shape_manual(values = shapes)

fig4d_msoil
ggsave(paste0("atherton_figure4c_16sdiversitymsoil_vs_dfb_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig4d_msoil, dpi = 300, width = 5.75, height = 3.5, units = "in")

fig4d_root <- ggplot(diversity_long_root, aes(x = log(tree_dist_to_boston_km +1), 
                                              y = value, color = diversity_type,
                                              fill = diversity_type)) + 
  geom_point(aes(shape = sample_type), alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", alpha = 0.12, 
              aes(ymax = after_stat(y + se * sqrt(length(y))), 
                  ymin = after_stat(y - se * sqrt(length(y))))) + 
  scale_color_manual(values = linear_pal) + 
  scale_fill_manual(values = linear_pal) + 
  theme_classic() + 
  labs(x = "log(Distance from Boston (km))", y = "Shannon diversity", 
       legend = NA) + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10), 
        legend.position = "none") + 
  ggtitle("Bacterial diversity belowground") +
  scale_shape_manual(values = shapes)

fig4d_root
ggsave(paste0("atherton_figure4c_16sdiversityroot_vs_dfb_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig4d_root, dpi = 300, width = 5.75, height = 3.5, units = "in")


fig4d_wlegend <- ggplot(diversity_long, aes(x = log(tree_dist_to_boston_km +1), 
                                            y = value, color = diversity_type, 
                                            fill = diversity_type)) + 
  geom_point(aes(shape = sample_type), alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", alpha = 0.12, 
              aes(ymax = after_stat(y + se * sqrt(length(y))), 
                  ymin = after_stat(y - se * sqrt(length(y))))) + 
  scale_color_manual(values = linear_pal) + 
  scale_fill_manual(values = linear_pal) + 
  theme_classic() + 
  labs(x = "log(Distance from Boston (km))", y = "Shannon diversity", shape = NA) + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_blank(), 
        plot.title = element_text(size = 10), 
        legend.position = "bottom", legend.box = "vertical", 
        legend.spacing = unit(0, "pt")) + 
  ggtitle("Bacterial diversity belowground") +
  scale_shape_manual(values = shapes) + 
  guides(color = guide_legend(order = 1),
         fill = guide_legend(order = 1),
         shape = guide_legend(order = 2))

fig4d_wlegend
ggsave(paste0("atherton_figure4d_16sdiversitybelowground_vs_dfb_wlegend_", 
              format(Sys.Date(), "%Y%m%d"), ".pdf"), device = cairo_pdf, 
       plot = fig4d_wlegend, dpi = 300, width = 5.75, height = 5, 
       units = "in")
