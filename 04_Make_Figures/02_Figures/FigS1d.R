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
setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/02_Clean_Data/05_Transform_Data/16S/Aitchison_Distance_Tables/")
bg_aitch_16s <- read.csv("atherton_16S_allsampletypes_aitchisondistance_20240925.csv", 
                         row.names = 1)
leaf_aitch_16s <- read.csv("atherton_16S_leaf_aitchisondistance_20241112.csv", 
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

# Normalize the data to be used in this figure
metadata_its$perc_wood_saprotroph <- log(metadata_its$perc_wood_saprotroph+1)

setwd("/projectnb/talbot-lab-data/Katies_data/M-BUDS/04_Make_Figures/02_Figures/")

### FIGURE S1C ##################################################################
# Select belowground data from ITS
metadata_soil_its <- metadata_its[which(metadata_its$sample_type != "Leaf"),]

### STATISTICS ECM VS WOOD DECOMPOSERS ###
model_ecm_wood <- lme(formula(perc_ecm ~ perc_wood_saprotroph), 
                      random = ~1|plot_name/tree_id, na.action = na.exclude, 
                      method = "ML", data = metadata_soil_its)
anova(model_ecm_wood)
r.squaredGLMM(model_ecm_wood)

figS1c <- ggplot(metadata_soil_its, aes(x = perc_ecm, 
                                        y = perc_wood_saprotroph)) + 
  geom_point(alpha = 0.5, col = "#003049", aes(shape = sample_type)) +
  xlab("ECM relative abundance (%)))") + 
  ylab("log(Wood decomposer\nrelative abundance (%))") +
  labs(shape = "") +
  geom_smooth(method = "lm", alpha = 0.12, col = "#003049", fill = "#003049", 
              aes(x = perc_ecm,
                  y = perc_wood_saprotroph,
                  col = Guild, fill = Guild,
                  ymax = after_stat(y + se * sqrt(length(y))),
                  ymin = after_stat(y - se * sqrt(length(y))))) +
  theme_classic() + 
  theme(text = element_text(family = "Arial"),
        axis.title.x = element_text(size = 10, color="black"), 
        axis.title.y = element_text(size = 10), 
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 1), 
        plot.title = element_text(size = 10),
        legend.position = "bottom") + 
  scale_shape_manual(values = shapes[2:5]) +
  ggtitle("Fungal wood decomposer and ECM relative abundance belowground")

figS1c
ggsave(paste0("atherton_figureS1c_wooddecomposersroots_vs_ecm_",
              format(Sys.Date(), "%Y%m%d"),".pdf"), 
       plot = figS1c, dpi = 300, width = 6, height = 3.5, units = "in", 
       device = cairo_pdf)
