library(dplyr)

# load data 
leaf_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/leaf_its_genus.csv", row.names = 1)
root_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/root_its_genus.csv", row.names = 1)
m_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/m_its_genus.csv", row.names = 1)
o_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/o_its_genus.csv", row.names = 1)

leaf_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/leaf_16s_genus.csv", row.names = 1)
root_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/root_16s_genus.csv", row.names = 1)
m_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/m_16s_genus.csv", row.names = 1)
o_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Collapse_Genus/o_16s_genus.csv", row.names = 1)

metadata_its <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Data/metadata_its.csv")
metadata_16s <- read.csv("/projectnb/talbot-lab-data/Katies_data/Street_Trees/Data/metadata_16s_200.csv")

metadata_its <- metadata_its[!is.na(metadata_its$sample_name),]
metadata_16s <- metadata_16s[!is.na(metadata_16s$sample_name),]

metadata_its <- metadata_its[which(metadata_its$sequences_dropped == "No"),]
metadata_16s <- metadata_16s[which(metadata_16s$sequences_dropped == "No"),]

# get names for each tree type and sample type
leaf_its_aa01_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "AA01")]
leaf_its_ab_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "AB")]
leaf_its_bb_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "BB")]
leaf_its_bh02_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "BH02")]
leaf_its_bm03_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "BM03")]
leaf_its_do_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "DO")]
leaf_its_hf04_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "HF04")]
leaf_its_hf06_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "HF06")]
leaf_its_hf069c5_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "C5")]
leaf_its_hf069c2_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "C2")]
leaf_its_hf069d3_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "D3")]
leaf_its_hf069d5_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "D5")]
leaf_its_hf069e3_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "E3")]
leaf_its_hw07_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "HW07")]
leaf_its_jp_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "JP")]
leaf_its_mc_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "MC")]
leaf_its_mp_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "MP")]
leaf_its_sb_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "SB")]
leaf_its_se_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "SE")]
leaf_its_sw08_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "SW08")]
leaf_its_wr_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Leaf" & metadata_its$plot_name == "WR")]

root_its_aa01_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "AA01")]
root_its_ab_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "AB")]
root_its_bb_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "BB")]
root_its_bh02_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "BH02")]
root_its_bm03_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "BM03")]
root_its_do_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "DO")]
root_its_hf04_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "HF04")]
root_its_hf06_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "HF06")]
root_its_hf069c5_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "C5")]
root_its_hf069c2_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "C2")]
root_its_hf069d3_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "D3")]
root_its_hf069d5_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "D5")]
root_its_hf069e3_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "E3")]
root_its_hw07_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "HW07")]
root_its_jp_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "JP")]
root_its_mc_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "MC")]
root_its_mp_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "MP")]
root_its_sb_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "SB")]
root_its_se_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "SE")]
root_its_sw08_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "SW08")]
root_its_wr_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Root" & metadata_its$plot_name == "WR")]

m_its_aa01_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "AA01")]
m_its_ab_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "AB")]
m_its_bb_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "BB")]
m_its_bh02_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "BH02")]
m_its_bm03_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "BM03")]
m_its_do_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "DO")]
m_its_hf04_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "HF04")]
m_its_hf06_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "HF06")]
m_its_hf069c5_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "C5")]
m_its_hf069c2_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "C2")]
m_its_hf069d3_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "D3")]
m_its_hf069d5_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "D5")]
m_its_hf069e3_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "E3")]
m_its_hw07_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "HW07")]
m_its_jp_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "JP")]
m_its_mc_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "MC")]
m_its_mp_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "MP")]
m_its_sb_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "SB")]
m_its_se_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "SE")]
m_its_sw08_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "SW08")]
m_its_wr_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "M" & metadata_its$plot_name == "WR")]

o_its_aa01_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "AA01")]
o_its_ab_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "AB")]
o_its_bb_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "BB")]
o_its_bh02_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "BH02")]
o_its_bm03_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "BM03")]
o_its_do_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "DO")]
o_its_hf04_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "HF04")]
o_its_hf06_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "HF06" & metadata_its$tree_type == "Forest Edge")]
o_its_hf069c5_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "C5")]
o_its_hf069c2_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "C2")]
o_its_hf069d3_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "D3")]
o_its_hf069d5_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "D5")]
o_its_hf069e3_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "E3")]
o_its_hw07_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "HW07")]
o_its_jp_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "JP")]
o_its_mc_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "MC")]
o_its_mp_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "MP")]
o_its_sb_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "SB")]
o_its_se_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "SE")]
o_its_sw08_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "SW08")]
o_its_wr_names <- metadata_its$sample_name[which(metadata_its$sample_type == "Soil" & metadata_its$soil_horizon == "O" & metadata_its$plot_name == "WR")]

leaf_16s_aa01_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "AA01")]
leaf_16s_ab_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "AB")]
leaf_16s_bb_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "BB")]
leaf_16s_bh02_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "BH02")]
leaf_16s_bm03_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "BM03")]
leaf_16s_do_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "DO")]
leaf_16s_hf04_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "HF04")]
leaf_16s_hf06_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "HF06")]
leaf_16s_hf069c5_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "C5")]
leaf_16s_hf069c2_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "C2")]
leaf_16s_hf069d3_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "D3")]
leaf_16s_hf069d5_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "D5")]
leaf_16s_hf069e3_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "E3")]
leaf_16s_hw07_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "HW07")]
leaf_16s_jp_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "JP")]
leaf_16s_mc_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "MC")]
leaf_16s_mp_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "MP")]
leaf_16s_sb_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "SB")]
leaf_16s_se_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "SE")]
leaf_16s_sw08_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "SW08")]
leaf_16s_wr_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Leaf" & metadata_16s$plot_name == "WR")]

root_16s_aa01_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "AA01")]
root_16s_ab_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "AB")]
root_16s_bb_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "BB")]
root_16s_bh02_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "BH02")]
root_16s_bm03_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "BM03")]
root_16s_do_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "DO")]
root_16s_hf04_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "HF04")]
root_16s_hf06_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "HF06")]
root_16s_hf069c5_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "C5")]
root_16s_hf069c2_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "C2")]
root_16s_hf069d3_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "D3")]
root_16s_hf069d5_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "D5")]
root_16s_hf069e3_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "E3")]
root_16s_hw07_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "HW07")]
root_16s_jp_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "JP")]
root_16s_mc_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "MC")]
root_16s_mp_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "MP")]
root_16s_sb_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "SB")]
root_16s_se_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "SE")]
root_16s_sw08_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "SW08")]
root_16s_wr_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Root" & metadata_16s$plot_name == "WR")]

m_16s_aa01_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "AA01")]
m_16s_ab_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "AB")]
m_16s_bb_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "BB")]
m_16s_bh02_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "BH02")]
m_16s_bm03_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "BM03")]
m_16s_do_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "DO")]
m_16s_hf04_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "HF04")]
m_16s_hf06_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "HF06")]
m_16s_hf069c5_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "C5")]
m_16s_hf069c2_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "C2")]
m_16s_hf069d3_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "D3")]
m_16s_hf069d5_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "D5")]
m_16s_hf069e3_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "E3")]
m_16s_hw07_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "HW07")]
m_16s_jp_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "JP")]
m_16s_mc_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "MC")]
m_16s_mp_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "MP")]
m_16s_sb_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "SB")]
m_16s_se_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "SE")]
m_16s_sw08_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "SW08")]
m_16s_wr_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "M" & metadata_16s$plot_name == "WR")]

o_16s_aa01_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "AA01")]
o_16s_ab_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "AB")]
o_16s_bb_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "BB")]
o_16s_bh02_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "BH02")]
o_16s_bm03_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "BM03")]
o_16s_do_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "DO")]
o_16s_hf04_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "HF04")]
o_16s_hf06_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "HF06")]
o_16s_hf069c5_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "C5")]
o_16s_hf069c2_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "C2")]
o_16s_hf069d3_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "D3")]
o_16s_hf069d5_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "D5")]
o_16s_hf069e3_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "E3")]
o_16s_hw07_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "HW07")]
o_16s_jp_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "JP")]
o_16s_mc_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "MC")]
o_16s_mp_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "MP")]
o_16s_sb_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "SB")]
o_16s_se_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "SE")]
o_16s_sw08_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "SW08")]
o_16s_wr_names <- metadata_16s$sample_name[which(metadata_16s$sample_type == "Soil" & metadata_16s$soil_horizon == "O" & metadata_16s$plot_name == "WR")]

# make matrices
leaf_its_aa01 <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_aa01_names)]
leaf_its_abbb <- leaf_its[,which(colnames(leaf_its) %in% c(leaf_its_ab_names, leaf_its_bb_names))]
leaf_its_bm03 <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_bm03_names)]
leaf_its_dosb <- leaf_its[,which(colnames(leaf_its) %in% c(leaf_its_do_names, leaf_its_sb_names))]
leaf_its_hf04 <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_hf04_names)]
leaf_its_hf06 <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_hf06_names)]
leaf_its_hf069_c <- leaf_its[,which(colnames(leaf_its) %in% c(leaf_its_hf069c5_names, leaf_its_hf069c2_names))]
leaf_its_hf069_d <- leaf_its[,which(colnames(leaf_its) %in% c(leaf_its_hf069d3_names, leaf_its_hf069d5_names))]
leaf_its_hf069_e <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_hf069e3_names)]
leaf_its_hw07 <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_hw07_names)]
leaf_its_jp <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_jp_names)]
leaf_its_mcse <- leaf_its[,which(colnames(leaf_its) %in% c(leaf_its_mc_names, leaf_its_se_names))]
leaf_its_mp <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_mp_names)]
leaf_its_sw08 <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_sw08_names)]
leaf_its_wr <- leaf_its[,which(colnames(leaf_its) %in% leaf_its_wr_names)]

root_its_aa01 <- root_its[,which(colnames(root_its) %in% root_its_aa01_names)]
root_its_abbb <- root_its[,which(colnames(root_its) %in% c(root_its_ab_names, root_its_bb_names))]
root_its_bm03 <- root_its[,which(colnames(root_its) %in% root_its_bm03_names)]
root_its_dosb <- root_its[,which(colnames(root_its) %in% c(root_its_do_names, root_its_sb_names))]
root_its_hf04 <- root_its[,which(colnames(root_its) %in% root_its_hf04_names)]
root_its_hf06 <- root_its[,which(colnames(root_its) %in% root_its_hf06_names)]
root_its_hf069_c <- root_its[,which(colnames(root_its) %in% c(root_its_hf069c5_names, root_its_hf069c2_names))]
root_its_hf069_d <- root_its[,which(colnames(root_its) %in% c(root_its_hf069d3_names, root_its_hf069d5_names))]
root_its_hf069_e <- root_its[,which(colnames(root_its) %in% root_its_hf069e3_names)]
root_its_hw07 <- root_its[,which(colnames(root_its) %in% root_its_hw07_names)]
root_its_jp <- root_its[,which(colnames(root_its) %in% root_its_jp_names)]
root_its_mcse <- root_its[,which(colnames(root_its) %in% c(root_its_mc_names, root_its_se_names))]
root_its_mp <- root_its[,which(colnames(root_its) %in% root_its_mp_names)]
root_its_sw08 <- root_its[,which(colnames(root_its) %in% root_its_sw08_names)]
root_its_wr <- root_its[,which(colnames(root_its) %in% root_its_wr_names)]

m_its_aa01 <- m_its[,which(colnames(m_its) %in% m_its_aa01_names)]
m_its_abbb <- m_its[,which(colnames(m_its) %in% c(m_its_ab_names, m_its_bb_names))]
m_its_bm03 <- m_its[,which(colnames(m_its) %in% m_its_bm03_names)]
m_its_dosb <- m_its[,which(colnames(m_its) %in% c(m_its_do_names, m_its_sb_names))]
m_its_hf04 <- m_its[,which(colnames(m_its) %in% m_its_hf04_names)]
m_its_hf06 <- m_its[,which(colnames(m_its) %in% m_its_hf06_names)]
m_its_hf069_c <- m_its[,which(colnames(m_its) %in% c(m_its_hf069c5_names, m_its_hf069c2_names))]
m_its_hf069_d <- m_its[,which(colnames(m_its) %in% c(m_its_hf069d3_names, m_its_hf069d5_names))]
m_its_hf069_e <- m_its[,which(colnames(m_its) %in% m_its_hf069e3_names)]
m_its_hw07 <- m_its[,which(colnames(m_its) %in% m_its_hw07_names)]
m_its_jp <- m_its[,which(colnames(m_its) %in% m_its_jp_names)]
m_its_mcse <- m_its[,which(colnames(m_its) %in% c(m_its_mc_names, m_its_se_names))]
m_its_mp <- m_its[,which(colnames(m_its) %in% m_its_mp_names)]
m_its_sw08 <- m_its[,which(colnames(m_its) %in% m_its_sw08_names)]
m_its_wr <- m_its[,which(colnames(m_its) %in% m_its_wr_names)]

o_its_aa01 <- o_its[,which(colnames(o_its) %in% o_its_aa01_names)]
o_its_abbb <- o_its[,which(colnames(o_its) %in% c(o_its_ab_names, o_its_bb_names))]
o_its_bm03 <- o_its[,which(colnames(o_its) %in% o_its_bm03_names)]
o_its_dosb <- o_its[,which(colnames(o_its) %in% c(o_its_do_names, o_its_sb_names))]
o_its_hf04 <- o_its[,which(colnames(o_its) %in% o_its_hf04_names)]
o_its_hf06 <- o_its[,which(colnames(o_its) %in% o_its_hf06_names)]
o_its_hf069_c <- o_its[,which(colnames(o_its) %in% c(o_its_hf069c5_names, o_its_hf069c2_names))]
o_its_hf069_d <- o_its[,which(colnames(o_its) %in% c(o_its_hf069d3_names, o_its_hf069d5_names))]
o_its_hf069_e <- o_its[,which(colnames(o_its) %in% o_its_hf069e3_names)]
o_its_hw07 <- o_its[,which(colnames(o_its) %in% o_its_hw07_names)]
o_its_jp <- o_its[,which(colnames(o_its) %in% o_its_jp_names)]
o_its_mcse <- o_its[,which(colnames(o_its) %in% c(o_its_mc_names, o_its_se_names))]
o_its_mp <- o_its[,which(colnames(o_its) %in% o_its_mp_names)]
o_its_sw08 <- o_its[,which(colnames(o_its) %in% o_its_sw08_names)]
o_its_wr <- o_its[,which(colnames(o_its) %in% o_its_wr_names)]

leaf_16s_aa01 <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_aa01_names)]
leaf_16s_abbb <- leaf_16s[,which(colnames(leaf_16s) %in% c(leaf_16s_ab_names, leaf_16s_bb_names))]
leaf_16s_bm03 <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_bm03_names)]
leaf_16s_dosb <- leaf_16s[,which(colnames(leaf_16s) %in% c(leaf_16s_do_names, leaf_16s_sb_names))]
leaf_16s_hf04 <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_hf04_names)]
leaf_16s_hf06 <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_hf06_names)]
leaf_16s_hf069_c <- leaf_16s[,which(colnames(leaf_16s) %in% c(leaf_16s_hf069c5_names, leaf_16s_hf069c2_names))]
leaf_16s_hf069_d <- leaf_16s[,which(colnames(leaf_16s) %in% c(leaf_16s_hf069d3_names, leaf_16s_hf069d5_names))]
leaf_16s_hf069_e <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_hf069e3_names)]
leaf_16s_hw07 <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_hw07_names)]
leaf_16s_jp <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_jp_names)]
leaf_16s_mcse <- leaf_16s[,which(colnames(leaf_16s) %in% c(leaf_16s_mc_names, leaf_16s_se_names))]
leaf_16s_mp <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_mp_names)]
leaf_16s_sw08 <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_sw08_names)]
leaf_16s_wr <- leaf_16s[,which(colnames(leaf_16s) %in% leaf_16s_wr_names)]

root_16s_aa01 <- root_16s[,which(colnames(root_16s) %in% root_16s_aa01_names)]
root_16s_abbb <- root_16s[,which(colnames(root_16s) %in% c(root_16s_ab_names, root_16s_bb_names))]
root_16s_bm03 <- root_16s[,which(colnames(root_16s) %in% root_16s_bm03_names)]
root_16s_dosb <- root_16s[,which(colnames(root_16s) %in% c(root_16s_do_names, root_16s_sb_names))]
root_16s_hf04 <- root_16s[,which(colnames(root_16s) %in% root_16s_hf04_names)]
root_16s_hf06 <- root_16s[,which(colnames(root_16s) %in% root_16s_hf06_names)]
root_16s_hf069_c <- root_16s[,which(colnames(root_16s) %in% c(root_16s_hf069c5_names, root_16s_hf069c2_names))]
root_16s_hf069_d <- root_16s[,which(colnames(root_16s) %in% c(root_16s_hf069d3_names, root_16s_hf069d5_names))]
root_16s_hf069_e <- root_16s[,which(colnames(root_16s) %in% root_16s_hf069e3_names)]
root_16s_hw07 <- root_16s[,which(colnames(root_16s) %in% root_16s_hw07_names)]
root_16s_jp <- root_16s[,which(colnames(root_16s) %in% root_16s_jp_names)]
root_16s_mcse <- root_16s[,which(colnames(root_16s) %in% c(root_16s_mc_names, root_16s_se_names))]
root_16s_mp <- root_16s[,which(colnames(root_16s) %in% root_16s_mp_names)]
root_16s_sw08 <- root_16s[,which(colnames(root_16s) %in% root_16s_sw08_names)]
root_16s_wr <- root_16s[,which(colnames(root_16s) %in% root_16s_wr_names)]

m_16s_aa01 <- m_16s[,which(colnames(m_16s) %in% m_16s_aa01_names)]
m_16s_abbb <- m_16s[,which(colnames(m_16s) %in% c(m_16s_ab_names, m_16s_bb_names))]
m_16s_bm03 <- m_16s[,which(colnames(m_16s) %in% m_16s_bm03_names)]
m_16s_dosb <- m_16s[,which(colnames(m_16s) %in% c(m_16s_do_names, m_16s_sb_names))]
m_16s_hf04 <- m_16s[,which(colnames(m_16s) %in% m_16s_hf04_names)]
m_16s_hf06 <- m_16s[,which(colnames(m_16s) %in% m_16s_hf06_names)]
m_16s_hf069_c <- m_16s[,which(colnames(m_16s) %in% c(m_16s_hf069c5_names, m_16s_hf069c2_names))]
m_16s_hf069_d <- m_16s[,which(colnames(m_16s) %in% c(m_16s_hf069d3_names, m_16s_hf069d5_names))]
m_16s_hf069_e <- m_16s[,which(colnames(m_16s) %in% m_16s_hf069e3_names)]
m_16s_hw07 <- m_16s[,which(colnames(m_16s) %in% m_16s_hw07_names)]
m_16s_jp <- m_16s[,which(colnames(m_16s) %in% m_16s_jp_names)]
m_16s_mcse <- m_16s[,which(colnames(m_16s) %in% c(m_16s_mc_names, m_16s_se_names))]
m_16s_mp <- m_16s[,which(colnames(m_16s) %in% m_16s_mp_names)]
m_16s_sw08 <- m_16s[,which(colnames(m_16s) %in% m_16s_sw08_names)]
m_16s_wr <- m_16s[,which(colnames(m_16s) %in% m_16s_wr_names)]

o_16s_aa01 <- o_16s[,which(colnames(o_16s) %in% o_16s_aa01_names)]
o_16s_abbb <- o_16s[,which(colnames(o_16s) %in% c(o_16s_ab_names, o_16s_bb_names))]
o_16s_bm03 <- o_16s[,which(colnames(o_16s) %in% o_16s_bm03_names)]
o_16s_dosb <- o_16s[,which(colnames(o_16s) %in% c(o_16s_do_names, o_16s_sb_names))]
o_16s_hf04 <- o_16s[,which(colnames(o_16s) %in% o_16s_hf04_names)]
o_16s_hf06 <- o_16s[,which(colnames(o_16s) %in% o_16s_hf06_names)]
o_16s_hf069_c <- o_16s[,which(colnames(o_16s) %in% c(o_16s_hf069c5_names, o_16s_hf069c2_names))]
o_16s_hf069_d <- o_16s[,which(colnames(o_16s) %in% c(o_16s_hf069d3_names, o_16s_hf069d5_names))]
o_16s_hf069_e <- o_16s[,which(colnames(o_16s) %in% o_16s_hf069e3_names)]
o_16s_hw07 <- o_16s[,which(colnames(o_16s) %in% o_16s_hw07_names)]
o_16s_jp <- o_16s[,which(colnames(o_16s) %in% o_16s_jp_names)]
o_16s_mcse <- o_16s[,which(colnames(o_16s) %in% c(o_16s_mc_names, o_16s_se_names))]
o_16s_mp <- o_16s[,which(colnames(o_16s) %in% o_16s_mp_names)]
o_16s_sw08 <- o_16s[,which(colnames(o_16s) %in% o_16s_sw08_names)]
o_16s_wr <- o_16s[,which(colnames(o_16s) %in% o_16s_wr_names)]

# make lists
leaf_its <- list(leaf_its_aa01, leaf_its_bm03, leaf_its_hf04, leaf_its_hf06, 
                 leaf_its_hf069_c, leaf_its_hf069_d, leaf_its_hw07, leaf_its_jp, 
                 leaf_its_sw08)

root_its <- list(root_its_aa01, root_its_bm03, root_its_dosb, root_its_hf06, 
                 root_its_hf069_c, root_its_hf069_d, root_its_hw07, root_its_jp, 
                 root_its_mcse, root_its_sw08, root_its_wr)

m_its <- list(m_its_aa01, m_its_abbb, m_its_bm03, m_its_dosb, m_its_hf04, 
              m_its_hf06, m_its_hf069_c, m_its_hf069_d, m_its_hf069_e, m_its_hw07, 
              m_its_jp, m_its_mcse, m_its_mp, m_its_wr)

o_its <- list(o_its_aa01, o_its_abbb, o_its_bm03, o_its_dosb, o_its_hf04, 
              o_its_hf069_c, o_its_hf069_d, o_its_hf069_e, o_its_hw07, o_its_jp, 
              o_its_mcse, o_its_mp, o_its_sw08, o_its_wr)

leaf_16s <- list(leaf_16s_aa01, leaf_16s_bm03, leaf_16s_hf04, leaf_16s_hf06, 
                 leaf_16s_hf069_c, leaf_16s_hf069_d, leaf_16s_hw07, leaf_16s_jp, 
                 leaf_16s_sw08)

root_16s <- list(root_16s_aa01, root_16s_bm03, root_16s_dosb, root_16s_hf06, 
                 root_16s_hf069_c, root_16s_hf069_d, root_16s_hw07, root_16s_jp, 
                 root_16s_mcse, root_16s_sw08, root_16s_wr)

m_16s <- list(m_16s_aa01, m_16s_abbb, m_16s_bm03, m_16s_dosb, m_16s_hf04, 
              m_16s_hf06, m_16s_hf069_c, m_16s_hf069_d, m_16s_hf069_e, m_16s_hw07, 
              m_16s_jp, m_16s_mcse, m_16s_mp, m_16s_wr)

o_16s <- list(o_16s_aa01, o_16s_abbb, o_16s_bm03, o_16s_dosb, o_16s_hf04, 
              o_16s_hf069_c, o_16s_hf069_d, o_16s_hf069_e, o_16s_hw07, o_16s_jp, 
              o_16s_mcse, o_16s_mp, o_16s_sw08, o_16s_wr)

# filter out genera with < 25 reads
for(i in 1:length(leaf_its)){
  x <- leaf_its[[i]]
  x <- x[which(rowSums(x) > 25),]
  leaf_its[[i]] <- x
}

for(i in 1:length(root_its)){
  x <- root_its[[i]]
  x <- x[which(rowSums(x) > 25),]
  root_its[[i]] <- x
}

for(i in 1:length(m_its)){
  x <- m_its[[i]]
  x <- x[which(rowSums(x) > 25),]
  m_its[[i]] <- x
}

for(i in 1:length(o_its)){
  x <- o_its[[i]]
  x <- x[which(rowSums(x) > 25),]
  o_its[[i]] <- x
}

for(i in 1:length(leaf_16s)){
  x <- leaf_16s[[i]]
  x <- x[which(rowSums(x) > 25),]
  leaf_16s[[i]] <- x
}

for(i in 1:length(root_16s)){
  x <- root_16s[[i]]
  x <- x[which(rowSums(x) > 25),]
  root_16s[[i]] <- x
}

for(i in 1:length(m_16s)){
  x <- m_16s[[i]]
  x <- x[which(rowSums(x) > 25),]
  m_16s[[i]] <- x
}

for(i in 1:length(o_16s)){
  x <- o_16s[[i]]
  x <- x[which(rowSums(x) > 25),]
  o_16s[[i]] <- x
}

# filter out genera present in < 20% of samples in a network
for(i in 1:length(leaf_its)){
  x <- leaf_its[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  leaf_its[[i]] <- x
}

for(i in 1:length(root_its)){
  x <- root_its[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  root_its[[i]] <- x
}

for(i in 1:length(m_its)){
  x <- m_its[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  m_its[[i]] <- x
}

for(i in 1:length(o_its)){
  x <- o_its[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  o_its[[i]] <- x
}

for(i in 1:length(leaf_16s)){
  x <- leaf_16s[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  leaf_16s[[i]] <- x
}

for(i in 1:length(root_16s)){
  x <- root_16s[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  root_16s[[i]] <- x
}

for(i in 1:length(m_16s)){
  x <- m_16s[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  m_16s[[i]] <- x
}

for(i in 1:length(o_16s)){
  x <- o_16s[[i]]
  x <- x[which(rowSums(x!=0) > (0.2 * ncol(x))),]
  o_16s[[i]] <- x
}

# filter out genera < 0.1% abundance
for(i in 1:length(leaf_its)){
  x <- leaf_its[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  leaf_its[[i]] <- x
}

for(i in 1:length(root_its)){
  x <- root_its[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  root_its[[i]] <- x
}

for(i in 1:length(m_its)){
  x <- m_its[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  m_its[[i]] <- x
}

for(i in 1:length(o_its)){
  x <- o_its[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  o_its[[i]] <- x
}

for(i in 1:length(leaf_16s)){
  x <- leaf_16s[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  leaf_16s[[i]] <- x
}

for(i in 1:length(root_16s)){
  x <- root_16s[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  root_16s[[i]] <- x
}

for(i in 1:length(m_16s)){
  x <- m_16s[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  m_16s[[i]] <- x
}

for(i in 1:length(o_16s)){
  x <- o_16s[[i]]
  x <- x[which(rowSums(x) > (0.001 * sum(rowSums(x)))),]
  o_16s[[i]] <- x
}

# save tables
for(i in 1:length(leaf_its)){
  rownames(leaf_its[[i]]) <- gsub("g__", "f__", rownames(leaf_its[[i]]))
  colnames(leaf_its[[i]]) <- gsub("_ITS.*", "", colnames(leaf_its[[i]]))
}

for(i in 1:length(root_its)){
  rownames(root_its[[i]]) <- gsub("g__", "f__", rownames(root_its[[i]]))
  colnames(root_its[[i]]) <- gsub("_ITS.*", "", colnames(root_its[[i]]))
}

for(i in 1:length(m_its)){
  rownames(m_its[[i]]) <- gsub("g__", "f__", rownames(m_its[[i]]))
  colnames(m_its[[i]]) <- gsub("_ITS.*", "", colnames(m_its[[i]]))
}

for(i in 1:length(o_its)){
  rownames(o_its[[i]]) <- gsub("g__", "f__", rownames(o_its[[i]]))
  colnames(o_its[[i]]) <- gsub("_ITS.*", "", colnames(o_its[[i]]))
}

for(i in 1:length(leaf_16s)){
  rownames(leaf_16s[[i]]) <- gsub("g__", "b__", rownames(leaf_16s[[i]]))
  colnames(leaf_16s[[i]]) <- gsub("_16S.*", "", colnames(leaf_16s[[i]]))
}

for(i in 1:length(root_16s)){
  rownames(root_16s[[i]]) <- gsub("g__", "b__", rownames(root_16s[[i]]))
  colnames(root_16s[[i]]) <- gsub("_16S.*", "", colnames(root_16s[[i]]))
}

for(i in 1:length(m_16s)){
  rownames(m_16s[[i]]) <- gsub("g__", "b__", rownames(m_16s[[i]]))
  colnames(m_16s[[i]]) <- gsub("_16S.*", "", colnames(m_16s[[i]]))
}

for(i in 1:length(o_16s)){
  rownames(o_16s[[i]]) <- gsub("g__", "b__", rownames(o_16s[[i]]))
  colnames(o_16s[[i]]) <- gsub("_16S.*", "", colnames(o_16s[[i]]))
}

saveRDS(leaf_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/leaf_its.RDS")
saveRDS(root_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/root_its.RDS")
saveRDS(m_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/m_its.RDS")
saveRDS(o_its, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/o_its.RDS")

saveRDS(leaf_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/leaf_16s.RDS")
saveRDS(root_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/root_16s.RDS")
saveRDS(m_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/m_16s.RDS")
saveRDS(o_16s, "/projectnb/talbot-lab-data/Katies_data/Street_Trees/Analysis/Networks/Data/Filtered/o_16s.RDS")