#======== PROJECT COM2LIFE ========
## Stackbarplots for species, help of L.Ezzat 

# Load libraries
library(phyloseq)
library(phyloseqCompanion)
library(tidyverse)
library(dplyr)
library(ggplot2)

# Load physeq object
physeq <- readRDS(here::here("Data",
                             "mon_objet_physeq_17_06.rds"))
# Filter physeq to keep only samples of fish
physeq_filtered <- subset_samples(physeq, !(grepl("BOI|JAB|VAI|SED|ADN", sample_names(physeq))))
rows_to_remove <- grepl("BOI|JAB|VAI|SED|ADN", rownames(sample_data(physeq_filtered)))

# keep the sam data
physeq_filtered@sam_data <- sample_data(physeq_filtered)[!rows_to_remove, ]

# Modify lakes's name
# Access samples name
sample_names_current <- sample_names(physeq_filtered)
# Modify names
sample_names_modified <- sample_names_current %>%
  gsub("CHA", "CSM", .) %>%
  gsub("CRJ", "CERL", .) %>% 
  gsub("TRI", "CERS", .) %>%
  gsub("CTL", "CRE", .) %>%
  gsub("GDP", "LGP", .) %>%
  gsub("VER", "VSS", .)
# Apply new names
sample_names(physeq_filtered) <- sample_names_modified
# Modify names for sam_data
sample_data <- sample_data(physeq_filtered)
# for lake 
if ("lake" %in% colnames(sample_data)) {
  sample_data[["lake"]] <- sample_data[["lake"]] %>%
    gsub("CHA", "CSM", .) %>%
    gsub("CRJ1", "CERL", .) %>% 
    gsub("CRJ2", "CERS", .) %>%
    gsub("CTL", "CRE", .) %>%
    gsub("GDP", "LGP", .) %>%
    gsub("VER", "VSS", .)
}
# for lake_origin
if ("lake_origin" %in% colnames(sample_data)) {
  sample_data[["lake_origin"]] <- sample_data[["lake_origin"]] %>%
    gsub("CHA", "CSM", .) %>%
    gsub("CRJ1", "CERL", .) %>% 
    gsub("CRJ2", "CERS", .) %>%
    gsub("CTL", "CRE", .) %>%
    gsub("GDP", "LGP", .) %>%
    gsub("VER", "VSS", .)
}
#for samples 
if ("samples" %in% colnames(sample_data)) {
  sample_data[["samples"]] <- sample_data[["samples"]] %>%
    gsub("CHA", "CSM", .) %>%
    gsub("CRJ", "CERL", .) %>% 
    gsub("TRI", "CERS", .) %>%
    gsub("CTL", "CRE", .) %>%
    gsub("GDP", "LGP", .) %>%
    gsub("VER", "VSS", .)
}
# Update sample_data in phyloseq object
sample_data(physeq_filtered) <- sample_data

#Save
path_to_my_object = here::here("Data","mon_objet_physeq_filtered.rds")
saveRDS(physeq_filtered, file = path_to_my_object)

metadata <- phyloseqCompanion::sample.data.frame(physeq_filtered)

# ====== PHYLUM LEVEL =======
# Most abundant phyla
phyla_sp <- tax_glom(physeq_filtered, taxrank = rank_names(physeq_filtered)[2], NArm = F) # agglomérate to phylum level
phyla_sp <- transform_sample_counts(phyla_sp, function(x) x / sum(x)) # relative abundance
TopASV_p <- names(sort(taxa_sums(phyla_sp), TRUE)[1:10]) # most abundant phylum
top10_water_p <- prune_species(TopASV_p, phyla_sp)
top10_water_p <- prune_taxa(taxa_sums(top10_water_p) > 0, top10_water_p)
top_phyla <- as.data.frame(tax_table(top10_water_p))

# seperate data by lakes
physeq_merge_by <- merge_samples(phyla_sp, "lake_origin")
sp_RA <- transform_sample_counts(physeq_merge_by, function(x) x / sum(x))
sp_melt <- psmelt(sp_RA) # create a datframe from phyloseq object
sp_melt$Phylum <- as.character(sp_melt$Phylum) # convert into character

# Keep 10 most abundant phylum
sumtot_sp <- sp_melt %>%
  dplyr::group_by(Phylum) %>%
  dplyr::summarize(sum = sum(Abundance)) %>%
  filter(Phylum %in% top_phyla$Phylum) %>%
  filter(!(Phylum %in% c("", " p__uncultured", NA, "Unknown Phylum")))

sp_melt$Phylum[!(sp_melt$Phylum %in% sumtot_sp$Phylum)] <- "Other"

sp_melt$totalAbundance <- sum(sp_melt$Abundance)

data_sp_mod <- sp_melt %>%
  group_by(Phylum, Sample) %>%
  summarise(Abundance = sum(Abundance)) %>%
  dplyr::distinct()

# Specify order of samples
data_sp_mod$Sample <- factor(data_sp_mod$Sample, levels = c("CERL_PER", "CERL_LEP", "CERS_PER","CERS_LEP", 
                                                            "CRE_PER", "CRE_LEP", "VSS_PER","VSS_LEP",
                                                            "LGP_PER", "LGP_LEP", "CSM_PER", "CSM_LEP"))

color_palette <- c(
  "#3C5488",
  "#4DBBD5",
  "#B0E0E6",
  "#8491B4",
  "#cab2d6",
  "#AD002A",
  "lightgrey",
  "#F39B7F",
  "#fdbf6f",
  "#B09C85"
)

# create barplot
barplot_for_sp <- ggplot2::ggplot(data = data_sp_mod, ggplot2::aes(x = Sample, y = Abundance, fill = Phylum)) +
  ggplot2::geom_bar(stat = "identity", position = "stack") +
  ggplot2::scale_fill_manual(values = color_palette) +
  ggplot2::theme_bw() +
  ggplot2::labs(y = "Relative Abundance") +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 12),     
    axis.title.y = ggplot2::element_text(size = 12),     
    axis.text.x = ggplot2::element_text(size = 10, angle = 45, hjust = 1), 
    axis.text.y = ggplot2::element_text(size = 10),      
    legend.title = ggplot2::element_text(size = 14),     
    legend.text = ggplot2::element_text(size = 12)       
  ) +
  ggplot2::theme(panel.border = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "black"))

barplot_for_sp

path_to_my_object = here::here("figures", "figure1.png")
ggplot2::ggsave(filename = path_to_my_object, plot = barplot_for_sp, device = "png")
