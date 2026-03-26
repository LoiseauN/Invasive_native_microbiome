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

#========== GENUS LEVEL =======
# Most abundant phyla
genus_sp <- tax_glom(physeq_filtered, taxrank = rank_names(physeq_filtered)[6], NArm = F) # agglomérate to genus level
genus_sp <- transform_sample_counts(genus_sp, function(x) x / sum(x)) # relative abundance
TopASV_g <- names(sort(taxa_sums(genus_sp), TRUE)[1:15]) # most abundant genus
top15_water_g <- prune_species(TopASV_g, genus_sp)
top15_water_g <- prune_taxa(taxa_sums(top15_water_g) > 0, top15_water_g)
top_genus <- as.data.frame(tax_table(top15_water_g))

# seperate data by lakes
physeq_merge_by <- merge_samples(genus_sp, "lake_origin")
sp_RA <- transform_sample_counts(physeq_merge_by, function(x) x / sum(x))
sp_melt <- psmelt(sp_RA) # create a datframe from phyloseq object
sp_melt$Genus <- as.character(sp_melt$Genus) # convert into character

# Keep 10 most abundant phylum
sumtot_sp <- sp_melt %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarize(sum = sum(Abundance)) %>%
  filter(Genus %in% top_genus$Genus) %>%
  filter(!(Genus %in% c("", " p__uncultured", NA, "Unknown Genus")))

sp_melt$Genus[!(sp_melt$Genus %in% sumtot_sp$Genus)] <- "Other"

sp_melt$totalAbundance <- sum(sp_melt$Abundance)

data_sp_mod <- sp_melt %>%
  group_by(Genus, Sample) %>%
  summarise(Abundance = sum(Abundance)) %>%
  dplyr::distinct()

# Specify order of samples
data_sp_mod$Sample <- factor(data_sp_mod$Sample, levels = c("CERL_PER", "CERL_LEP", "CERS_PER", 
                                                            "CERS_LEP", "CRE_PER", "CRE_LEP", "VSS_PER", 
                                                            "VSS_LEP", "LGP_PER", "LGP_LEP", "CSM_PER", "CSM_LEP"))

color_palette <- c("#fdbf6f",
                   "#a276c7",
                   "#E74C3C",
                   "#ba8cdc",
                   "#8491B4",
                   "#707ca2",
                   "#f7b75d",
                   "#364B9A",
                   "lightgrey",
                   "#eba63a",
                   "#a085ae",
                   "#c57300"
)

labels <- c(
  expression(italic("Aeromonas")),
  expression(italic("Candidatus Bacilloplasma")),
  expression(italic("Cetobacterium")),
  expression(italic("Clostridium sensu stricto 1")),
  expression(italic("Cyanobium PCC-6307")),
  expression(italic("LD29")),
  expression(italic("Legionella")),
  expression(italic("Mycobacterium")),
  "Other",
  expression(italic("Plesiomonas")),
  expression(italic("Romboutsia")),
  expression(italic("Silvianigrella"))
)

# create barplot
barplot_for_sp_genus <- ggplot2::ggplot(data = data_sp_mod, ggplot2::aes(x = Sample, y = Abundance, fill = Genus)) +
  ggplot2::geom_bar(stat = "identity", position = "stack") +
  ggplot2::scale_fill_manual(values = color_palette,
                             labels = labels) +
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

path_to_my_object = here::here("figures", "figSM2.png")
ggplot2::ggsave(filename = path_to_my_object, plot = barplot_for_sp_genus, device = "png")
