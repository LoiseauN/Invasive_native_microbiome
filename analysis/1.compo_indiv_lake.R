#======== PROJECT COM2LIFE ========
## Stackbarplots for species


# Load physeq object
physeq <- readRDS(here::here("data",
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
path_to_my_object = here::here("data","mon_objet_physeq_filtered.rds")
saveRDS(physeq_filtered, file = path_to_my_object)

metadata <- sample.data.frame(physeq_filtered)

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

color_palette <- c("#A50021", "#D82632", "#F76D5E", "#FFAD72","#FFE099", "#FFFFBF", 
                   "#E0FFFF","#AAF7FF", "#72D8FF", "#3FA0FF","#264CFF")


# create barplot
barplot_for_sp <- ggplot(data = data_sp_mod, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

barplot_for_sp

path_to_my_object = here::here("figures", "barplot_for_sp.png")
ggplot2::ggsave(filename = path_to_my_object, plot = barplot_for_sp, device = "png")

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

color_palette <- c("#A50021", "#D82632", "#F72836", "#FF7857","#FFAc75", "#FFD699", 
                   "#FFF2BD","#EBFFFF", "#BDF9FF", "#99EBFF", "#75D3FF", "#57B0ff")


# create barplot
barplot_for_sp_genus <- ggplot(data = data_sp_mod, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

barplot_for_sp_genus

path_to_my_object = here::here("figures", "barplot_for_sp_genus.png")
ggplot2::ggsave(filename = path_to_my_object, plot = barplot_for_sp_genus, device = "png")
