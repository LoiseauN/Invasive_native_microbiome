#======== PROJECT COM2LIFE ========
## Stackbarplots for water, help of L.Ezzat 

# Load physeq object
physeq <- readRDS(here::here("data", "mon_objet_physeq_13_02.rds"))

# Filter samples to keep only samples from water
physeq_lake <- subset_samples(physeq, !(grepl("BOI|JAB|VAI|SED|PER|GAR", sample_names(physeq))))
rows_to_remove <- grepl("BOI|JAB|VAI|SED|PER|GAR", rownames(sample_data(physeq_lake)))

# keep the sam data
physeq_lake@sam_data <- sample_data(physeq_lake)[!rows_to_remove, ]

# Modify lakes's name
# Access samples name
sample_names_current <- sample_names(physeq_lake)
# Modify names
sample_names_modified <- sample_names_current %>%
  gsub("CHA", "CSM", .) %>%
  gsub("CRJ", "CERL", .) %>% 
  gsub("TRI", "CERS", .) %>%
  gsub("CTL", "CRE", .) %>%
  gsub("GDP", "LGP", .) %>%
  gsub("VER", "VSS", .)
# Apply new names
sample_names(physeq_lake) <- sample_names_modified
# Modify names for sam_data
sample_data <- sample_data(physeq_lake)
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
sample_data(physeq_lake) <- sample_data

path_to_my_object = here::here("data", "physeq_lake.rds")
saveRDS(physeq_lake, file = path_to_my_object)

#========== PHYLUM LEVEL =======
#Most abundant phyla
phyla_lake <- tax_glom(physeq_lake, taxrank = rank_names(physeq_lake)[2], NArm = F) 
phyla_lake <- transform_sample_counts(phyla_lake, function(x) x / sum(x)) 
TopASV_p <- names(sort(taxa_sums(phyla_lake), TRUE)[1:6]) 
top6_water_p <- prune_species(TopASV_p, phyla_lake)
top6_water_p <- prune_taxa(taxa_sums(top6_water_p) > 0, top6_water_p)
top_phyla <- as.data.frame(tax_table(top6_water_p))

# Separate data and caluclate abundace
physeq_merge_by <- merge_samples(phyla_lake, "lake")
lake_RA <- transform_sample_counts(physeq_merge_by, function(x) x / sum(x))
lake_melt <- psmelt(lake_RA) 
lake_melt$Phylum <- as.character(lake_melt$Phylum)

# Keep 6 most abundant phyla
sumtot_lake <- lake_melt %>%
  dplyr::group_by(Phylum) %>%
  dplyr::summarize(sum = sum(Abundance)) %>%
  dplyr::filter(Phylum %in% top_phyla$Phylum) %>%
  dplyr::filter(!(Phylum %in% c("", " p__uncultured", NA, "Unknown Phylum")))

lake_melt$Phylum[!(lake_melt$Phylum %in% sumtot_lake$Phylum)] <- "Other"

lake_melt$totalAbundance <- sum(lake_melt$Abundance)

data_lake_mod <- lake_melt %>%
  dplyr::group_by(Phylum, Sample) %>%
  dplyr::summarise(Abundance = sum(Abundance)) %>%
  dplyr::distinct()

# Specify order of lakes
data_sp_mod$Sample <- factor(data_sp_mod$Sample, levels = c("CERL", "CERS", "CRE", "VSS","LGP", "CSM"))

color_palette <- c("#439EB7", "#E28B55", "#DCB64D", "#4CA198", "#835B82", "#645135", "#FF966B")

# Create barplot
barplot_for_all_lakes <- ggplot(data = data_lake_mod, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

barplot_for_all_lakes


path_to_my_object = here::here("Figures", "barplot_for_all_lakes.png")
ggsave(filename = path_to_my_object, plot = barplot_for_all_lakes, device = "png")

#========== GENUS LEVEL =======
# Most abundant phyla
genus_sp <- tax_glom(physeq_lake, taxrank = rank_names(physeq_lake)[6], NArm = F) # agglomérate to genus level
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
data_sp_mod$Sample <- factor(data_sp_mod$Sample, levels = c("CERL_WATER", "CERS_WATER", "CRE_WATER", 
                                                            "VSS_WATER", "LGP_WATER", "CSM_WATER"))

color_palette <- c("#A50021", "#D82632", "#F72836", "#FF7857","#FFAc75", "#FFD699", 
                   "#FFF2BD","#EBFFFF", "#BDF9FF", "#99EBFF", "#75D3FF", "#8CB2FF","#5991FF","#57B0ff")


# create barplot
barplot_lake_genus <- ggplot(data = data_sp_mod, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

barplot_lake_genus

path_to_my_object = here::here("Figures", "barplot_lake_genus.png")
ggplot2::ggsave(filename = path_to_my_object, plot = barplot_lake_genus, device = "png")
