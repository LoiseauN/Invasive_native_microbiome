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
#saveRDS(physeq_filtered, file = path_to_my_object)

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


data_sp_mod$lake <- sub("\\_.*", "", data_sp_mod$Sample )
data_sp_mod$species <- sub('.*_', '', data_sp_mod$Sample)

data_sp_mod$species <-  gsub("PER", "Perca fluviatilis", data_sp_mod$species)
data_sp_mod$species <-  gsub("LEP", "Lepomis gibbosus", data_sp_mod$species)

custom_order <- c("CSM", "LGP", "VSS", "CRE", "CERS", "CERL")

species_order <- c("Perca fluviatilis", "Lepomis gibbosus")

data_sp_mod$species <- factor(data_sp_mod$species, levels = species_order)


barplot_for_sp <- ggplot2::ggplot(data = data_sp_mod,
                                  ggplot2::aes(x = lake, y = Abundance, fill = Phylum)) +
  ggplot2::geom_bar(stat = "identity", position = "stack") +
  ggplot2::scale_x_discrete(limits = custom_order) +
  ggplot2::scale_fill_manual(values = color_palette) +
  ggplot2::theme_bw() +
  ggplot2::labs(y = "Relative Abundance") +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 12),
    axis.title.y = ggplot2::element_text(size = 12),
    axis.text.x  = ggplot2::element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y  = ggplot2::element_text(size = 10),
    legend.title = ggplot2::element_text(size = 14),
    legend.text  = ggplot2::element_text(size = 12),
    panel.border     = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.line        = ggplot2::element_line(colour = "black"),
    strip.text       = ggplot2::element_text(face = "italic", size = 11)  # italic species names
  ) +
  ggplot2::facet_wrap(ggplot2::vars(species))



barplot_for_sp

path_to_my_object = here::here("figures", "figure1.png")
ggplot2::ggsave(filename = path_to_my_object, plot = barplot_for_sp, device = "png")





## IN SUP


#======== PROJECT COM2LIFE ========
## Genus-level stacked barplots: fish (by species) and water, side by side, ONE shared legend

library(phyloseq)
library(tidyverse)

physeq <- readRDS(here::here("Data", "mon_objet_physeq_17_06.rds"))

# ---- 1. Subsets ------------------------------------------------------------
physeq_fish  <- prune_samples(!grepl("BOI|JAB|VAI|SED|ADN", sample_names(physeq)), physeq)
physeq_water <- prune_samples(!grepl("BOI|JAB|VAI|SED|PER|LEP|GAR", sample_names(physeq)), physeq)
physeq_fish  <- prune_taxa(taxa_sums(physeq_fish)  > 0, physeq_fish)
physeq_water <- prune_taxa(taxa_sums(physeq_water) > 0, physeq_water)

# ---- 2. Lake renaming (same rules as before, written once) -----------------
rename_sample <- function(x) {
  x |> gsub("CHA", "CSM", x = _) |> gsub("CRJ", "CERL", x = _) |>
    gsub("TRI", "CERS", x = _) |> gsub("CTL", "CRE", x = _) |>
    gsub("GDP", "LGP", x = _) |> gsub("VER", "VSS", x = _)
}
rename_lake <- function(x) {
  x |> gsub("CHA", "CSM", x = _) |> gsub("CRJ1", "CERL", x = _) |>
    gsub("CRJ2", "CERS", x = _) |> gsub("CTL", "CRE", x = _) |>
    gsub("GDP", "LGP", x = _) |> gsub("VER", "VSS", x = _)
}
prep_physeq <- function(ps) {
  sample_names(ps) <- rename_sample(sample_names(ps))
  sd <- data.frame(sample_data(ps), check.names = FALSE)
  for (col in c("lake", "lake_origin"))
    if (col %in% colnames(sd)) sd[[col]] <- rename_lake(as.character(sd[[col]]))
  if ("samples" %in% colnames(sd)) sd$samples <- rename_sample(as.character(sd$samples))
  sd$lake_grp <- sub("_.*", "", as.character(sd[["lake"]]))
  sample_data(ps) <- sample_data(sd)
  ps
}
physeq_fish  <- prep_physeq(physeq_fish)
physeq_water <- prep_physeq(physeq_water)

# ---- 3. Genus table (relative abundance per group) -------------------------
unclassified <- c("", "uncultured", "Unknown", "unknown", "Unknown Genus",
                  "unclassified", "Incertae_Sedis", "Incertae Sedis")

genus_table <- function(ps, group_var) {
  rank_genus <- if ("Genus" %in% rank_names(ps)) "Genus" else rank_names(ps)[6]
  g <- tax_glom(ps, taxrank = rank_genus, NArm = FALSE)
  g <- transform_sample_counts(g, function(x) x / sum(x))
  g <- merge_samples(g, group_var)
  g <- transform_sample_counts(g, function(x) x / sum(x))
  psmelt(g) |>
    mutate(Genus = trimws(sub("^\\s*g__", "", as.character(.data[[rank_genus]]))),
           Genus = ifelse(is.na(Genus) | Genus %in% unclassified |
                            grepl("^uncultured|^unclassified|^Unknown", Genus), "Other", Genus),
           group = as.character(Sample)) |>
    group_by(Genus, group) |>
    summarise(Abundance = sum(Abundance), .groups = "drop")
}

panel_levels <- c("Perca fluviatilis", "Lepomis gibbosus", "Water")

fish_tab <- genus_table(physeq_fish, "lake_origin") |>
  mutate(lake  = sub("_.*", "", group),
         panel = recode(sub(".*_", "", group),
                        PER = "Perca fluviatilis", LEP = "Lepomis gibbosus")) |>
  filter(panel %in% panel_levels)          # keeps only PER / LEP groups

water_tab <- genus_table(physeq_water, "lake_grp") |>
  mutate(lake = group, panel = "Water")

all_tab <- bind_rows(fish_tab, water_tab)

# ---- 4. ONE set of top genera for both datasets -> same legend & colors ----
top_genera <- all_tab |>
  filter(Genus != "Other") |>
  group_by(Genus) |> summarise(tot = sum(Abundance)) |>
  slice_max(tot, n = 10, with_ties = FALSE) |> pull(Genus)

all_tab <- all_tab |>
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other")) |>
  group_by(Genus, lake, panel) |>
  summarise(Abundance = sum(Abundance), .groups = "drop") |>
  mutate(Genus = factor(Genus, levels = c(top_genera, "Other")),
         panel = factor(panel, levels = panel_levels))

# ---- 5. Colors: "Other" is the only grey ------------------------------------
color_palette <- c("#3C5488", "#4DBBD5", "#B0E0E6", "#8491B4", "#cab2d6",
                   "#AD002A", "#F39B7F", "#fdbf6f", "#B09C85", "#7E6148")
fill_values <- c(setNames(rep_len(color_palette, length(top_genera)), top_genera),
                 Other = "lightgrey")

# ---- 6. Plot: 3 panels in one row, single legend ---------------------------
custom_order <- c("CSM", "LGP", "VSS", "CRE", "CERS", "CERL")

barplot_genus <- ggplot(all_tab, aes(x = lake, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(vars(panel), nrow = 1) +
  scale_x_discrete(limits = custom_order) +
  scale_fill_manual(values = fill_values, drop = FALSE) +
  labs(x = NULL, y = "Relative Abundance") +
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 12),
    axis.text.x  = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 10),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = "right",
    strip.text   = element_text(size = 12, face = "italic"),
    strip.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

barplot_genus

ggsave(here::here("figures", "figure_genus_fish_water.png"), plot = barplot_genus,
       device = "png", width = 14, height = 5, dpi = 300)