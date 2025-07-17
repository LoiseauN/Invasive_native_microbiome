# Load physeq object
physeq_water_fish <- readRDS(here::here("Data",
                             "mon_objet_physeq_filtered_all.rds"))

# Filter physeq to keep only samples of fish
physeq_water_fish <- subset_samples(physeq, !(grepl("BOI|JAB|VAI|SED", sample_names(physeq))))
rows_to_remove <- grepl("BOI|JAB|VAI|SED", rownames(sample_data(physeq_water_fish)))

# keep the sam data
physeq_water_fish@sam_data <- sample_data(physeq_water_fish)[!rows_to_remove, ]

path_to_my_object = here::here("Data","physeq_water_fish.rds")
saveRDS(physeq_water_fish, file = path_to_my_object)


library(VennDiagram)
# ====== VENN DIAGRAM ======
physeq_water_fish <- transform_sample_counts(physeq_water_fish, function(x) x / sum(x))
otu_table <- as.data.frame(physeq_water_fish@otu_table)
tax_table <- as.data.frame(physeq_water_fish@tax_table)

# Filtrer les échantillons par groupe
lep_samples <- rownames(otu_table)[grepl("GAR", rownames(otu_table))]
per_samples <- rownames(otu_table)[grepl("PER", rownames(otu_table))]
water_samples <- rownames(otu_table)[grepl("ADN", rownames(otu_table))]

# Extraire les ASVs pour chaque groupe
asvs_lep <- colnames(otu_table)[colSums(otu_table[lep_samples, ] > 0) > 0]
asvs_per <- colnames(otu_table)[colSums(otu_table[per_samples, ] > 0) > 0]
asvs_water <- colnames(otu_table)[colSums(otu_table[water_samples, ] > 0) > 0]

# Calculer les intersections manuellement pour vérifier
inter_lep_per <- intersect(asvs_lep, asvs_per)
inter_per_water <- intersect(asvs_per, asvs_water)
inter_lep_water <- intersect(asvs_lep, asvs_water)
inter_all <- intersect(inter_lep_per, asvs_water)

# Afficher les intersections pour vérifier
cat("Intersections entre lep et per:", inter_lep_per, "\n")
cat("Intersections entre per et water:", inter_per_water, "\n")
cat("Intersections entre lep et water:", inter_lep_water, "\n")
cat("Intersections entre les trois:", inter_all, "\n")

# Vérifiez les longueurs des intersections
cat("Nombre d'ASVs partagés entre lep et per:", length(inter_lep_per), "\n")
cat("Nombre d'ASVs partagés entre per et water:", length(inter_per_water), "\n")
cat("Nombre d'ASVs partagés entre lep et water:", length(inter_lep_water), "\n")
cat("Nombre d'ASVs partagés entre les trois groupes:", length(inter_all), "\n")
# Nombre d'ASVs partagés uniquement entre lep et water, mais pas avec per
n13_only = length(inter_lep_water) - length(inter_all)

# Nombre d'ASVs partagés uniquement entre lep et per, mais pas avec water
n12_only = length(inter_lep_per) - length(inter_all)

# Nombre d'ASVs partagés uniquement entre per et water, mais pas avec lep
n23_only = length(inter_per_water) - length(inter_all)

# Nombre total d'ASVs dans chaque groupe
area1 = length(asvs_lep)
area2 = length(asvs_per)
area3 = length(asvs_water)

# Dessiner le diagramme de Venn
venn_plot <- draw.triple.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  n12 = n12_only + length(inter_all),  # ASVs partagés entre lep et per (inclut ceux partagés avec water)
  n23 = n23_only + length(inter_all),  # ASVs partagés entre per et water (inclut ceux partagés avec lep)
  n13 = n13_only + length(inter_all),  # ASVs partagés entre lep et water (inclut ceux partagés avec per)
  n123 = length(inter_all),  # ASVs partagés entre les trois groupes
  category = c("lep", "per", "water"),
  fill = c("red", "blue", "green"),
  euler.d = TRUE
)

 

# Trouver les ASVs partagés entre GAR et PER
asvs_shared_fish <- intersect(asvs_lep, asvs_per)

# Afficher les ASVs partagés
asvs_shared_fish

# Récupérer les informations taxonomiques des ASVs partagés
taxo_shared <- tax_table[asvs_shared_fish, ]

# Extraire les informations de genre
genus_shared <- taxo_shared[, "Genus"]

# Afficher les genres des ASVs partagés
genus_shared
# Convertir les informations de genre en data frame
genus_shared <- data.frame(
  ASV = asvs_shared_fish,
  Genus = as.vector(genus_shared)
)

# Filtrer pour supprimer les genres NA
genus_shared <- genus_shared %>%
  filter(!is.na(Genus))

# Récupérer les abondances des ASVs partagés dans les échantillons GAR et PER
otu_shared_lep <- otu_table[lep_samples, asvs_shared_fish]
otu_shared_per <- otu_table[per_samples, asvs_shared_fish]

# Convertir en dataframes
df_lep <- as.data.frame(otu_shared_lep)
df_per <- as.data.frame(otu_shared_per)

# Ajouter les noms d'échantillons comme colonne
df_lep$Sample <- rownames(df_lep)
df_per$Sample <- rownames(df_per)

# Fusionner les dataframes GAR et PER
df_combined <- bind_rows(
  df_lep %>% mutate(Group = "GAR"),
  df_per %>% mutate(Group = "PER")
)

# Transformer les données en format long pour ggplot2
df_long <- df_combined %>%
  tidyr::pivot_longer(cols = -c(Sample, Group), names_to = "ASV", values_to = "Abundance")

# Ajouter l'information de genre
df_long <- df_long %>%
  left_join(genus_shared, by = "ASV")

# Agréger les abondances par genre et groupe
df_genus_abundance <- df_long %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE)) %>%
  ungroup()
# Filtrer pour supprimer les genres NA
df_genus_abundance <- df_genus_abundance %>%
  filter(!is.na(Genus))

df_genus_abundance <- df_genus_abundance %>%
  mutate_at(vars(-Genus), ~ . / sum(., na.rm = TRUE))

# Filtrer pour ne conserver que les genres avec une abondance > 1%
df_genus_abundance <- df_genus_abundance %>%
  filter(TotalAbundance >= 0.01)

# Créer le graphique des abondances par genre
ggplot(df_genus_abundance, aes(x = Genus, y = TotalAbundance)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Abondances des genres partagés entre GAR et PER",
       x = "Genus",
       y = "Abondance totale") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



