#======== PROJECT COM2LIFE ========
## permanova on beta diversity matrix to see if there is lake effect or species effect

# =========== PERCA FLUVIATILIS ============
#create a physeq_per for only perca fluviatilis
physeq_per <- subset_samples(physeq_filtered, !(grepl("GAR", sample_names(physeq_filtered))))
# Save my phyloseq object 
path_to_my_object = here::here("Data","mon_objet_physeq_per.rds")
saveRDS(physeq_per, file = path_to_my_object)

metadata_per <- metadata[grep("PER", metadata$origin), ]

#'@taxonomic_diversity
comm = physeq_per@otu_table

# Beta diversity
#q = 0 (species richness)
tax_q0_beta <- hillR::hill_taxa_parti_pairwise(comm,q = 0) 
tax_q0_beta$beta_diss <- 1 - tax_q0_beta$region_similarity

#q = 1 (shannon entropy)
tax_q1_beta <- hillR::hill_taxa_parti_pairwise(comm,q = 1) 
tax_q1_beta$beta_diss <- 1 - tax_q1_beta$region_similarity


#'@Phylogenetic_diversity
tree= physeq_per@phy_tree

# Beta diversity
#q = 0 (species richness)
phylo_q0_beta <- hillR::hill_phylo_parti_pairwise(comm,tree,q = 0)
phylo_q0_beta$beta_diss <- 1 - phylo_q0_beta$region_similarity

#q = 1 (shannon entropy)
phylo_q1_beta <- hillR::hill_phylo_parti_pairwise(comm,tree,q = 1) 
phylo_q1_beta$beta_diss <- 1 - phylo_q1_beta$region_similarity

#keep only diversity indices
phylo_q0_beta <- phylo_q0_beta %>%
  dplyr::select(-q, -PD_gamma, -PD_alpha, -PD_beta, -local_similarity, -region_similarity) %>%
  dplyr::rename(sample_b = site1, sample_a = site2, phylo_q0 = beta_diss)
phylo_q1_beta <- phylo_q1_beta %>%
  dplyr::select(-q, -PD_gamma, -PD_alpha, -PD_beta, -local_similarity, -region_similarity) %>%
  dplyr::rename(sample_b = site1, sample_a = site2, phylo_q1 = beta_diss)

tax_q0_beta <- tax_q0_beta %>%
  dplyr::select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
  dplyr::rename(sample_b = site1, sample_a = site2, taxo_q0 = beta_diss)
tax_q1_beta <- tax_q1_beta %>%
  dplyr::select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
  dplyr::rename(sample_b = site1, sample_a = site2, taxo_q1 = beta_diss)

# transform the dissimilarity indices in matrix
taxo_q0_beta <- as.data.frame(tax_q0_beta)
taxo_q1_beta <- as.data.frame(tax_q1_beta) 
phylo_q0_beta <- as.data.frame(phylo_q0_beta) 
phylo_q1_beta <- as.data.frame(phylo_q1_beta)

samples <- unique(c(taxo_q0_beta$sample_a, taxo_q0_beta$sample_b))
tax_q0_beta_matrix <- matrix(0, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))
for (i in 1:nrow(taxo_q0_beta)) {
  sample1 <- taxo_q0_beta[i, "sample_a"]
  sample2 <- taxo_q0_beta[i, "sample_b"]
  distance <- taxo_q0_beta[i, "taxo_q0"]
  tax_q0_beta_matrix[sample1, sample2] <- distance
  tax_q0_beta_matrix[sample2, sample1] <- distance
}

samples <- unique(c(taxo_q1_beta$sample_a, taxo_q1_beta$sample_b))
tax_q1_beta_matrix <- matrix(0, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))
for (i in 1:nrow(taxo_q1_beta)) {
  sample1 <- taxo_q1_beta[i, "sample_a"]
  sample2 <- taxo_q1_beta[i, "sample_b"]
  distance <- taxo_q1_beta[i, "taxo_q1"]
  tax_q1_beta_matrix[sample1, sample2] <- distance
  tax_q1_beta_matrix[sample2, sample1] <- distance
}

samples <- unique(c(phylo_q0_beta$sample_a, phylo_q0_beta$sample_b))
phylo_q0_beta_matrix <- matrix(0, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))
for (i in 1:nrow(phylo_q0_beta)) {
  sample1 <- phylo_q0_beta[i, "sample_a"]
  sample2 <- phylo_q0_beta[i, "sample_b"]
  distance <- phylo_q0_beta[i, "phylo_q0"]
  phylo_q0_beta_matrix[sample1, sample2] <- distance
  phylo_q0_beta_matrix[sample2, sample1] <- distance
}

samples <- unique(c(phylo_q1_beta$sample_a, phylo_q1_beta$sample_b))
phylo_q1_beta_matrix <- matrix(0, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))
for (i in 1:nrow(phylo_q1_beta)) {
  sample1 <- phylo_q1_beta[i, "sample_a"]
  sample2 <- phylo_q1_beta[i, "sample_b"]
  distance <- phylo_q1_beta[i, "phylo_q1"]
  phylo_q1_beta_matrix[sample1, sample2] <- distance
  phylo_q1_beta_matrix[sample2, sample1] <- distance
}

# Save my phyloseq object 
path_to_my_object = here::here("Data","matrix_phylo_q0_per.rds")
saveRDS(phylo_q0_beta_matrix, file = path_to_my_object)
path_to_my_object = here::here("Data","matrix_phylo_q1_per.rds")
saveRDS(phylo_q1_beta_matrix, file = path_to_my_object)
path_to_my_object = here::here("Data","matrix_taxo_q0_per.rds")
saveRDS(tax_q0_beta_matrix, file = path_to_my_object)
path_to_my_object = here::here("Data","matrix_taxo_q1_per.rds")
saveRDS(tax_q1_beta_matrix, file = path_to_my_object)

# ========== LEPOMIS GIBBOSUS ===============
#create a physeq_per for only perca fluviatilis
physeq_lep <- subset_samples(physeq_filtered, !(grepl("PER", sample_names(physeq_filtered))))
# Save my phyloseq object 
path_to_my_object = here::here("Data","mon_objet_physeq_lep.rds")
saveRDS(physeq_lep, file = path_to_my_object)

metadata_lep <- metadata[grep("GAR", rownames(metadata)), ]

#'@taxonomic_diversity
comm = physeq_lep@otu_table

# Beta diversity
#q = 0 (species richness)
tax_q0_beta <- hillR::hill_taxa_parti_pairwise(comm,q = 0) 
tax_q0_beta$beta_diss <- 1 - tax_q0_beta$region_similarity

#q = 1 (shannon entropy)
tax_q1_beta <- hillR::hill_taxa_parti_pairwise(comm,q = 1) 
tax_q1_beta$beta_diss <- 1 - tax_q1_beta$region_similarity


#'@Phylogenetic_diversity
tree= physeq_lep@phy_tree

# Beta diversity
#q = 0 (species richness)
phylo_q0_beta <- hillR::hill_phylo_parti_pairwise(comm,tree,q = 0)
phylo_q0_beta$beta_diss <- 1 - phylo_q0_beta$region_similarity

#q = 1 (shannon entropy)
phylo_q1_beta <- hillR::hill_phylo_parti_pairwise(comm,tree,q = 1) 
phylo_q1_beta$beta_diss <- 1 - phylo_q1_beta$region_similarity

## supprimer les autres colonnes sur la beta garde que beta_diss
phylo_q0_beta <- phylo_q0_beta %>%
  select(-q, -PD_gamma, -PD_alpha, -PD_beta, -local_similarity, -region_similarity) %>%
  rename(sample_b = site1, sample_a = site2, phylo_q0 = beta_diss)
phylo_q1_beta <- phylo_q1_beta %>%
  select(-q, -PD_gamma, -PD_alpha, -PD_beta, -local_similarity, -region_similarity) %>%
  rename(sample_b = site1, sample_a = site2, phylo_q1 = beta_diss)

tax_q0_beta <- tax_q0_beta %>%
  select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
  rename(sample_b = site1, sample_a = site2, taxo_q0 = beta_diss)
tax_q1_beta <- tax_q1_beta %>%
  select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
  rename(sample_b = site1, sample_a = site2, taxo_q1 = beta_diss)

# transform the dissimilarity indice in matrix
taxo_q0_beta <- as.data.frame(tax_q0_beta)
taxo_q1_beta <- as.data.frame(tax_q1_beta) 
phylo_q0_beta <- as.data.frame(phylo_q0_beta) 
phylo_q1_beta <- as.data.frame(phylo_q1_beta)

samples <- unique(c(taxo_q0_beta$sample_a, taxo_q0_beta$sample_b))
tax_q0_beta_matrix <- matrix(0, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))
for (i in 1:nrow(taxo_q0_beta)) {
  sample1 <- taxo_q0_beta[i, "sample_a"]
  sample2 <- taxo_q0_beta[i, "sample_b"]
  distance <- taxo_q0_beta[i, "taxo_q0"]
  tax_q0_beta_matrix[sample1, sample2] <- distance
  tax_q0_beta_matrix[sample2, sample1] <- distance
}

samples <- unique(c(taxo_q1_beta$sample_a, taxo_q1_beta$sample_b))
tax_q1_beta_matrix <- matrix(0, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))
for (i in 1:nrow(taxo_q1_beta)) {
  sample1 <- taxo_q1_beta[i, "sample_a"]
  sample2 <- taxo_q1_beta[i, "sample_b"]
  distance <- taxo_q1_beta[i, "taxo_q1"]
  tax_q1_beta_matrix[sample1, sample2] <- distance
  tax_q1_beta_matrix[sample2, sample1] <- distance
}

samples <- unique(c(phylo_q0_beta$sample_a, phylo_q0_beta$sample_b))
phylo_q0_beta_matrix <- matrix(0, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))
for (i in 1:nrow(phylo_q0_beta)) {
  sample1 <- phylo_q0_beta[i, "sample_a"]
  sample2 <- phylo_q0_beta[i, "sample_b"]
  distance <- phylo_q0_beta[i, "phylo_q0"]
  phylo_q0_beta_matrix[sample1, sample2] <- distance
  phylo_q0_beta_matrix[sample2, sample1] <- distance
}

samples <- unique(c(phylo_q1_beta$sample_a, phylo_q1_beta$sample_b))
phylo_q1_beta_matrix <- matrix(0, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))
for (i in 1:nrow(phylo_q1_beta)) {
  sample1 <- phylo_q1_beta[i, "sample_a"]
  sample2 <- phylo_q1_beta[i, "sample_b"]
  distance <- phylo_q1_beta[i, "phylo_q1"]
  phylo_q1_beta_matrix[sample1, sample2] <- distance
  phylo_q1_beta_matrix[sample2, sample1] <- distance
}

#Save
path_to_my_object = here::here("Data","matrix_phylo_q0_lep.rds")
saveRDS(phylo_q0_beta_matrix, file = path_to_my_object)
path_to_my_object = here::here("Data","matrix_phylo_q1_lep.rds")
saveRDS(phylo_q1_beta_matrix, file = path_to_my_object)
path_to_my_object = here::here("Data","matrix_taxo_q0_lep.rds")
saveRDS(tax_q0_beta_matrix, file = path_to_my_object)
path_to_my_object = here::here("Data","matrix_taxo_q1_lep.rds")
saveRDS(tax_q1_beta_matrix, file = path_to_my_object)

# load data 
metadata <- readRDS(here::here("data",
                               "metadata.rds"))

# ===== PERMANOVA - LAKE EFFECT - PERCA FLUVIATILIS =====
metadata_per <- metadata[grep("PER", metadata$origin), ]


permanova_tax_q0_per <- vegan::adonis2(tax_q0_matrix_per ~ lake,
                                       data = metadata_per,
                                       perm = 1000)
permanova_tax_q0_per

permanova_tax_q1_per <- vegan::adonis2(tax_q1_matrix_per ~ lake,
                                       data = metadata_per,
                                       perm = 1000)
permanova_tax_q1_per

permanova_phylo_q0_per <- vegan::adonis2(phylo_q0_matrix_per ~ lake,
                                         data = metadata_per,
                                         perm = 1000)
permanova_phylo_q0_per

permanova_phylo_q1_per <- vegan::adonis2(phylo_q1_matrix_per ~ lake,
                                         data = metadata_per,
                                         perm = 1000)
permanova_phylo_q1_per

## Betadisp()
# === TAXO Q0 ===
tax_q0_dist_per <- dist(tax_q0_matrix_per)

disp <- betadisper(tax_q0_dist_per, metadata_per$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

# === TAXO Q1 ===
tax_q1_dist_per <- dist(tax_q1_matrix_per)

disp <- betadisper(tax_q1_dist_per, metadata_per$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

# === PHYLO Q0 ===
phylo_q0_dist_per <- dist(phylo_q0_matrix_per)

disp <- betadisper(phylo_q0_dist_per, metadata_per$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

# === PHYLO Q1 ===
phylo_q1_dist_per <- dist(phylo_q1_matrix_per)

disp <- betadisper(phylo_q1_dist_per, metadata_per$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

# ===== PERMANOVA - LAKE EFFECT - LEPOMIS GIBBOSUS =====
metadata_lep <- metadata[grep("LEP", metadata$origin), ]


lepmanova_tax_q0_lep <- vegan::adonis2(tax_q0_matrix_lep ~ lake,
                                       data = metadata_lep,
                                       lepm = 1000)
lepmanova_tax_q0_lep

lepmanova_tax_q1_lep <- vegan::adonis2(tax_q1_matrix_lep ~ lake,
                                       data = metadata_lep,
                                       lepm = 1000)
lepmanova_tax_q1_lep

lepmanova_phylo_q0_lep <- vegan::adonis2(phylo_q0_matrix_lep ~ lake,
                                         data = metadata_lep,
                                         lepm = 1000)
lepmanova_phylo_q0_lep

lepmanova_phylo_q1_lep <- vegan::adonis2(phylo_q1_matrix_lep ~ lake,
                                         data = metadata_lep,
                                         lepm = 1000)
lepmanova_phylo_q1_lep

## Betadisp()
# === TAXO Q0 ===
tax_q0_dist_lep <- dist(tax_q0_matrix_lep)

disp <- betadisper(tax_q0_dist_lep, metadata_lep$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

# === TAXO Q1 ===
tax_q1_dist_lep <- dist(tax_q1_matrix_lep)

disp <- betadisper(tax_q1_dist_lep, metadata_lep$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

# === PHYLO Q0 ===
phylo_q0_dist_lep <- dist(phylo_q0_matrix_lep)

disp <- betadisper(phylo_q0_dist_lep, metadata_lep$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

# === PHYLO Q1 ===
phylo_q1_dist_lep <- dist(phylo_q1_matrix_lep)

disp <- betadisper(phylo_q1_dist_lep, metadata_lep$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

## No parwise adonis here, as lake variance homogeneity shows no significant differences