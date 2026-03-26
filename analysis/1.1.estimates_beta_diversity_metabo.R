#======== PROJECT COM2LIFE ========
## Estimates beta diversity based on hill numbers for metabolite data  
## First load the libraries and physeq_object. Calculates beta dveristy on phylogenetic (qo, q1) taxonomic (q0, q1)


# ======= ALL METABO =======
# load data
all_metabo <- readRDS(here::here("data",
                                 "all_metabo_filtered.rds"))
# change with the right lakes name
all_metabo <- all_metabo %>%
  mutate(new_rowname = rownames(all_metabo) %>%
           gsub("CHA", "CSM", .) %>%
           gsub("GDP", "LGP", .) %>%
           gsub("VER", "VSS", .) %>%
           gsub("CRJ", "CERL", .) %>%
           gsub("TRI", "CERS", .) %>%
           gsub("CTL", "CRE", .))

rownames(all_metabo) <- all_metabo$new_rowname
all_metabo <- all_metabo[, -which(names(all_metabo) == "new_rowname")]

all_metabo <- all_metabo %>%
  mutate(Lake = gsub("CHA", "CSM", Lake) %>%
           gsub("GDP", "LGP", .) %>%
           gsub("VER", "VSS", .) %>%
           gsub("CRJ", "CERL", .) %>%
           gsub("TRI", "CERS", .) %>%
           gsub("CTL", "CRE", .))

# Save
path_to_my_object = here::here("data","all_metabo_filtered.rds")
saveRDS(all_metabo, file = path_to_my_object)

# Delete character columns 
all_metabo <- all_metabo %>%
  select(-Species, -Lake)

# Convert all remaining columns to numeric
all_metabo <- all_metabo %>%
  mutate(across(everything(), as.numeric))

#'@taxonomic_diversity
comm = all_metabo

# Beta diversity
#q = 0 (species richness)
tax_q0_beta <- hillR::hill_taxa_parti_pairwise(comm,q = 0) #output='matrix' for distance matric in output
tax_q0_beta$beta_diss <- 1 - tax_q0_beta$region_similarity

#q = 1 (shannon entropy)
tax_q1_beta <- hillR::hill_taxa_parti_pairwise(comm,q = 1) #output='matrix' for distance matric in output
tax_q1_beta$beta_diss <- 1 - tax_q1_beta$region_similarity

## keep only beta_diss
tax_q0_beta <- tax_q0_beta %>%
  dplyr::select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
  dplyr::rename(sample_b = site1, sample_a = site2, taxo_q0 = beta_diss)
#save
path_to_my_object = here::here("data","tax_q0_beta_metabo.rds")
saveRDS(tax_q0_beta, file = path_to_my_object)

tax_q1_beta <- tax_q1_beta %>%
  dplyr::select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
  dplyr::rename(sample_b = site1, sample_a = site2, taxo_q1 = beta_diss)
#save
path_to_my_object = here::here("data","tax_q1_beta_metabo.rds")
saveRDS(tax_q1_beta, file = path_to_my_object)

# transform the dissimilarity indice in matrix
taxo_q0_beta <- as.data.frame(tax_q0_beta)
taxo_q1_beta <- as.data.frame(tax_q1_beta) 

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

# ======= ANNOTATED METABO =======
# load data
annotated_metabo <- readRDS(here::here("Data",
                                       "annotated_metabo_filtered.rds"))
# change with the right lakes name
annotated_metabo <- annotated_metabo %>%
  mutate(new_rowname = rownames(annotated_metabo) %>%
           gsub("CHA", "CSM", .) %>%
           gsub("GDP", "LGP", .) %>%
           gsub("VER", "VSS", .) %>%
           gsub("CRJ", "CERL", .) %>%
           gsub("TRI", "CERS", .) %>%
           gsub("CTL", "CRE", .))

rownames(annotated_metabo) <- annotated_metabo$new_rowname
annotated_metabo <- annotated_metabo[, -which(names(annotated_metabo) == "new_rowname")]

# Save
path_to_my_object = here::here("data","annotated_metabo_filtered.rds")
saveRDS(annotated_metabo, file = path_to_my_object)

# Convert all remaining columns to numeric
annotated_metabo <- annotated_metabo %>%
  mutate(across(everything(), as.numeric))

#'@taxonomic_diversity
comm = annotated_metabo

# Beta diversity
#q = 0 (species richness)
tax_q0_beta <- hillR::hill_taxa_parti_pairwise(comm,q = 0) #output='matrix' for distance matric in output
tax_q0_beta$beta_diss <- 1 - tax_q0_beta$region_similarity

#q = 1 (shannon entropy)
tax_q1_beta <- hillR::hill_taxa_parti_pairwise(comm,q = 1) #output='matrix' for distance matric in output
tax_q1_beta$beta_diss <- 1 - tax_q1_beta$region_similarity

## keep only beta_diss
tax_q0_beta <- tax_q0_beta %>%
  dplyr::select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
  dplyr::rename(sample_b = site1, sample_a = site2, taxo_q0 = beta_diss)
#save
path_to_my_object = here::here("data","tax_q0_beta_annotated_metabo.rds")
saveRDS(tax_q0_beta, file = path_to_my_object)

tax_q1_beta <- tax_q1_beta %>%
  dplyr::select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
  dplyr::rename(sample_b = site1, sample_a = site2, taxo_q1 = beta_diss)
#save
path_to_my_object = here::here("data","tax_q1_beta_annotated_metabo.rds")
saveRDS(tax_q1_beta, file = path_to_my_object)

# transform the dissimilarity indice in matrix
taxo_q0_beta <- as.data.frame(tax_q0_beta_annotated_metabo)
taxo_q1_beta <- as.data.frame(tax_q1_beta_annotated_metabo) 

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
