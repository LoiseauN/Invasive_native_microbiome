#======== PROJECT COM2LIFE ========
## Estimates beta diversity based on hill numbers 
## First load the libraries and physeq_object. Calculates beta dveristy on phylogenetic (qo, q1) taxonomic (q0, q1)

#Load libraries
library(dplyr)
library(tidyverse) 
library(hillR)

# Load my physeq filtered object 
physeq_filtered <- readRDS(here::here("Data",
                                      "mon_objet_physeq_filtered.rds"))

metadata <- physeq_filtered@sam_data %>%
  `class<-`(NULL) %>%               
  `attr<-`("package", NULL) %>%     
  as.data.frame() 

# Save
path_to_my_object = here::here("Data","metadata.rds")
saveRDS(metadata, file = path_to_my_object)

#'@taxonomic_diversity
comm = physeq_filtered@otu_table

# Beta diversity
#q = 0 (species richness)
tax_q0_beta <- hillR::hill_taxa_parti_pairwise(comm,q = 0) #output='matrix' for distance matric in output
tax_q0_beta$beta_diss <- 1 - tax_q0_beta$region_similarity

#q = 1 (shannon entropy)
tax_q1_beta <- hillR::hill_taxa_parti_pairwise(comm,q = 1) #output='matrix' for distance matric in output
tax_q1_beta$beta_diss <- 1 - tax_q1_beta$region_similarity


#'@Phylogenetic_diversity
tree= physeq_filtered@phy_tree

# Beta diversity
#q = 0 (species richness)
phylo_q0_beta <- hillR::hill_phylo_parti_pairwise(comm,tree,q = 0)#output='matrix' for distance matric in output
phylo_q0_beta$beta_diss <- 1 - phylo_q0_beta$region_similarity

#q = 1 (shannon entropy)
phylo_q1_beta <- hillR::hill_phylo_parti_pairwise(comm,tree,q = 1) #output='matrix' for distance matric in output
phylo_q1_beta$beta_diss <- 1 - phylo_q1_beta$region_similarity

## keep only beta_diss
phylo_q0_beta <- phylo_q0_beta %>%
  select(-q, -PD_gamma, -PD_alpha, -PD_beta, -local_similarity, -region_similarity) %>%
  rename(sample_b = site1, sample_a = site2, phylo_q0 = beta_diss)
#save 
path_to_my_object = here::here("Data","phylo_q0_beta.rds")
saveRDS(phylo_q0_beta, file = path_to_my_object)

phylo_q1_beta <- phylo_q1_beta %>%
  select(-q, -PD_gamma, -PD_alpha, -PD_beta, -local_similarity, -region_similarity) %>%
  rename(sample_b = site1, sample_a = site2, phylo_q1 = beta_diss)
#save
path_to_my_object = here::here("Data","phylo_q1_beta.rds")
saveRDS(phylo_q1_beta, file = path_to_my_object)
        
tax_q0_beta <- tax_q0_beta %>%
  select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
  rename(sample_b = site1, sample_a = site2, taxo_q0 = beta_diss)
#save
path_to_my_object = here::here("Data","tax_q0_beta.rds")
saveRDS(tax_q0_beta, file = path_to_my_object)

tax_q1_beta <- tax_q1_beta %>%
  select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
  rename(sample_b = site1, sample_a = site2, taxo_q1 = beta_diss)
#save
path_to_my_object = here::here("Data","tax_q1_beta.rds")
saveRDS(tax_q1_beta, file = path_to_my_object)

#Create dataframe for analysis of dissimilarity
beta_diss <- phylo_q0_beta %>%
  merge(phylo_q1_beta, by = c('sample_a', 'sample_b')) %>%
  merge(tax_q0_beta, by = c('sample_a', 'sample_b')) %>%
  merge(tax_q1_beta, by = c('sample_a', 'sample_b'))

beta_diss$spl_comb <- sapply(1:nrow(beta_diss), function(x) { 
  tmp <- sort(c(as.character(beta_diss$sample_a)[x], 
                as.character(beta_diss$sample_b)[x]))
  paste0(tmp[1], "__", tmp[2])
})
beta_diss$spl_comb <- gsub("\\.", "-", beta_diss$spl_comb)
beta_diss$sp_comb <- sapply(beta_diss$spl_comb, function(x) {
  sapply(str_split_fixed(x, "__", 2), function(xx) {
    metadata[metadata$samples == xx, "origin"]
  }) %>% sort() %>% paste0(collapse = "__")
})
beta_diss$reg_comb <- sapply(beta_diss$spl_comb, function(x) {
  sapply(str_split_fixed(x, "__", 2), function(xx) {
    metadata[metadata$samples == xx, "lake"]
  }) %>% sort() %>% paste0(collapse = "__")
})

# reformat and clean the df2 ----
beta_diss <- beta_diss %>% 
  separate(sp_comb, c("species_a", "species_b"), sep = "__", remove = FALSE) %>%
  separate(reg_comb, c("region_a", "region_b"), sep = "__", remove = FALSE) %>%
  dplyr::select(sample_a, sample_b, species_a, species_b, region_a, region_b, 
                spl_comb, sp_comb, reg_comb,
                everything())

beta_diss$sp_lev <- rep("intra_species", nrow(beta_diss))
beta_diss$sp_lev[beta_diss$species_a != beta_diss$species_b] <- "inter_species"
beta_diss$reg_lev <- rep("intra_region", nrow(beta_diss))
beta_diss$reg_lev[beta_diss$region_a != beta_diss$region_b] <- "inter_region"

# change factor levels
for (i in names(beta_diss)[!(names(beta_diss) %in% c("taxo_q0", "taxo_q1", "phylo_q0", "phylo_q1"))]) {
  beta_diss[, i] <- factor(beta_diss[,i], levels = sort(unique(beta_diss[,i])))
}

path_to_my_object = here::here("Data","beta_diss.rds")
saveRDS(beta_diss, file = path_to_my_object)
