#======== PROJECT COM2LIFE ========
## Estimates alpha diversity based on hill numbers 
## First load the libraries and physeq_object. Calculates alpha dveristy on phylogenetic (qo, q1) taxonomic (q0, q1)
## create fonction to plot the inidces.

#create a file 
dir.create("outputs/hill", recursive = TRUE)

#libraries
library(hillR)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(dplyr)
library(cowplot)

# Load my physeq object 
physeq_filtered <- readRDS(here::here("data",
                             "mon_objet_physeq_filtered.rds"))

# Estimates alpha diversity with Hill numbers
#'@taxonomic_diversity
comm = physeq_filtered@otu_table

#q = 0 (species richness)
tax_q0_alpha <- hillR::hill_taxa(comm,q = 0)

#q = 1 (shannon entropy)
tax_q1_alpha <- hillR::hill_taxa(comm, q = 1)

#'@Phylogenetic_diversity
tree= physeq_filtered@phy_tree

#q = 0 (species richness)
phylo_q0_alpha <- hillR::hill_phylo(comm,tree,q = 0)

#q = 1 (shannon entropy)
phylo_q1_alpha <- hillR::hill_phylo(comm,tree,q = 1)

# Have a physeq_al^ha for alpha plot
physeq_alpha <- physeq_filtered
## put every new df in physeq@sam_data
physeq_alpha@sam_data$taxo_q0 <- tax_q0_alpha
physeq_alpha@sam_data$taxo_q1 <- tax_q1_alpha
physeq_alpha@sam_data$phylo_q0 <- phylo_q0_alpha
physeq_alpha@sam_data$phylo_q1<- phylo_q1_alpha

# Save my phyloseq object 
path_to_my_object = here::here("Data","mon_objet_physeq_alpha.rds")
saveRDS(physeq_alpha, file = path_to_my_object)


#create metadata_alpha object
metadata_alpha <- data.frame(sample_data(physeq_alpha))
# Définir l'ordre souhaité des niveaux
order <- c("CERL", "CERS", "CRE", "VSS", "LGP", "CSM")
# Réordonner la variable lake en fonction de cet ordre
metadata_alpha$lake <- factor(metadata_alpha$lake, levels = order)
# Save my metadata object 
path_to_my_object = here::here("Data","metadata_alpha.rds")
saveRDS(metadata_alpha, file = path_to_my_object)
