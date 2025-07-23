#======== PROJECT COM2LIFE ========
## Estimates beta diversity based on hill numbers for metabolite data  
## First load the libraries and physeq_object. Calculates beta diversity on phylogenetic (qo, q1) taxonomic (q0, q1)

## =============== PCOA and DbRDA ================== ##
# Load objects 
physeq_filtered <- readRDS(here::here("Data",
                             "mon_objet_physeq_filtered.rds"))
phylo_q0_beta <- readRDS(here::here("Data",
                                      "phylo_q0_beta.rds"))
phylo_q1_beta <- readRDS(here::here("Data",
                                    "phylo_q1_beta.rds"))
taxo_q0_beta <- readRDS(here::here("Data",
                                    "tax_q0_beta.rds"))
taxo_q1_beta <- readRDS(here::here("Data",
                                    "tax_q1_beta.rds"))

# library
library(hillR)
library(tidyverse)
library(ape)
library(gridExtra)
library(ggpubr)
library(ggplot2)

# transform the dissimilarity indice in matrix
taxo_q0_beta <- as.data.frame(taxo_q0_beta)
taxo_q1_beta <- as.data.frame(taxo_q1_beta) 
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

## =========== Make the PCoA ==================== ##
lake_colors <- c(
  "CSM" = "#608F3D",
  "CERL" = "#41AEBD",
  "CERS" = "#97E9D5",
  "CRE" = "#F4DE3A",
  "LGP" = "#A2CF49",
  "VSS" = "#FCB11C"
)

plot_pcoa <- function(matrix_dist, title) {
  sample_names <- rownames(matrix_dist)
  lake_names <- substr(sample_names, 1, 4)
  lake_names[lake_names == "CRE."] <- "CRE"
  lake_names[lake_names == "VSS."] <- "VSS"
  lake_names[lake_names == "LGP."] <- "LGP"
  lake_names[lake_names == "CSM."] <- "CSM"
  
  # Perform PCoA
  pcoa_result <- pcoa(matrix_dist)
  
  # Assign lake_names to pcoa_result$lake_names
  pcoa_result$lake_names <- lake_names
  
  # Extract PCoA coordinates
  pcoa_coordinates <- pcoa_result$vectors
  
  # Define point shapes based on sample names
  point_shapes <- ifelse(grepl("PER", rownames(matrix_dist)), 16, 1)
  
  # Create the plot using ggplot2
  plot <- ggplot(data = pcoa_coordinates, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(color = lake_names, shape = factor(point_shapes)), size = 3) +
    xlab("Dim 1") +
    ylab("Dim 2") +
    labs(title = title) +
    scale_color_manual(values = lake_colors, name = "Lakes") +  
    scale_shape_manual(values = c(16, 1),
                       labels = c("Perca fluviatilis", "Lepomis gibbosus")) +  
    guides(shape = guide_legend(title = "Species",     
                                override.aes = list(shape = c(16, 1)))) +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"), 
      panel.background = element_blank()
    )
  
  return(plot)
}

# Generate PCoA plots for taxonomic q0 and q1 dissimilarities
plot_taxo_q0_pcoa <- plot_pcoa(tax_q0_beta_matrix, "Taxonomic dissimilarity (q0)")
plot_taxo_q1_pcoa <- plot_pcoa(tax_q1_beta_matrix, "Taxonomic dissimilarity (q1)")

# Generate PCoA plots for phylogenetic q0 and q1 dissimilarities
plot_phylo_q0_pcoa <- plot_pcoa(phylo_q0_beta_matrix, "Phylogenetic dissimilarity (q0)")
plot_phylo_q1_pcoa <- plot_pcoa(phylo_q1_beta_matrix, "Phylogenetic dissimilarity (q1)")

# Combine plots
pcoa_diss <- ggarrange(plot_taxo_q0_pcoa, plot_taxo_q1_pcoa, plot_phylo_q0_pcoa, plot_phylo_q1_pcoa,
                       labels = c("A", "B", "C", "D"),
                       ncol = 2, nrow = 3, 
                       common.legend = TRUE, legend = "right")

path_to_my_object = here::here("Figures","hill", "pcoa_diss_all.png")
ggsave(filename = path_to_my_object, plot = pcoa_diss, device = "png")


## =========== Make the dbRDA ==================== ##
library(vegan)
library(ade4)
library(tibble)
library(ggplot2)

#create metadata
metadata <- physeq_filtered@sam_data %>%
  `class<-`(NULL) %>%               
  `attr<-`("package", NULL) %>%     
  as.data.frame() %>%
  dplyr::select(-PROK_1_median, -PROK_2_median, -LARGE_CELLS_median, -VLP_TOT_median, -CELLS_TOT_median, -Secchi) %>%
  dplyr::rename(Temperature = Temperature_median,
                Chla = Chla_median, 
                Salinity = Salinity_median,
                O2 = Oxygen_median,
                pH = pH_median)

metadata <- column_to_rownames(metadata, var = "samples")

custom_order <- c("CSM", "CERL","CERS", "CRE","LGP", "VSS")
metadata$lake <- factor(metadata$lake, levels = custom_order)

#parameters for dbRDA
colors <-  c(
  "CSM" = "#608F3D",
  "CERL" = "#41AEBD",
  "CERS" = "#97E9D5",
  "CRE" = "#F4DE3A",
  "LGP" = "#A2CF49",
  "VSS" = "#FCB11C"
)
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


# ========== Perform dbRDA for PERCA FLUVIATILIS ===============
dbrda_plot <- function(matrix_dist, title) {
  dbrda_result <- dbrda(matrix_dist ~ Temperature + Salinity + O2 + Chla + pH, data = metadata_per)
  
  site_scores <- as.data.frame(scores(dbrda_result, display = "sites"))
  selected_vars <- metadata_per[, c("Temperature", "Salinity", "O2", "pH", "Chla")]
  env_fit <- envfit(dbrda_result, env = selected_vars)
  arrow_data <- as.data.frame(scores(env_fit, display = "vectors"))
  
  arrow_data <- arrow_data %>%
    dplyr::mutate(
      arrow_end_x = dbRDA1,   
      arrow_end_y = dbRDA2   
    )
  
  dbRDA_plot <- ggplot(data = site_scores, aes(x = dbRDA1, y = dbRDA2, color = metadata_per$lake)) +
    geom_point(size = 3) +
    geom_segment(data = arrow_data, 
                 aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
                 arrow = arrow(length = unit(0.5, "cm")),
                 color = "red") +
    scale_color_manual(values = colors, 
                       labels = names(colors),
                       name = "Lakes") +
    xlab("dbRDA1") + 
    ylab("dbRDA2") + 
    labs(title = title) +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"), 
      panel.background = element_blank(),
      axis.ticks = element_line(colour = "black"),
      legend.position = "right") +
    annotate(geom = "segment", x = 0, y = 0, xend = max(site_scores$dbRDA1), yend = 0,
             linetype = "dashed", color = "black") +
    annotate(geom = "segment", x = 0, y = 0, xend = 0, yend = max(site_scores$dbRDA2),
             linetype = "dashed", color = "black") +
    geom_text(data = arrow_data,
              aes(x = arrow_end_x, y = arrow_end_y, label = rownames(arrow_data)),
              color = "red", size = 4,
              hjust = 0.8, vjust = 0.7) 
  
  return(dbRDA_plot)
}

# Perform DB-RDA and plot for taxonomic q0 and q1 dissimilarities
dbRDA_taxo_q0 <- dbrda_plot(tax_q0_beta_matrix, "Taxonomic (q0)")
dbRDA_taxo_q1 <- dbrda_plot(tax_q1_beta_matrix, "Taxonomic (q1)")

# Perform DB-RDA and plot for phylogenetic q0 and q1 dissimilarities
dbRDA_phylo_q0 <- dbrda_plot(phylo_q0_beta_matrix, "Phylogenetic (q0)")
dbRDA_phylo_q1 <- dbrda_plot(phylo_q1_beta_matrix, "Phylogenetic (q1)")

# Combine plots
dbRDA_per <- ggarrange(dbRDA_taxo_q0, dbRDA_taxo_q1, dbRDA_phylo_q0, dbRDA_phylo_q1,
                       labels = c("A", "B", "C", "D"),
                       ncol = 2, nrow = 2, 
                       common.legend = TRUE, legend = "right")

path_to_my_object = here::here("Figures","hill", "dbRDA_per.png")
ggsave(filename = path_to_my_object, plot = dbRDA_per, device = "png")


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

# ========== Perform dbRDA for LEPOMIS GIBBOSUS ===============
dbrda_plot <- function(matrix_dist, title) {
  dbrda_result <- dbrda(matrix_dist ~ Temperature + Salinity + O2 + Chla + pH, data = metadata_lep)
  
  site_scores <- as.data.frame(scores(dbrda_result, display = "sites"))
  selected_vars <- metadata_lep[, c("Temperature", "Salinity", "O2", "pH", "Chla")]
  env_fit <- envfit(dbrda_result, env = selected_vars)
  arrow_data <- as.data.frame(scores(env_fit, display = "vectors"))
  
  arrow_data <- arrow_data %>%
    dplyr::mutate(
      arrow_end_x = dbRDA1,   
      arrow_end_y = dbRDA2   
    )
  
  dbRDA_plot <- ggplot(data = site_scores, aes(x = dbRDA1, y = dbRDA2, color = metadata_lep$lake)) +
    geom_point(size = 3) +
    geom_segment(data = arrow_data, 
                 aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
                 arrow = arrow(length = unit(0.5, "cm")),
                 color = "red") +
    scale_color_manual(values = colors, 
                       labels = names(colors),
                       name = "Lakes") +
    xlab("dbRDA1") + 
    ylab("dbRDA2") + 
    labs(title = title) +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"), 
      panel.background = element_blank(),
      axis.ticks = element_line(colour = "black"),
      legend.position = "right") +
    annotate(geom = "segment", x = 0, y = 0, xend = max(site_scores$dbRDA1), yend = 0,
             linetype = "dashed", color = "black") +
    annotate(geom = "segment", x = 0, y = 0, xend = 0, yend = max(site_scores$dbRDA2),
             linetype = "dashed", color = "black") +
    geom_text(data = arrow_data,
              aes(x = arrow_end_x, y = arrow_end_y, label = rownames(arrow_data)),
              color = "red", size = 4,
              hjust = 0.8, vjust = 0.7) 
  
  return(dbRDA_plot)
}

# Perform DB-RDA and plot for taxonomic q0 and q1 dissimilarities
dbRDA_taxo_q0 <- dbrda_plot(tax_q0_beta_matrix, "Taxonomic (q0)")
dbRDA_taxo_q1 <- dbrda_plot(tax_q1_beta_matrix, "Taxonomic (q1)")

# Perform DB-RDA and plot for phylogenetic q0 and q1 dissimilarities
dbRDA_phylo_q0 <- dbrda_plot(phylo_q0_beta_matrix, "Phylogenetic (q0)")
dbRDA_phylo_q1 <- dbrda_plot(phylo_q1_beta_matrix, "Phylogenetic (q1)")

# Combine plots
dbRDA_lep <- ggarrange(dbRDA_taxo_q0, dbRDA_taxo_q1, dbRDA_phylo_q0, dbRDA_phylo_q1,
                       labels = c("A", "B", "C", "D"),
                       ncol = 2, nrow = 2, 
                       common.legend = TRUE, legend = "right")

path_to_my_object = here::here("Figures","hill", "dbRDA_lep.png")
ggsave(filename = path_to_my_object, plot = dbRDA_lep, device = "png")
