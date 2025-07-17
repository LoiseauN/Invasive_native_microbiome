#======== PROJECT COM2LIFE ========
## Estimates beta diversity based on hill numbers for metabolite data  
## First load the libraries and physeq_object. Calculates beta dveristy on phylogenetic (qo, q1) taxonomic (q0, q1)

# load libraries 
library(dplyr)
library(spaa)
library(hillR)
library(ape)
library(ggpubr)

# ======= ALL METABO =======
# load data
all_metabo <- readRDS(here::here("Data",
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
path_to_my_object = here::here("Data","all_metabo_filtered.rds")
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
  select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
  rename(sample_b = site1, sample_a = site2, taxo_q0 = beta_diss)
#save
path_to_my_object = here::here("Data","tax_q0_beta_metabo.rds")
saveRDS(tax_q0_beta, file = path_to_my_object)

tax_q1_beta <- tax_q1_beta %>%
  select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
  rename(sample_b = site1, sample_a = site2, taxo_q1 = beta_diss)
#save
path_to_my_object = here::here("Data","tax_q1_beta_metabo.rds")
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
  lake_names[lake_names == "CRE_"] <- "CRE"
  lake_names[lake_names == "VSS_"] <- "VSS"
  lake_names[lake_names == "LGP_"] <- "LGP"
  lake_names[lake_names == "CSM_"] <- "CSM"
  
  pcoa_result <- pcoa(matrix_dist)
  pcoa_result$lake_names <- lake_names
  
  pcoa_coordinates <- pcoa_result$vectors
  
  point_shapes <- ifelse(grepl("_P", rownames(matrix_dist)), 16, 1)
  
  plot <- ggplot(data = pcoa_coordinates, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(color = lake_names, shape = factor(point_shapes)), size = 3) +
    xlab("Dim 1") +
    ylab("Dim 2") +
    labs(title = title) +
    scale_color_manual(values = lake_colors, name = "Lakes") +  
    scale_shape_manual(values = c(16, 1),
                       labels = c(expression(italic("Perca fluviatilis")), expression(italic("Lepomis gibbosus")))) +  
    guides(shape = guide_legend(title = "Species",     
                                override.aes = list(shape = c(16, 1)))) +
    theme(
      legend.text = element_text(size = 10.5),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"), 
      panel.background = element_blank()
    )
  
  return(plot)
}

# Generate PCoA plots for taxonomic q0 and q1 dissimilarities
plot_taxo_q0_metabo <- plot_pcoa(tax_q0_beta_matrix, "Metabolite dissimilarity (q0)")
plot_taxo_q1_metabo <- plot_pcoa(tax_q1_beta_matrix, "Metabolite dissimilarity (q1)")

# Combine plots
pcoa_diss_all_metabo <- ggarrange(plot_taxo_q0_metabo, plot_taxo_q1_metabo,
                       labels = c("A", "B"),
                       ncol = 2, nrow = 1, 
                       common.legend = TRUE, legend = "right")

path_to_my_object = here::here("Figures","hill", "pcoa_diss_all_metabo.png")
ggsave(filename = path_to_my_object, plot = pcoa_diss_all_metabo, device = "png")


# ===== PCoA FOR PERCA AND LEPOMIS ON DIFFERENT PLOT =====
library(hillR)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(vegan)

#create object for Perca  
metabo_per <- all_metabo[grepl("_P", rownames(all_metabo)), ]
# Save my phyloseq object 
path_to_my_object = here::here("Data","metabo_per.rds")
saveRDS(metabo_per, file = path_to_my_object)
#cretae objetc for Lepomis
metabo_lep <- all_metabo[grepl("_G", rownames(all_metabo)), ]
# Save my phyloseq object 
path_to_my_object = here::here("Data","metabo_lep.rds")
saveRDS(metabo_lep, file = path_to_my_object)

# Function to generate dissimilarity matrices and plot PCoA
generate_and_plot <- function(comm, species_name, prefix) {
  # Dissimilarity matrices
  tax_q0_beta <- hill_taxa_parti_pairwise(comm, q = 0)
  tax_q0_beta$beta_diss <- 1 - tax_q0_beta$region_similarity
  
  tax_q1_beta <- hill_taxa_parti_pairwise(comm, q = 1)
  tax_q1_beta$beta_diss <- 1 - tax_q1_beta$region_similarity
  
  # Keep only beta_diss
  tax_q0_beta <- tax_q0_beta %>%
    select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
    rename(sample_b = site1, sample_a = site2, taxo_q0 = beta_diss)
  
  tax_q1_beta <- tax_q1_beta %>%
    select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
    rename(sample_b = site1, sample_a = site2, taxo_q1 = beta_diss)
  
  # Transform dissimilarity indices into matrix
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
  
  # Generate PCoA plots
  plot_q0 <- plot_pcoa(tax_q0_beta_matrix, paste(species_name, " - Metabolite dissimilarity (q0)"))
  plot_q1 <- plot_pcoa(tax_q1_beta_matrix, paste(species_name, " - Metabolite dissimilarity (q1)"))
  
  # Return plots
  return(list(plot_q0, plot_q1))
}

# Function to plot PCoA
plot_pcoa <- function(matrix_dist, title) {
  sample_names <- rownames(matrix_dist)
  lake_names <- substr(sample_names, 1, 4)
  lake_names[lake_names == "CRE_"] <- "CRE"
  lake_names[lake_names == "VSS_"] <- "VSS"
  lake_names[lake_names == "LGP_"] <- "LGP"
  lake_names[lake_names == "CSM_"] <- "CSM"
  
  pcoa_result <- pcoa(matrix_dist)
  pcoa_result$lake_names <- lake_names
  
  pcoa_coordinates <- pcoa_result$vectors
  
  plot <- ggplot(data = pcoa_coordinates, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(color = lake_names), size = 3) +
    xlab("Dim 1") +
    ylab("Dim 2") +
    labs(title = title) +
    scale_color_manual(values = lake_colors, name = "Lakes") +  
    theme(
      legend.text = element_text(size = 10.5),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"), 
      panel.background = element_blank()
    )
  
  return(plot)
}

# Apply functions
species_data <- list(
  list(species_name = "Perca fluviatilis", prefix = "per", comm = metabo_per),
  list(species_name = "Lepomis gibbosus", prefix = "lep", comm = metabo_lep)
)

plots <- list()
for (species in species_data) {
  species_name <- species$species_name
  prefix <- species$prefix
  comm <- species$comm
  
  plots[[length(plots) + 1]] <- generate_and_plot(comm, species_name, prefix)
}

# Unpack plots
plot_taxo_q0_metabo_per <- plots[[1]][[1]]
plot_taxo_q1_metabo_per <- plots[[1]][[2]]
plot_taxo_q0_metabo_lep <- plots[[2]][[1]]
plot_taxo_q1_metabo_lep <- plots[[2]][[2]]

# Combine plots
pcoa_diss_all_metabo_sp <- ggarrange(plot_taxo_q0_metabo_per, plot_taxo_q1_metabo_per,
                                     plot_taxo_q0_metabo_lep, plot_taxo_q1_metabo_lep,
                                     labels = c("A", "B", "C", "D"),
                                     ncol = 2, nrow = 2,
                                     common.legend = TRUE, legend = "right")

# Save combined plot
path_to_my_object <- here::here("Figures", "hill", "pcoa_diss_all_metabo_sp_separate.png")
ggsave(filename = path_to_my_object, plot = pcoa_diss_all_metabo_sp, device = "png")


# ======= ANNOTATED METABO =======
# load data
annotated_metabo <- readRDS(here::here("Data",
                                 "metabo_annotated_filtered.rds"))
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
path_to_my_object = here::here("Data","annotated_metabo_filtered.rds")
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
  select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
  rename(sample_b = site1, sample_a = site2, taxo_q0 = beta_diss)
#save
path_to_my_object = here::here("Data","tax_q0_beta_annotated_metabo.rds")
saveRDS(tax_q0_beta, file = path_to_my_object)

tax_q1_beta <- tax_q1_beta %>%
  select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
  rename(sample_b = site1, sample_a = site2, taxo_q1 = beta_diss)
#save
path_to_my_object = here::here("Data","tax_q1_beta_annotated_metabo.rds")
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
  lake_names[lake_names == "CRE_"] <- "CRE"
  lake_names[lake_names == "VSS_"] <- "VSS"
  lake_names[lake_names == "LGP_"] <- "LGP"
  lake_names[lake_names == "CSM_"] <- "CSM"
  
  pcoa_result <- pcoa(matrix_dist)
  pcoa_result$lake_names <- lake_names
  
  pcoa_coordinates <- pcoa_result$vectors
  
  point_shapes <- ifelse(grepl("_P", rownames(matrix_dist)), 16, 1)
  
  plot <- ggplot(data = pcoa_coordinates, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(color = lake_names, shape = factor(point_shapes)), size = 3) +
    xlab("Dim 1") +
    ylab("Dim 2") +
    labs(title = title) +
    scale_color_manual(values = lake_colors, name = "Lakes") +  
    scale_shape_manual(values = c(16, 1),
                       labels = c(expression(italic("Perca fluviatilis")), expression(italic("Lepomis gibbosus")))) +  
    guides(shape = guide_legend(title = "Species",     
                                override.aes = list(shape = c(16, 1)))) +
    theme(
      legend.text = element_text(size = 10.5),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"), 
      panel.background = element_blank()
    )
  
  return(plot)
}

# Generate PCoA plots for taxonomic q0 and q1 dissimilarities
plot_taxo_q0_metabo <- plot_pcoa(tax_q0_beta_matrix, "Metabolite dissimilarity (q0)")
plot_taxo_q1_metabo <- plot_pcoa(tax_q1_beta_matrix, "Metabolite dissimilarity (q1)")

# Combine plots
pcoa_diss_annotated_metabo <- ggarrange(plot_taxo_q0_metabo, plot_taxo_q1_metabo,
                                  labels = c("A", "B"),
                                  ncol = 2, nrow = 1, 
                                  common.legend = TRUE, legend = "right")

path_to_my_object = here::here("Figures","hill", "pcoa_diss_annotated_metabo.png")
ggsave(filename = path_to_my_object, plot = pcoa_diss_annotated_metabo, device = "png")


# ===== PCoA FOR PERCA AND LEPOMIS ON DIFFERENT PLOT FOR ANNOTATED METABO =====
library(hillR)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(vegan)

#create object for Perca  
metabo_per <- annotated_metabo[grepl("_P", rownames(annotated_metabo)), ]
# Save my phyloseq object 
path_to_my_object = here::here("Data","annotated_metabo_per.rds")
saveRDS(metabo_per, file = path_to_my_object)
#cretae objetc for Lepomis
metabo_lep <- annotated_metabo[grepl("_G", rownames(annotated_metabo)), ]
# Save my phyloseq object 
path_to_my_object = here::here("Data","annotated_metabo_lep.rds")
saveRDS(metabo_lep, file = path_to_my_object)

# Function to generate dissimilarity matrices and plot PCoA
generate_and_plot <- function(comm, species_name, prefix) {
  # Dissimilarity matrices
  tax_q0_beta <- hill_taxa_parti_pairwise(comm, q = 0)
  tax_q0_beta$beta_diss <- 1 - tax_q0_beta$region_similarity
  
  tax_q1_beta <- hill_taxa_parti_pairwise(comm, q = 1)
  tax_q1_beta$beta_diss <- 1 - tax_q1_beta$region_similarity
  
  # Keep only beta_diss
  tax_q0_beta <- tax_q0_beta %>%
    select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
    rename(sample_b = site1, sample_a = site2, taxo_q0 = beta_diss)
  
  tax_q1_beta <- tax_q1_beta %>%
    select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
    rename(sample_b = site1, sample_a = site2, taxo_q1 = beta_diss)
  
  # Transform dissimilarity indices into matrix
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
  
  # Generate PCoA plots
  plot_q0 <- plot_pcoa(tax_q0_beta_matrix, paste(species_name, " - Annotated metabolite dissimilarity (q0)"))
  plot_q1 <- plot_pcoa(tax_q1_beta_matrix, paste(species_name, " - Annotated metabolite dissimilarity (q1)"))
  
  # Return plots
  return(list(plot_q0, plot_q1))
}

# Function to plot PCoA
plot_pcoa <- function(matrix_dist, title) {
  sample_names <- rownames(matrix_dist)
  lake_names <- substr(sample_names, 1, 4)
  lake_names[lake_names == "CRE_"] <- "CRE"
  lake_names[lake_names == "VSS_"] <- "VSS"
  lake_names[lake_names == "LGP_"] <- "LGP"
  lake_names[lake_names == "CSM_"] <- "CSM"
  
  pcoa_result <- pcoa(matrix_dist)
  pcoa_result$lake_names <- lake_names
  
  pcoa_coordinates <- pcoa_result$vectors
  
  plot <- ggplot(data = pcoa_coordinates, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(color = lake_names), size = 3) +
    xlab("Dim 1") +
    ylab("Dim 2") +
    labs(title = title) +
    scale_color_manual(values = lake_colors, name = "Lakes") +  
    theme(
      legend.text = element_text(size = 10.5),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"), 
      panel.background = element_blank()
    )
  
  return(plot)
}

# Apply functions
species_data <- list(
  list(species_name = "Perca fluviatilis", prefix = "per", comm = metabo_per),
  list(species_name = "Lepomis gibbosus", prefix = "lep", comm = metabo_lep)
)

plots <- list()
for (species in species_data) {
  species_name <- species$species_name
  prefix <- species$prefix
  comm <- species$comm
  
  plots[[length(plots) + 1]] <- generate_and_plot(comm, species_name, prefix)
}

# Unpack plots
plot_taxo_q0_annotated_metabo_per <- plots[[1]][[1]]
plot_taxo_q1_annotated_metabo_per <- plots[[1]][[2]]
plot_taxo_q0_annotated_metabo_lep <- plots[[2]][[1]]
plot_taxo_q1_annotated_metabo_lep <- plots[[2]][[2]]

# Combine plots
pcoa_diss_annotated_metabo_sp <- ggarrange(plot_taxo_q0_annotated_metabo_per, plot_taxo_q1_annotated_metabo_per,
                                     plot_taxo_q0_annotated_metabo_lep, plot_taxo_q1_annotated_metabo_lep,
                                     labels = c("A", "B", "C", "D"),
                                     ncol = 2, nrow = 2,
                                     common.legend = TRUE, legend = "right")

# Save combined plot
path_to_my_object <- here::here("Figures", "hill", "pcoa_diss_anotated_metabo_sp_separate.png")
ggsave(filename = path_to_my_object, plot = pcoa_diss_annotated_metabo_sp, device = "png")
