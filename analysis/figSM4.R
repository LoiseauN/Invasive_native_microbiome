#======== PROJECT COM2LIFE ======

# Make PCoA for annotated metabo. Species separate.

# ===== PCoA FOR PERCA AND LEPOMIS ON DIFFERENT PLOT FOR ANNOTATED METABO =====
library(hillR)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(vegan)
library(ape)
library(dplyr)
library(here)

# Load the full annotated metabolite table (adjust the file name to your source object)
# annotated_metabo <- readRDS(here::here("Data", "annotated_metabo.rds"))

# ===== LAKE COLORS =====
lake_colors <- c(
  CSM  = "#00332AFF",
  LGP  = "#035236FF",
  VSS  = "#227548FF",
  CRE  = "#60A360FF",
  CERS = "#C5D163FF",
  CERL = "#F2F26DFF"
)

# Create object for Perca
metabo_per <- annotated_metabo[grepl("_P", rownames(annotated_metabo)), ]
path_to_my_object <- here::here("Data", "annotated_metabo_per.rds")
saveRDS(metabo_per, file = path_to_my_object)

# Create object for Lepomis
metabo_lep <- annotated_metabo[grepl("_G", rownames(annotated_metabo)), ]
path_to_my_object <- here::here("Data", "annotated_metabo_lep.rds")
saveRDS(metabo_lep, file = path_to_my_object)

# Function to plot PCoA (defined BEFORE generate_and_plot, which calls it)
plot_pcoa <- function(matrix_dist, title) {
  sample_names <- rownames(matrix_dist)
  lake_names <- substr(sample_names, 1, 4)
  lake_names <- sub("_$", "", lake_names)   # CRE_ -> CRE, VSS_ -> VSS, etc.
  
  pcoa_result <- ape::pcoa(matrix_dist)
  
  # Data frame with coordinates + lake names
  pcoa_coordinates <- as.data.frame(pcoa_result$vectors)
  pcoa_coordinates$lake_names <- factor(lake_names, levels = names(lake_colors))
  
  plot <- ggplot2::ggplot(
    data = pcoa_coordinates,
    ggplot2::aes(x = Axis.1, y = Axis.2)
  ) +
    ggplot2::geom_point(ggplot2::aes(color = lake_names), size = 3) +
    ggplot2::xlab("Dim 1") +
    ggplot2::ylab("Dim 2") +
    ggplot2::labs(title = title) +
    ggplot2::scale_color_manual(values = lake_colors, name = "Lakes", drop = FALSE) +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 10.5),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      panel.background = ggplot2::element_blank()
    )
  
  return(plot)
}

# Function to generate dissimilarity matrices and plot PCoA
generate_and_plot <- function(comm, species_name, prefix) {
  # Dissimilarity matrices
  tax_q0_beta <- hillR::hill_taxa_parti_pairwise(comm, q = 0)
  tax_q0_beta$beta_diss <- 1 - tax_q0_beta$region_similarity
  
  tax_q1_beta <- hillR::hill_taxa_parti_pairwise(comm, q = 1)
  tax_q1_beta$beta_diss <- 1 - tax_q1_beta$region_similarity
  
  # Keep only beta_diss
  tax_q0_beta <- tax_q0_beta %>%
    dplyr::select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
    dplyr::rename(sample_b = site1, sample_a = site2, taxo_q0 = beta_diss)
  
  tax_q1_beta <- tax_q1_beta %>%
    dplyr::select(-q, -TD_gamma, -TD_alpha, -TD_beta, -local_similarity, -region_similarity) %>%
    dplyr::rename(sample_b = site1, sample_a = site2, taxo_q1 = beta_diss)
  
  # Transform dissimilarity indices into matrix
  taxo_q0_beta <- as.data.frame(tax_q0_beta)
  taxo_q1_beta <- as.data.frame(tax_q1_beta)
  
  samples <- unique(c(taxo_q0_beta$sample_a, taxo_q0_beta$sample_b))
  tax_q0_beta_matrix <- matrix(0, nrow = length(samples), ncol = length(samples),
                               dimnames = list(samples, samples))
  for (i in 1:nrow(taxo_q0_beta)) {
    sample1  <- taxo_q0_beta[i, "sample_a"]
    sample2  <- taxo_q0_beta[i, "sample_b"]
    distance <- taxo_q0_beta[i, "taxo_q0"]
    tax_q0_beta_matrix[sample1, sample2] <- distance
    tax_q0_beta_matrix[sample2, sample1] <- distance
  }
  
  samples <- unique(c(taxo_q1_beta$sample_a, taxo_q1_beta$sample_b))
  tax_q1_beta_matrix <- matrix(0, nrow = length(samples), ncol = length(samples),
                               dimnames = list(samples, samples))
  for (i in 1:nrow(taxo_q1_beta)) {
    sample1  <- taxo_q1_beta[i, "sample_a"]
    sample2  <- taxo_q1_beta[i, "sample_b"]
    distance <- taxo_q1_beta[i, "taxo_q1"]
    tax_q1_beta_matrix[sample1, sample2] <- distance
    tax_q1_beta_matrix[sample2, sample1] <- distance
  }
  
  # Generate PCoA plots
  plot_q0 <- plot_pcoa(tax_q0_beta_matrix,
                       paste(species_name, "- Annotated metabolite dissimilarity (q0)"))
  plot_q1 <- plot_pcoa(tax_q1_beta_matrix,
                       paste(species_name, "- Annotated metabolite dissimilarity (q1)"))
  
  # Return plots
  return(list(plot_q0, plot_q1))
}

# Apply functions
species_data <- list(
  list(species_name = "Perca fluviatilis", prefix = "per", comm = metabo_per),
  list(species_name = "Lepomis gibbosus",  prefix = "lep", comm = metabo_lep)
)

plots <- list()
for (species in species_data) {
  species_name <- species$species_name
  prefix       <- species$prefix
  comm         <- species$comm
  
  plots[[length(plots) + 1]] <- generate_and_plot(comm, species_name, prefix)
}

# Unpack plots
plot_taxo_q0_annotated_metabo_per <- plots[[1]][[1]]
plot_taxo_q1_annotated_metabo_per <- plots[[1]][[2]]
plot_taxo_q0_annotated_metabo_lep <- plots[[2]][[1]]
plot_taxo_q1_annotated_metabo_lep <- plots[[2]][[2]]

# Combine plots
pcoa_diss_annotated_metabo_sp <- ggpubr::ggarrange(
  plot_taxo_q0_annotated_metabo_per, plot_taxo_q1_annotated_metabo_per,
  plot_taxo_q0_annotated_metabo_lep, plot_taxo_q1_annotated_metabo_lep,
  labels = c("A", "B", "C", "D"),
  ncol = 2, nrow = 2,
  common.legend = TRUE, legend = "right"
)

# Save combined plot
path_to_my_object <- here::here("figures", "figSM4.png")
ggplot2::ggsave(
  filename = path_to_my_object,
  plot = pcoa_diss_annotated_metabo_sp,
  device = "png",
  width = 297,
  height = 210,
  units = "mm",
  dpi = 300
)
