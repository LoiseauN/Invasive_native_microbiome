#======== PROJECT COM2LIFE ========
## Estimates beta diversity based on hill numbers for metabolite data  
## First load the libraries and physeq_object. Calculates beta diversity on phylogenetic (qo, q1) taxonomic (q0, q1)

## =============== PCOA and DbRDA ================== ##
# Load objects 
physeq_filtered <- readRDS(here::here("data",
                                      "mon_objet_physeq_filtered.rds"))
phylo_q0_beta <- readRDS(here::here("data",
                                    "phylo_q0_beta.rds"))
phylo_q1_beta <- readRDS(here::here("data",
                                    "phylo_q1_beta.rds"))
taxo_q0_beta <- readRDS(here::here("data",
                                   "tax_q0_beta.rds"))
taxo_q1_beta <- readRDS(here::here("data",
                                   "tax_q1_beta.rds"))



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
# Form more CHLa to less CSM LGP VSS CRE CERS CERL
pal <- harrypotter::hp(n = 6, option = "DracoMalfoy")

lake_colors <- c(
  CSM ="#00332AFF",
  LGP = "#035236FF",
  VSS ="#227548FF",
  CRE ="#60A360FF",
  CERS = "#C5D163FF",
  CERL = "#F2F26DFF")

  #lake_colors <- c(
#  "CSM" = "#608F3D",
#  "CERL" = "#41AEBD",
#  "CERS" = "#97E9D5",
#  "CRE" = "#F4DE3A",
#  "LGP" = "#A2CF49",
#  "VSS" = "#FCB11C"
#)
plot_pcoa <- function(matrix_dist, title, lake_order) {
  sample_names <- rownames(matrix_dist)
  lake_names <- substr(sample_names, 1, 4)
  lake_names[lake_names == "CRE."] <- "CRE"
  lake_names[lake_names == "VSS."] <- "VSS"
  lake_names[lake_names == "LGP."] <- "LGP"
  lake_names[lake_names == "CSM."] <- "CSM"
  
  # Order lake_names according to the chosen level order
  lake_names <- factor(lake_names, levels = lake_order)
  
  # Perform PCoA
  pcoa_result <- ape::pcoa(matrix_dist)
  pcoa_result$lake_names <- lake_names
  pcoa_coordinates <- pcoa_result$vectors
  
  point_shapes <- ifelse(grepl("PER", rownames(matrix_dist)), 16, 1)
  
  plot <- ggplot2::ggplot(data = pcoa_coordinates, ggplot2::aes(x = Axis.1, y = Axis.2)) +
    ggplot2::geom_point(ggplot2::aes(color = lake_names, shape = factor(point_shapes)), size = 3) +
    ggplot2::xlab("Dim 1") +
    ggplot2::ylab("Dim 2") +
    ggplot2::labs(title = title) +
    ggplot2::scale_color_manual(values = lake_colors, name = "Lakes", limits = lake_order) +  
    ggplot2::scale_shape_manual(values = c(16, 1),
                                labels = c("Perca fluviatilis", "Lepomis gibbosus")) +  
    ggplot2::guides(shape = ggplot2::guide_legend(title = "Species",     
                                                  override.aes = list(shape = c(16, 1)))) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"), 
      panel.background = ggplot2::element_blank()
    )
  
  return(plot)
}

lake_order <- c("CSM", "LGP", "VSS", "CRE", "CERS", "CERL")  # your chosen order

plot_taxo_q0_pcoa   <- plot_pcoa(tax_q0_beta_matrix,   "Taxonomic dissimilarity (q0)",   lake_order)
plot_taxo_q1_pcoa   <- plot_pcoa(tax_q1_beta_matrix,   "Taxonomic dissimilarity (q1)",   lake_order)
plot_phylo_q0_pcoa  <- plot_pcoa(phylo_q0_beta_matrix, "Phylogenetic dissimilarity (q0)",lake_order)
plot_phylo_q1_pcoa  <- plot_pcoa(phylo_q1_beta_matrix, "Phylogenetic dissimilarity (q1)",lake_order)



## Add metabolo
taxo_q0_beta <- as.data.frame(tax_q0_beta_metabo)


samples <- unique(c(taxo_q0_beta$sample_a, taxo_q0_beta$sample_b))
tax_q0_beta_matrix <- matrix(0, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))
for (i in 1:nrow(taxo_q0_beta)) {
  sample1 <- taxo_q0_beta[i, "sample_a"]
  sample2 <- taxo_q0_beta[i, "sample_b"]
  distance <- taxo_q0_beta[i, "taxo_q0"]
  tax_q0_beta_matrix[sample1, sample2] <- distance
  tax_q0_beta_matrix[sample2, sample1] <- distance
}


# transform the dissimilarity indice in matrix
tax_q0_beta_metabo <- readRDS(here::here("Data","tax_q0_beta_metabo.rds"))
taxo_q0_beta <- as.data.frame(tax_q0_beta_metabo)


samples <- unique(c(taxo_q0_beta$sample_a, taxo_q0_beta$sample_b))
tax_q0_beta_matrix <- matrix(0, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))
for (i in 1:nrow(taxo_q0_beta)) {
  sample1 <- taxo_q0_beta[i, "sample_a"]
  sample2 <- taxo_q0_beta[i, "sample_b"]
  distance <- taxo_q0_beta[i, "taxo_q0"]
  tax_q0_beta_matrix[sample1, sample2] <- distance
  tax_q0_beta_matrix[sample2, sample1] <- distance
}


# transform the dissimilarity indice in matrix
taxo_q0_beta_metabo <- as.data.frame(tax_q0_beta_matrix)

rownames(taxo_q0_beta_metabo) <- gsub("_G[0-9]+_", "_GAR_", rownames(taxo_q0_beta_metabo))
rownames(taxo_q0_beta_metabo) <- gsub("_P[0-9]+_", "_PER_", rownames(taxo_q0_beta_metabo))

colnames(taxo_q0_beta_metabo) <- gsub("_G[0-9]+_", "_GAR_", colnames(taxo_q0_beta_metabo))
colnames(taxo_q0_beta_metabo) <- gsub("_P[0-9]+_", "_PER_", colnames(taxo_q0_beta_metabo))


plot_pcoa_metabo <- function(matrix_dist, title, lake_order) {
  
  # matrix_dist =  tax_q0_beta_matrix_metabo
  
  
  sample_names <- rownames(matrix_dist)
  lake_names <- substr(sample_names, 1, 4)
  lake_names[lake_names == "CRE_"] <- "CRE"
  lake_names[lake_names == "VSS_"] <- "VSS"
  lake_names[lake_names == "LGP_"] <- "LGP"
  lake_names[lake_names == "CSM_"] <- "CSM"
  
  # Order lake_names according to the chosen level order
  lake_names <- factor(lake_names, levels = lake_order)
  
  # Perform PCoA
  pcoa_result <- ape::pcoa(matrix_dist)
  pcoa_result$lake_names <- lake_names
  pcoa_coordinates <- pcoa_result$vectors
  
  point_shapes <- ifelse(grepl("PER", rownames(matrix_dist)), 16, 1)
  
  plot <- ggplot2::ggplot(data = pcoa_coordinates, ggplot2::aes(x = Axis.1, y = Axis.2)) +
    ggplot2::geom_point(ggplot2::aes(color = lake_names, shape = factor(point_shapes)), size = 3) +
    ggplot2::xlab("Dim 1") +
    ggplot2::ylab("Dim 2") +
    ggplot2::labs(title = title) +
    ggplot2::scale_color_manual(values = lake_colors, name = "Lakes", limits = lake_order) +  
    ggplot2::scale_shape_manual(values = c(16, 1),
                                labels = c("Perca fluviatilis", "Lepomis gibbosus")) +  
    ggplot2::guides(shape = ggplot2::guide_legend(title = "Species",     
                                                  override.aes = list(shape = c(16, 1)))) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"), 
      panel.background = ggplot2::element_blank()
    )
  
  return(plot)
}

lake_order <- c("CSM", "LGP", "VSS", "CRE", "CERS", "CERL")  # your chosen order

plot_metabo_q0_pcoa   <- plot_pcoa_metabo(taxo_q0_beta_metabo,   "Metabolic dissimilarity (q0)",   lake_order)


## METABO 
#Combine plots
# Extract one legend (from a plot that contains all the legend entries)
legend_common <- ggpubr::get_legend(plot_taxo_q0_pcoa)

# Top: A-D in a 2 x 2 grid, no legends
top <- ggpubr::ggarrange(
  plot_taxo_q0_pcoa, plot_taxo_q1_pcoa,
  plot_phylo_q0_pcoa, plot_phylo_q1_pcoa,
  labels = c("A", "B", "C", "D"),
  ncol = 2, nrow = 2,
  legend = "none"
)

# Bottom: E centered (empty panels on each side)
bottom <- ggpubr::ggarrange(
  NULL, plot_metabo_q0_pcoa, NULL,
  labels = c("", "E", ""),
  ncol = 3, nrow = 1,
  widths = c(1, 2, 1),      # E = 2/4 of the width, same as one panel above
  legend = "none"
)

# Combine the grid and the row of E (2 rows for A-D, 1 row for E)
figure <- ggpubr::ggarrange(top, bottom, ncol = 1, nrow = 2, heights = c(2, 1))

# Add the common legend on the right
pcoa_diss <- ggpubr::ggarrange(
  figure, ggpubr::as_ggplot(legend_common),
  ncol = 2, widths = c(1, 0.15)
)

path_to_my_object = here::here("figures", "figure2_All.png")
ggplot2::ggsave(filename = path_to_my_object, plot = pcoa_diss, device = "png",
                width = 10, height = 7)




