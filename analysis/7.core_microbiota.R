#======== PROJECT COM2LIFE ========
## Estimates beta diversity on core microbiota of indivduals.
## Delete ASVs that are in less than 10% of samples

# Load libraries
library(phyloseq)
library(dplyr)
library(tidyverse)
library(hillR)
library(ggpubr)
library(ggplot2)
library(cowplot)

# load physeq object
physeq_filtered <- readRDS(here::here("data",
                             "mon_objet_physeq_filtered.rds"))

# Transform counts in relative abundance ici 
popfish_ra <- transform_sample_counts(physeq_filtered, function(x) x / sum(x))

# Calculate relative percentages for each ASVs in each sample
relative_abundance <- transform_sample_counts(popfish_ra, function(x) x / sum(x))

# delete ASVs that are in less than 10% of samples
filtered_abundance <- prune_taxa(taxa_sums(relative_abundance) >= 0.10, relative_abundance)

# cretae metadata object from filtered_abundance
metadata <- filtered_abundance@sam_data %>%
  `class<-`(NULL) %>%               
  `attr<-`("package", NULL) %>%     
  as.data.frame() 


# estimates beta diversity with Hill numbers
#'@taxonomic_diversity
comm = filtered_abundance@otu_table

# Beta diversity
#q = 0 (species richness)
tax_q0_beta <- hillR::hill_taxa_parti_pairwise(comm,q = 0) #output='matrix' for distance matric in output
tax_q0_beta$beta_diss <- 1 - tax_q0_beta$region_similarity

#q = 1 (shannon entropy)
tax_q1_beta <- hillR::hill_taxa_parti_pairwise(comm,q = 1) #output='matrix' for distance matric in output
tax_q1_beta$beta_diss <- 1 - tax_q1_beta$region_similarity


#'@Phylogenetic_diversity
tree= filtered_abundance@phy_tree

# Beta diversity
#q = 0 (species richness)
phylo_q0_beta <- hillR::hill_phylo_parti_pairwise(comm,tree,q = 0)#output='matrix' for distance matric in output
phylo_q0_beta$beta_diss <- 1 - phylo_q0_beta$region_similarity

#q = 1 (shannon entropy)
phylo_q1_beta <- hillR::hill_phylo_parti_pairwise(comm,tree,q = 1) #output='matrix' for distance matric in output
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

beta_diss <- phylo_q0_beta %>%
  merge(phylo_q1_beta, by = c('sample_a', 'sample_b')) %>%
  merge(tax_q0_beta, by = c('sample_a', 'sample_b')) %>%
  merge(tax_q1_beta, by = c('sample_a', 'sample_b'))

######################## ANALYSIS OF DISSIMILARITY #############################
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

## plot 
my_palette <- c("#41AEBD","#97E9D5","#F4DE3A", "#FCB11C","#A2CF49","#608F3D") 
region_colors <-  c(
  "CSM" = "#608F3D",
  "CERL" = "#41AEBD",
  "CERS" = "#97E9D5",
  "CRE" = "#F4DE3A",
  "LGP" = "#A2CF49",
  "VSS" = "#FCB11C"
)
custom_order <- c("CERL", "CERS", "CRE", "VSS", "LGP", "CSM")

# =================== PHYLO_Q0 ===================
# Filtrer les données pour LEP, inter_species et intra_region
filtered_data <- subset(beta_diss, species_a == "LEP" & species_b == "LEP" & sp_lev == "intra_species" & reg_lev == "intra_region")
filtered_data$region_a <- factor(filtered_data$region_a, levels = custom_order)
# Créer le graphique avec ggplot
a <- ggplot(filtered_data, aes(x = region_a, y = as.numeric(as.character(phylo_q0)), fill = region_a)) +
  geom_boxplot(color = "black", size = 0.2) + 
  labs(x = NULL, y = NULL) +
  ggtitle("Lepomis gibbosus") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 13.5, face = "bold"),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(y = bquote(bold("Au sein des lacs"))) + 
  scale_fill_manual(values = my_palette, name = "Lacs")  

# Filtrer les données pour LEP, inter_species et intra_region
filtered_data2 <- subset(beta_diss, species_a == "PER" & species_b == "PER" & sp_lev == "intra_species" & reg_lev == "intra_region")
filtered_data2$region_a <- factor(filtered_data2$region_a, levels = custom_order)
# Créer le graphique avec ggplot
b <- ggplot(filtered_data2, aes(x = region_a, y = as.numeric(as.character(phylo_q0)), fill = region_a)) +
  geom_boxplot(color = "black", size = 0.2) + 
  labs( x = NULL, y= NULL) +
  ggtitle("Perca fluviatilis") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12.5, face = "bold",  hjust = 0.5),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_fill_manual(values = my_palette, name = "Lacs")  


# Filtrer les données pour les interactions spécifiques
filtered_data3 <- subset(beta_diss, species_a == "LEP" & species_b == "LEP" & sp_lev == "intra_species" & reg_lev == "inter_region")
# Créer le graphique avec ggplot
c <- ggplot(filtered_data3, aes(x = reg_comb, y = as.numeric(as.character(phylo_q0))))  +
  geom_boxplot(size = 0.6, aes(fill = region_b, color = region_a)) +
  labs(x =NULL, y = NULL) +
  ggtitle("Lepomis gibbosus") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12.5, face = "bold",  hjust = 0.5),
        axis.title.y = element_text(size = 13.5, face = "bold"),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right")+
  scale_fill_manual(values = region_colors, name = "Lacs") +  
  scale_color_manual(values = region_colors, name = "Comparaison lacs") +  
  labs(y = bquote(bold("Inter-lacs"))) +  
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) 


# Filtrer les données pour les interactions spécifiques
filtered_data4 <- subset(beta_diss, species_a == "PER" & species_b == "PER" & sp_lev == "intra_species" & reg_lev == "inter_region")
# Créer le graphique avec ggplot
d <- ggplot(filtered_data4, aes(x = reg_comb, y = as.numeric(as.character(phylo_q0))))  +
  geom_boxplot(size = 0.6, aes(fill = region_b, color = region_a)) +
  labs(x =NULL, y = NULL) +
  ggtitle("Perca fluviatilis") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12.5, face = "bold",  hjust = 0.5),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right")+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_fill_manual(values = region_colors, name = "Lacs") +  
  scale_color_manual(values = region_colors, name = "Comparaison lacs")

# Créez la grille des graphiques avec des hauteurs spécifiées
phylo_0 <- plot_grid(
  a + theme(legend.position = "none"),
  b + theme(legend.position = "none"),
  c + theme(legend.position = "none"),
  d + theme(legend.position = "none"),  
  labels = c("A", "B", "C", "D"),
  ncol = 2
)

# Combinez le titre avec le reste du plot
boxplot_diss_phylo_q0 <- cowplot::plot_grid(
  phylo_0,
  ncol = 1,
  rel_heights = c(0.15, 20) 
)

path_to_my_object = here::here("Figures","hill", "CM_boxplot_diss_phylo_q0.png")
ggsave(filename = path_to_my_object, plot = boxplot_diss_phylo_q0, device = "png")


# =================== PHYLO_Q1 ===================
# Filtrer les données pour LEP, inter_species et intra_region
filtered_data <- subset(beta_diss, species_a == "LEP" & species_b == "LEP" & sp_lev == "intra_species" & reg_lev == "intra_region")
filtered_data$region_a <- factor(filtered_data$region_a, levels = custom_order)

# Créer le graphique avec ggplot
a <- ggplot(filtered_data, aes(x = region_a, y = as.numeric(as.character(phylo_q1)), fill = region_a)) +
  geom_boxplot(color = "black", size = 0.2) +
  labs(x = NULL, y = NULL) +
  ggtitle("Lepomis gibbosus") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 13.5, face = "bold"),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 1, by = 0.2)) +
  labs(y = bquote(bold("Au sein des lacs"))) +  
  scale_fill_manual(values = my_palette, name = "Lacs")

# Filtrer les données pour LEP, inter_species et intra_region
filtered_data2 <- subset(beta_diss, species_a == "PER" & species_b == "PER" & sp_lev == "intra_species" & reg_lev == "intra_region")
filtered_data2$region_a <- factor(filtered_data2$region_a, levels = custom_order)

# Créer le graphique avec ggplot
b <- ggplot(filtered_data2, aes(x = region_a, y = as.numeric(as.character(phylo_q1)), fill = region_a)) +
  geom_boxplot(color = "black", size = 0.2) +
  labs(x = NULL, y = NULL) +
  ggtitle("Perca fluviatilis") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 1, by = 0.2)) +
  scale_fill_manual(values = my_palette, name = "Lacs")

# Filtrer les données pour les interactions spécifiques
filtered_data3 <- subset(beta_diss, species_a == "LEP" & species_b == "LEP" & sp_lev == "intra_species" & reg_lev == "inter_region")

# Créer le graphique avec ggplot
c <- ggplot(filtered_data3, aes(x = reg_comb, y = as.numeric(as.character(phylo_q1)))) +
  geom_boxplot(size = 0.6, aes(fill = region_b, color = region_a)) +
  labs(x = NULL, y = NULL) +
  ggtitle("Lepomis gibbosus") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 13.5, face = "bold"),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 1, by = 0.2)) +
  labs(y = bquote(bold("Inter-lacs"))) +  
  scale_fill_manual(values = region_colors, name = "Lacs") + 
  scale_color_manual(values = region_colors, name = "Comparaison lacs")


# Filtrer les données pour les interactions spécifiques
filtered_data4 <- subset(beta_diss, species_a == "PER" & species_b == "PER" & sp_lev == "intra_species" & reg_lev == "inter_region")

# Créer le graphique avec ggplot
d <- ggplot(filtered_data4, aes(x = reg_comb, y = as.numeric(as.character(phylo_q1)))) +
  geom_boxplot(size = 0.6, aes(fill = region_b, color = region_a)) +
  labs(x = NULL, y = NULL) +
  ggtitle("Perca fluviatilis") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 1, by = 0.2)) +
  scale_fill_manual(values = region_colors, name = "Lacs") + 
  scale_color_manual(values = region_colors, name = "Comparaison lacs") 

# Créez la grille des graphiques avec des hauteurs spécifiées
phylo_1 <- cowplot::plot_grid(
  a + theme(legend.position = "none"),
  b + theme(legend.position = "none"),
  c + theme(legend.position = "none"),
  d + theme(legend.position = "none"),  # Conserve la légende
  labels = c("A", "B", "C", "D"),
  ncol = 2
)

# Combinez le titre avec le reste du plot
boxplot_diss_phylo_q1 <- cowplot::plot_grid(
  phylo_1,
  ncol = 1,
  rel_heights = c(0.15, 20)  
)

path_to_my_object = here::here("Figures","hill", "CM_boxplot_diss_phylo_q1.png")
ggsave(filename = path_to_my_object, plot = boxplot_diss_phylo_q1, device = "png")


# =================== TAXO_q0 ===================
# Filtrer les données pour LEP, inter_species et intra_region
filtered_data <- subset(beta_diss, species_a == "LEP" & species_b == "LEP" & sp_lev == "intra_species" & reg_lev == "intra_region")
filtered_data$region_a <- factor(filtered_data$region_a, levels = custom_order)

# Créer le graphique avec ggplot
a <- ggplot(filtered_data, aes(x = region_a, y = as.numeric(as.character(taxo_q0)), fill = region_a)) +
  geom_boxplot(color = "black", size = 0.2) +
  labs(x = NULL, y = NULL) +
  ggtitle("Lepomis gibbosus") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 13.5, face = "bold"),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(y = bquote(bold("Au sein des lacs"))) +  
  scale_fill_manual(values = my_palette, name = "Lacs")

# Filtrer les données pour LEP, inter_species et intra_region
filtered_data2 <- subset(beta_diss, species_a == "PER" & species_b == "PER" & sp_lev == "intra_species" & reg_lev == "intra_region")
filtered_data2$region_a <- factor(filtered_data2$region_a, levels = custom_order)

# Créer le graphique avec ggplot
b <-  ggplot(filtered_data2, aes(x = region_a, y = as.numeric(as.character(taxo_q0)), fill = region_a)) +
  geom_boxplot(color = "black", size = 0.2) +
  labs(x = NULL, y = NULL) +
  ggtitle("Perca fluviatilis") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_fill_manual(values = my_palette, name = "Lacs")

# Filtrer les données pour les interactions spécifiques
filtered_data3 <- subset(beta_diss, species_a == "LEP" & species_b == "LEP" & sp_lev == "intra_species" & reg_lev == "inter_region")
filtered_data3$region_a <- factor(filtered_data3$region_a, levels = custom_order)

# Créer le graphique avec ggplot
c <-  ggplot(filtered_data3, aes(x = reg_comb, y = as.numeric(as.character(taxo_q0)))) +
  geom_boxplot(size = 0.6, aes(fill = region_b, color = region_a)) +
  labs(x = NULL, y = NULL) +
  ggtitle("Lepomis gibbosus") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 13.5, face = "bold"),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  scale_y_continuous(limits = c(0.4, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(y = bquote(bold("Inter-lacs"))) +  
  scale_fill_manual(values = region_colors, name = "Lacs") + 
  scale_color_manual(values = region_colors, name = "Comparaison lacs")

# Filtrer les données pour les interactions spécifiques
filtered_data4 <- subset(beta_diss, species_a == "PER" & species_b == "PER" & sp_lev == "intra_species" & reg_lev == "inter_region")
filtered_data4$region_a <- factor(filtered_data4$region_a, levels = custom_order)

# Créer le graphique avec ggplot
d <- ggplot(filtered_data4, aes(x = reg_comb, y = as.numeric(as.character(taxo_q0)))) +
  geom_boxplot(size = 0.6, aes(fill = region_b, color = region_a)) +
  labs(x = NULL, y = NULL) +
  ggtitle("Perca fluviatilis") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  scale_y_continuous(limits = c(0.4, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_fill_manual(values = region_colors, name = "Lacs") + 
  scale_color_manual(values = region_colors, name = "Comparaison lacs")

# Organiser les graphiques dans une grille avec un titre commun
# Créez la grille des graphiques avec des hauteurs spécifiées
taxo_0 <- cowplot::plot_grid(
  a + theme(legend.position = "none"),
  b + theme(legend.position = "none"),
  c + theme(legend.position = "none"),
  d + theme(legend.position = "none"),  
  labels = c("A", "B", "C", "D"),
  ncol = 2
)

# Combinez le titre avec le reste du plot
boxplot_diss_taxo_q0 <- cowplot::plot_grid(
  taxo_0,
  ncol = 1,
  rel_heights = c(0.15, 20)  
)

path_to_my_object = here::here("Figures","hill", "CM_boxplot_diss_taxo_q0.png")
ggsave(filename = path_to_my_object, plot = boxplot_diss_taxo_q0, device = "png")


# =================== TAXO_q1 ===================
# Filtrer les données pour LEP, inter_species et intra_region
filtered_data <- subset(beta_diss, species_a == "LEP" & species_b == "LEP" & sp_lev == "intra_species" & reg_lev == "intra_region")
filtered_data$region_a <- factor(filtered_data$region_a, levels = custom_order)

# Créer le graphique avec ggplot
a <- ggplot(filtered_data, aes(x = region_a, y = as.numeric(as.character(taxo_q1)), fill = region_a)) +
  geom_boxplot(color = "black", size = 0.2) +
  labs(x = NULL, y = NULL) +
  ggtitle("Lepomis gibbosus") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 13.5, face = "bold"),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(y = bquote(bold("Au sein des lacs"))) +  
  scale_fill_manual(values = my_palette, name = "Lacs")

# Filtrer les données pour LEP, inter_species et intra_region
filtered_data2 <- subset(beta_diss, species_a == "PER" & species_b == "PER" & sp_lev == "intra_species" & reg_lev == "intra_region")
filtered_data2$region_a <- factor(filtered_data2$region_a, levels = custom_order)

# Créer le graphique avec ggplot
b <- ggplot(filtered_data2, aes(x = region_a, y = as.numeric(as.character(taxo_q1)), fill = region_a)) +
  geom_boxplot(color = "black", size = 0.2) +
  labs(x = NULL, y = NULL) +
  ggtitle("Perca fluviatilis") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_fill_manual(values = my_palette, name = "Lacs")

# Filtrer les données pour les interactions spécifiques
filtered_data3 <- subset(beta_diss, species_a == "LEP" & species_b == "LEP" & sp_lev == "intra_species" & reg_lev == "inter_region")
filtered_data3$region_a <- factor(filtered_data3$region_a, levels = custom_order)

# Créer le graphique avec ggplot
c <- ggplot(filtered_data3, aes(x = reg_comb, y = as.numeric(as.character(taxo_q1)))) +
  geom_boxplot(size = 0.6, aes(fill = region_b, color = region_a)) +
  labs(x = NULL, y = NULL) +
  ggtitle("Lepomis gibbosus") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 13.5, face = "bold"),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(y = bquote(bold("Inter-lacs"))) +  
  scale_fill_manual(values = region_colors, name = "Lacs") + 
  scale_color_manual(values = region_colors, name = "Comparaison lacs") 

# Filtrer les données pour les interactions spécifiques
filtered_data4 <- subset(beta_diss, species_a == "PER" & species_b == "PER" & sp_lev == "intra_species" & reg_lev == "inter_region")
filtered_data4$region_a <- factor(filtered_data4$region_a, levels = custom_order)

# Créer le graphique avec ggplot
d <- ggplot(filtered_data4, aes(x = reg_comb, y = as.numeric(as.character(taxo_q1)))) +
  geom_boxplot(size = 0.6, aes(fill = region_b, color = region_a)) +
  labs(x = NULL, y = NULL) +
  ggtitle("Perca fluviatilis") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 11),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  scale_y_continuous(limits = c(0.2, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_fill_manual(values = region_colors, name = "Lacs") + 
  scale_color_manual(values = region_colors, name = "Comparaison lacs")

# Organiser les graphiques dans une grille avec un titre commun
taxo_1 <- cowplot::plot_grid(
  a + theme(legend.position = "none"),
  b + theme(legend.position = "none"),
  c + theme(legend.position = "none"),
  d + theme(legend.position = "none"),  
  labels = c("A", "B", "C", "D"),
  ncol = 2
)

# Combinez le titre avec le reste du plot
boxplot_diss_taxo_q1 <- cowplot::plot_grid(
  taxo_1,
  ncol = 1,
  rel_heights = c(0.15, 20) 
)

path_to_my_object = here::here("Figures","hill", "CM_boxplot_diss_taxo_q1.png")
ggsave(filename = path_to_my_object, plot = boxplot_diss_taxo_q1, device = "png")


# ================ INTER SPECIES INTRA REGION : PHYLO_Q0, PHYLO_Q1, TAXO_Q0, TAXO_Q1 ===========================
filtered_data5 <- subset(beta_diss, sp_lev == "inter_species" & reg_lev == "intra_region")
filtered_data5$region_a <- factor(filtered_data5$region_a, levels = custom_order)
# Créer le graphique avec ggplot
e <- ggplot(filtered_data5, aes(x = region_a, y = as.numeric(as.character(phylo_q0)))) +
  geom_boxplot(fill = my_palette, color = "black") +
  labs(x = NULL, y = NULL) +
  ggtitle("Phylogenetic dissimilarity (q0)") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, face = "bold"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(values = my_palette)



f <- ggplot(filtered_data5, aes(x = region_a, y = as.numeric(as.character(phylo_q1)))) +
  geom_boxplot(fill = my_palette, color = "black") +
  labs(x = NULL, y = NULL) +
  ggtitle("Phylogenetic dissimilarity (q1)") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, face = "bold"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(values = my_palette)


g <- ggplot(filtered_data5, aes(x = region_a, y = as.numeric(as.character(taxo_q0)))) +
  geom_boxplot(fill = my_palette, color = "black") +
  labs(x = NULL, y = NULL) +
  ggtitle("Taxonomic dissimilarity (q0)") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, face = "bold"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(values = my_palette)


h <- ggplot(filtered_data5, aes(x = region_a, y = as.numeric(as.character(taxo_q1)))) +
  geom_boxplot(fill = my_palette, color = "black") +
  labs(x = NULL, y = NULL) +
  ggtitle("Taxonomic dissimilarity (q1)") +
  theme_minimal() + 
  theme(plot.title = element_text(size = 10, face = "bold"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(values = my_palette)

# Organiser les graphiques dans une grille avec un titre commun
inter_species <- plot_grid(g, h, e, f, labels = c("A", "B", "C", "D"), ncol = 2)


# Combinez le titre avec le reste du plot
boxplot_diss_inter_sp <- cowplot::plot_grid(
  inter_species,
  ncol = 1,
  rel_heights = c(0.15, 20)  
)

path_to_my_object = here::here("Figures","hill", "CM_boxplot_diss_inter_sp.png")
ggsave(filename = path_to_my_object, plot = boxplot_diss_inter_sp, device = "png")
