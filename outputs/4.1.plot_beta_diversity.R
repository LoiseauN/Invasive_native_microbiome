#======== PROJECT COM2LIFE ========
## Estimates beta diversity based on hill numbers 
## load the object and plot the beta diversity 

# load objects
beta_diss <- readRDS(here::here("Data",
                                      "beta_diss.rds"))
metadata <- readRDS(here::here("Data",
                                      "metadata.rds"))

## load libraries
library(ggpubr)
library(ggplot2)
library(cowplot)

#parameters for boxplot
species_list <- c("LEP", "PER")
species_labels <- c("LEP" = "Lepomis gibbosus", "PER" = "Perca fluviatilis")
sp_level <- "intra_species"
reg_levels <- c("intra_region", "inter_region")
reg_comb_var <- "reg_comb"
my_palette <- c("#41AEBD","#97E9D5","#F4DE3A", "#FCB11C","#A2CF49","#608F3D") 
region_colors <- c(
  "CSM" = "#608F3D",
  "CERL" = "#41AEBD",
  "CERS" = "#97E9D5",
  "CRE" = "#F4DE3A",
  "LGP" = "#A2CF49",
  "VSS" = "#FCB11C"
)
custom_order <- c("CERL", "CERS", "CRE", "VSS", "LGP", "CSM")

# Function to generate phylogenetic diversity plots for different species and regions
diversity_plots <- function(data, species_list, species_labels, indice_var, sp_level, reg_levels, reg_comb_var, custom_order, my_palette, region_colors) {
  
  # List
  plot_list_intra <- list()
  plot_list_inter <- list()
  
  # Loop
  for (species in species_list) {
    
    # Filter data 
    filtered_data_intra <- subset(data, species_a == species & species_b == species & sp_lev == sp_level & reg_lev == reg_levels[1])
    filtered_data_intra$region_a <- factor(filtered_data_intra$region_a, levels = custom_order)
    
    # Create the plot
    plot_intra <- ggplot(filtered_data_intra, aes(x = region_a, y = as.numeric(as.character(get(indice_var))), fill = region_a)) +
      geom_boxplot(color = "black", size = 0.2) + 
      labs(x = NULL, y = NULL) +
      ggtitle(species_labels[species]) +
      theme_minimal() + 
      theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
            axis.text.y = element_text(size = 11),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position = "none") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      labs(y = bquote(bold("Within lakes"))) + 
      scale_fill_manual(values = my_palette, name = "Lacs")
    
    plot_list_intra[[species]] <- plot_intra
    
    # Filter data to keep inter region
    filtered_data_inter <- subset(data, species_a == species & species_b == species & sp_lev == sp_level & reg_lev == reg_levels[2])
    
    # Create the plot
    plot_inter <- ggplot(filtered_data_inter, aes(x = get(reg_comb_var), y = as.numeric(as.character(get(indice_var))))) +
      geom_boxplot(size = 0.6, aes(fill = region_b, color = region_a)) +
      labs(x = NULL, y = NULL) +
      ggtitle(species_labels[species]) +
      theme_minimal() +
      theme(plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5),
            axis.text.y = element_text(size = 11),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            legend.position = "none") +
      scale_fill_manual(values = region_colors, name = "Lacs") +  
      scale_color_manual(values = region_colors, name = "Comparaison lacs") +
      labs(y = bquote(bold("Between lakes"))) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
    
    plot_list_inter[[species]] <- plot_inter
  }
  
  # Combined plot intra region
  intra_row <- plot_grid(
    plot_list_intra[[species_list[1]]], plot_list_intra[[species_list[2]]],
    labels = c("A", "B"),
    ncol = 2
  )
  
  # Combine plot inter region
  inter_row <- plot_grid(
    plot_list_inter[[species_list[1]]], plot_list_inter[[species_list[2]]],
    labels = c("C", "D"),
    ncol = 2
  )
  
  combined_plots <- plot_grid(
    intra_row,
    inter_row,
    ncol = 1,
    rel_heights = c(1, 1) 
  )
  
  return(combined_plots)
}


# =================== TAXO_q0 ===================
indice_var <- "taxo_q0"
# Generate the combined plot for taxonomic q0 diversity
boxplot_diss_taxo_q0 <- diversity_plots(beta_diss, species_list, species_labels, indice_var, sp_level, reg_levels, reg_comb_var, custom_order, my_palette, region_colors)

path_to_my_object = here::here("Figures","hill", "boxplot_diss_taxo_q0.png")
ggsave(filename = path_to_my_object, plot = boxplot_diss_taxo_q0, device = "png")

# =================== TAXO_q1 ===================
indice_var <- "taxo_q1"
# Generate the combined plot for taxonomic q1 diversity
boxplot_diss_taxo_q1 <- diversity_plots(beta_diss, species_list, species_labels, indice_var, sp_level, reg_levels, reg_comb_var, custom_order, my_palette, region_colors)

path_to_my_object = here::here("Figures","hill", "boxplot_diss_taxo_q1.png")
ggsave(filename = path_to_my_object, plot = boxplot_diss_taxo_q1, device = "png")

# =================== PHYLO_Q0 ===================
indice_var <- "phylo_q0"
# Generate the combined plot for phylogenetic q0 diversity
boxplot_diss_phylo_q0 <- diversity_plots(beta_diss, species_list, species_labels, indice_var, sp_level, reg_levels, reg_comb_var, custom_order, my_palette, region_colors)

path_to_my_object = here::here("Figures","hill", "boxplot_diss_phylo_q0.png")
ggsave(filename = path_to_my_object, plot = boxplot_diss_phylo_q0, device = "png")

# =================== PHYLO_Q1 ===================
indice_var <- "phylo_q1"
# Generate the combined plot for phylogenetic q1 diversity
boxplot_diss_phylo_q1 <- diversity_plots(beta_diss, species_list, species_labels, indice_var, sp_level, reg_levels, reg_comb_var, custom_order, my_palette, region_colors)

path_to_my_object = here::here("Figures","hill", "boxplot_diss_phylo_q1.png")
ggsave(filename = path_to_my_object, plot = boxplot_diss_phylo_q1, device = "png")



# ================ INTER SPECIES INTRA REGION : PHYLO_Q0, PHYLO_Q1, TAXO_Q0, TAXO_Q1 ===========================
#Filter the data to keep the comparison between species and intra region
filtered_data5 <- subset(beta_diss, sp_lev == "inter_species" & reg_lev == "intra_region")
filtered_data5$region_a <- factor(filtered_data5$region_a, levels = custom_order)
# =================== PHYLO_Q0 ===================
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

# =================== PHYLO_Q1 ===================
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

# =================== TAXO_q0 ===================
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

# =================== TAXO_q1 ===================
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


inter_species <- plot_grid(g, h, e, f, labels = c("A", "B", "C", "D"), ncol = 2)

boxplot_diss_inter_sp <- cowplot::plot_grid(
  inter_species,
  ncol = 1,
  rel_heights = c(0.15, 20)  
)

path_to_my_object = here::here("Figures","hill", "boxplot_diss_inter_sp.png")
ggsave(filename = path_to_my_object, plot = boxplot_diss_inter_sp, device = "png")
