#======== PROJECT COM2LIFE ========
## Estimates beta diversity based on hill numbers 
## load the object and plot the beta diversity 

# load objects
beta_diss <- readRDS(here::here("data",
                                "beta_diss.rds"))
metadata <- readRDS(here::here("data",
                               "metadata.rds"))

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
    plot_intra <- ggplot2::ggplot(filtered_data_intra, ggplot2::aes(x = region_a, y = as.numeric(as.character(get(indice_var))), fill = region_a)) +
      ggplot2::geom_boxplot(color = "black", size = 0.2) + 
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::ggtitle(species_labels[species]) +
      ggplot2::theme_minimal() + 
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12.5, face = "bold", hjust = 0.5),
                     axis.text.x = ggplot2::element_text(size = 11, angle = 45, vjust = 0.5),
                     axis.text.y = ggplot2::element_text(size = 11),
                     panel.border = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black"),
                     legend.position = "none") +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      ggplot2::labs(y = bquote(bold("Within lakes"))) + 
      ggplot2::scale_fill_manual(values = my_palette, name = "Lacs")
    
    plot_list_intra[[species]] <- plot_intra
    
    # Filter data to keep inter region
    filtered_data_inter <- subset(data, species_a == species & species_b == species & sp_lev == sp_level & reg_lev == reg_levels[2])
    
    # Create the plot
    plot_inter <- ggplot2::ggplot(filtered_data_inter, ggplot2::aes(x = get(reg_comb_var), y = as.numeric(as.character(get(indice_var))))) +
      ggplot2::geom_boxplot(size = 0.6, ggplot2::aes(fill = region_b, color = region_a)) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::ggtitle(species_labels[species]) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12.5, face = "bold", hjust = 0.5),
                     axis.text.x = ggplot2::element_text(size = 11, angle = 45, vjust = 0.5),
                     axis.text.y = ggplot2::element_text(size = 11),
                     panel.border = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black"),
                     legend.position = "none") +
      ggplot2::scale_fill_manual(values = region_colors, name = "Lacs") +  
      ggplot2::scale_color_manual(values = region_colors, name = "Comparaison lacs") +
      ggplot2::labs(y = bquote(bold("Between lakes"))) +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
    
    plot_list_inter[[species]] <- plot_inter
  }
  
  # Combined plot intra region
  intra_row <- cowplot::plot_grid(
    plot_list_intra[[species_list[1]]], plot_list_intra[[species_list[2]]],
    labels = c("A", "B"),
    ncol = 2
  )
  
  # Combine plot inter region
  inter_row <- cowplot::plot_grid(
    plot_list_inter[[species_list[1]]], plot_list_inter[[species_list[2]]],
    labels = c("C", "D"),
    ncol = 2
  )
  
  combined_plots <- cowplot::plot_grid(
    intra_row,
    inter_row,
    ncol = 1,
    rel_heights = c(1, 1) 
  )
  
  return(combined_plots)
}


# =================== TAXO_q1 ===================
indice_var <- "taxo_q1"
# Generate the combined plot for taxonomic q1 diversity
boxplot_diss_taxo_q1 <- diversity_plots(beta_diss, species_list, species_labels, indice_var, sp_level, reg_levels, reg_comb_var, custom_order, my_palette, region_colors)

path_to_my_object = here::here("figures","figure3AB.png")
ggplot2::ggsave(filename = path_to_my_object, plot = boxplot_diss_taxo_q1, device = "png")

# ================ INTER SPECIES INTRA REGION : TAXO_Q1 ===========================
#Filter the data to keep the comparison between species and intra region
filtered_data5 <- subset(beta_diss, sp_lev == "inter_species" & reg_lev == "intra_region")
filtered_data5$region_a <- factor(filtered_data5$region_a, levels = custom_order)

# =================== TAXO_q1 ===================
h <- ggplot2::ggplot(filtered_data5, ggplot2::aes(x = region_a, y = as.numeric(as.character(taxo_q1)))) +
  ggplot2::geom_boxplot(fill = my_palette, color = "black") +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::ggtitle("Taxonomic dissimilarity (q1)") +
  ggplot2::theme_minimal() + 
  ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold"),
                 panel.border = ggplot2::element_blank(),
                 pantel.grid.major = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(colour = "black")) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  ggplot2::scale_color_manual(values = my_palette)

path_to_my_object = here::here("figures","figure3C.png")
ggplot2::ggsave(filename = path_to_my_object, plot = h, device = "png")
