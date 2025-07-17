#======== PROJECT COM2LIFE ========
## Estimates alpha diversity based on hill numbers 
## load the object and plot the alpha diversity

#create a directory
dir.create("Figures/hill", recursive = TRUE)

# Load my physeq object 
physeq_alpha <- readRDS(here::here("Data",
                                      "mon_objet_physeq_alpha.rds"))
metadata_alpha <- readRDS(here::here("Data",
                                   "metadata_alpha.rds"))

# Alpha diversity lakes and sp jitter
# Define parameters
lake_var <- "lake"
diversity_vars <- c("taxo_q0", "taxo_q1", "phylo_q0", "phylo_q1")
plot_titles <- c("Taxonomic (q0)", "Taxonomic (q1)", 
                 "Phylogenetic (q0)", "Phylogenetic (q1)")
boxplot_colors <- c("#608F3D","#41AEBD","#97E9D5","#F4DE3A","#A2CF49","#FCB11C")
jitter_colors <- c("#FFE800", "#A90C38")
jitter_labels <- c("PER" = "Perca fluviatilis", "LEP" = "Lepomis gibbosus")

# Function to generate the plots
# Function to generate lake alpha diversity plots
plot_alpha_lakes <- function(data, lake_var, diversity_vars, plot_titles, boxplot_colors, jitter_colors, jitter_labels) {
  
  # List to store individual plots
  plot_list <- list()
  
  # Loop through the diversity_vars to create plots
  for (i in seq_along(diversity_vars)) {
    p <- ggplot(data, aes_string(x = lake_var, y = diversity_vars[i])) +
      geom_boxplot(alpha = 0.8,
                   fill = boxplot_colors,
                   color = boxplot_colors) +
      geom_jitter(aes_string(colour = "origin"), position = position_jitter(0.07), cex = 2) +
      scale_color_manual(values = jitter_colors, labels = jitter_labels) +
      labs(color = "Species") +
      labs(x = NULL, y = NULL) + 
      theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")) + 
      ggtitle(plot_titles[i])
    
    plot_list[[i]] <- p
  }
  
  # Extract legend from the first plot
  legend <- get_legend(plot_list[[1]])
  
  # Combine the plots into a grid
  combined_plot <- plot_grid(
    plot_grid(plot_list[[1]] + theme(legend.position = "none"), 
              plot_list[[2]] + theme(legend.position = "none"),
              plot_list[[3]] + theme(legend.position = "none"),
              plot_list[[4]] + theme(legend.position = "none"), 
              ncol = 2),
    legend,
    ncol = 2, rel_widths = c(0.9, 0.15)
  )
  
  return(combined_plot)
}

# Generate the combined plot
alpha_plot_lake <- plot_alpha_lakes(metadata_alpha, lake_var, diversity_vars, plot_titles, boxplot_colors, jitter_colors, jitter_labels)

path_to_my_object = here::here("Figures","hill", "alpha_plot_lake.png")
ggsave(filename = path_to_my_object, plot = alpha_plot_lake, device = "png")


# alpha diversity species
# Define parameters for species alpha diversity plots
species_var <- "origin"
diversity_vars_species <- c("taxo_q0", "taxo_q1", "phylo_q0", "phylo_q1")
plot_titles_species <-c("Taxonomic (q0)", "Taxonomic (q1)", 
                        "Phylogenetic (q0)", "Phylogenetic (q1)")
fill_colors_species <- c("#FFE800", "#A90C38")

# Function to generate species alpha diversity plots
species_alpha_plots <- function(data, species_var, diversity_vars, plot_titles, fill_colors, border_color = "black") {
  
  # List to store individual plots
  plot_list <- list()
  
  # Loop through the diversity_vars to create plots
  for (i in seq_along(diversity_vars)) {
    p <- ggplot(data, aes_string(x = species_var, y = diversity_vars[i])) +
      geom_boxplot(alpha = 1,
                   fill = fill_colors,
                   color = border_color) +
      labs(x = NULL, y = NULL) + 
      theme(plot.title = element_text(hjust = 0, size = 12, face = "bold"),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")) + 
      ggtitle(plot_titles[i])
    
    plot_list[[i]] <- p
  }
  
  # Combine the plots into a grid
  combined_plot <- plot_grid(
    plot_grid(plot_list[[1]], plot_list[[2]], 
              plot_list[[3]], plot_list[[4]], 
              ncol = 2),
    ncol = 1
  )
  
  return(combined_plot)
}

alpha_plot_species <- species_alpha_plots(metadata_alpha, species_var, diversity_vars_species, plot_titles_species, fill_colors_species)

path_to_my_object = here::here("Figures","hill", "alpha_plot_species.png")
ggsave(filename = path_to_my_object, plot = alpha_plot_species, device = "png")


#Check the distribution and quantiles of alpha indices
indices_normality <- function(rich, nrow, ncol) {
  ### p-value < 0.05 means data failed normality test
  par(mfrow = c(nrow, ncol))
  
  for (i in names(rich)) {
    shap <- shapiro.test(rich[, i])
    qqnorm(rich[, i], main = i, sub = shap$p.value)
    qqline(rich[, i])
  }
  
  par(mfrow = c(1, 1))
}

metadata_alpha |>
  dplyr::select(taxo_q0,
                taxo_q1,
                phylo_q0,
                phylo_q1) |>
  indices_normality(nrow = 3, ncol = 2)

# Check homogeneity of variance between groups
# H0= equality of variances in the different populations
stats::bartlett.test(taxo_q0 ~ origin, metadata_alpha)
#Bartlett test of homogeneity of variances
#data:  taxo_q0 by origin
#Bartlett's K-squared = 11.34, df = 1, p-value = 0.0007587
stats::bartlett.test(taxo_q1 ~ origin, metadata_alpha)
#Bartlett test of homogeneity of variances
#data:  taxo_q1 by origin
#Bartlett's K-squared = 12.861, df = 1, p-value = 0.0003354
stats::bartlett.test(phylo_q0 ~ origin, metadata_alpha)
#Bartlett test of homogeneity of variances
#data:  phylo_q0 by origin
#Bartlett's K-squared = 0.46592, df = 1, p-value = 0.4949
stats::bartlett.test(phylo_q1 ~ origin, metadata_alpha)
#Bartlett test of homogeneity of variances
#data:  phylo_q1 by origin
#Bartlett's K-squared = 2.3783, df = 1, p-value = 0.123

# Custom function to create boxplot with pairwise comparison
create_boxplot_with_comparison <- function(data, x_var, y_var, title) {
  # Create boxplot
  graph <- ggplot(data, aes_string(x = x_var, y = y_var)) + 
    geom_boxplot(alpha = 1, fill = c("#FFE800", "#A90C38"), color = "black") +
    labs(x = NULL, y = NULL) + 
    geom_jitter(aes(colour = lake),
                position = position_jitter(0.02) ,
                cex=2.2)+
    scale_color_manual(values =c("#608F3D","#41AEBD","#97E9D5","#F4DE3A","#A2CF49","#FCB11C"))+
    labs(color = "Lakes") +
    theme(
      plot.title = element_text(hjust = 0, size = 12, face = "bold"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) + 
    ggtitle(title) +
    stat_summary(fun = mean, geom = "point", shape = 17, size = 3, color = "white")
  
  # Perform pairwise comparison
  pairwise_test <- compare_means(as.formula(paste(y_var, "~", x_var)), data, method = "wilcox.test")
  
  # Add p-value on graph
  graph + stat_pvalue_manual(
    pairwise_test,
    y.position = max(data[[y_var]]) * 1.1,
    label = "p.adj = {p.adj}",
    color = "blue",
    fontface = "bold",
    linetype = 1,
    tip.length = 0.01
  )
}


# Create graph
graph_taxo_q0 <- create_boxplot_with_comparison(metadata_alpha, "origin", "taxo_q0", "Taxonomic (q0)")
graph_taxo_q1 <- create_boxplot_with_comparison(metadata_alpha, "origin", "taxo_q1", "Taxonomic (q1)")
graph_phylo_q0 <- create_boxplot_with_comparison(metadata_alpha, "origin", "phylo_q0", "Phylogenetic (q0)")
graph_phylo_q1 <- create_boxplot_with_comparison(metadata_alpha, "origin", "phylo_q1", "Phylogenetic (q1)")

# Arrange all graphs in a grid on a single page with a common legend
common_legend <- cowplot::get_legend(graph_taxo_q0 + theme(legend.position = "right"))

# Combine graphs and legend
alpha_pvalue_sp <- grid.arrange(
  arrangeGrob(
    graph_taxo_q0 + theme(legend.position="none"),
    graph_taxo_q1 + theme(legend.position="none"),
    graph_phylo_q0 + theme(legend.position="none"),
    graph_phylo_q1 + theme(legend.position="none"),
    ncol = 2
  ),
  common_legend,
  ncol = 2,
  widths = c(1, 0.2))

path_to_my_object = here::here("Figures","hill", "alpha_pvalue_sp.png")
ggsave(filename = path_to_my_object, plot = alpha_pvalue_sp, device = "png")
