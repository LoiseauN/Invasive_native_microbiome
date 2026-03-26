#======== PROJECT COM2LIFE ========
## Graphical representation of GDM models
## Pourcentage of predictors for each model and each species

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyverse)

# ==== PERCA FLUVIATILIS =========
# Data for gdm.1
gdm.1 <- data.frame(
  model = "taxo_q0",
  predictor = c("Chla", "O2", "Salinity", "Temperature", "Spatial", "pH"),
  coefficient = c(0.881, 0.53, 0.378, 0.251, 0.22, 0.074),
  stringsAsFactors = FALSE
)

# Data for gdm.2
gdm.2 <- data.frame(
  model = "taxo_q1",
  predictor = c("Chla", "O2", "Temperature", "Spatial", "Salinity", "pH"),
  coefficient = c(1.091, 0.584, 0.566, 0.484, 0.481, 0),
  stringsAsFactors = FALSE
)

# Data for gdm.3
gdm.3 <- data.frame(
  model = "phylo_q0",
  predictor = c("Spatial", "Chla", "O2", "Salinity", "Temperature", "pH"),
  coefficient = c(0.337, 0.334, 0.081, 0.079, 0.03, 0),
  stringsAsFactors = FALSE
)

# Data for gdm.4
gdm.4 <- data.frame(
  model = "phylo_q1",
  predictor = c("Chla", "Spatial", "O2", "Temperature", "Salinity", "pH"),
  coefficient = c(0.156, 0.133, 0.111, 0.048, 0, 0),
  stringsAsFactors = FALSE
)

# Concatenated dataframe
combined_data <- rbind(gdm.1, gdm.2, gdm.3, gdm.4)

# Calculated each porcentage
combined_data <- combined_data %>%
  group_by(model) %>%
  mutate(total_coefficient = sum(coefficient)) %>%
  mutate(percentage = (coefficient / total_coefficient) * 100) %>%
  ungroup()

# Convert model into factor
combined_data$model <- factor(combined_data$model, levels = c("taxo_q0", "taxo_q1", "phylo_q0", "phylo_q1"))

# create barplot
results_gdm_per <- ggplot2::ggplot(combined_data, ggplot2::aes(x = model, y = percentage, fill = predictor)) +
  ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.7) +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Models", y = "Percentage (%)", fill = "Predictors") +
  ggplot2::scale_fill_brewer(palette = "Set3")

path_to_my_object = here::here("figures","figure4A.png")
ggplot2::ggsave(filename = path_to_my_object, plot = results_gdm_per, device = "png")

#===== LEPOMIS GIBBOSUS =======
# Data for gdm.1
gdm.1 <- data.frame(
  model = "taxo_q0",
  predictor = c("Chla", "Salinity", "pH", "O2", "Temperature", "Spatial"),
  coefficient = c(0.618, 0.447, 0.222, 0.191, 0.062,0),
  stringsAsFactors = FALSE
)

# Data for gdm.2
gdm.2 <- data.frame(
  model = "taxo_q1",
  predictor = c("Chla", "Salinity", "O2", "pH", "Temperature", "Spatial"),
  coefficient = c(1.261, 0.522, 0.08, 0.047, 0.041, 0),
  stringsAsFactors = FALSE
)

# Data for gdm.3
gdm.3 <- data.frame(
  model = "phylo_q0",
  predictor = c("pH", "Salinity", "Temperature", "Spatial", "O2", "Chla"),
  coefficient = c(0.467, 0.177, 0.07, 0.054, 0, 0),
  stringsAsFactors = FALSE
)

# Data for gdm.4
gdm.4 <- data.frame(
  model = "phylo_q1",
  predictor = c("Chla", "Salinity", "Spatial", "Temperature", "pH", "O2"),
  coefficient = c(0.113, 0.091, 0.006, 0.014, 0.011, 0.008),
  stringsAsFactors = FALSE
)

# Concatenated dataframe
combined_data <- rbind(gdm.1, gdm.2, gdm.3, gdm.4)

# Calculated each porcentage
combined_data <- combined_data %>%
  group_by(model) %>%
  mutate(total_coefficient = sum(coefficient)) %>%
  mutate(percentage = (coefficient / total_coefficient) * 100) %>%
  ungroup()

# Convert model into factor
combined_data$model <- factor(combined_data$model, levels = c("taxo_q0", "taxo_q1", "phylo_q0", "phylo_q1"))

# Plot
results_gdm_lep <- ggplot2::ggplot(combined_data, ggplot2::aes(x = model, y = percentage, fill = predictor)) +
  ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.7) +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Models", y = "Percentage (%)", fill = "Predictors") +
  ggplot2::scale_fill_brewer(palette = "Set3")

# Save
path_to_my_object = here::here("figures","figure4B.png")
ggplot2::ggsave(filename = path_to_my_object, plot = results_gdm_lep, device = "png")

# Create plot
plot_left <- results_gdm_per +
  ggplot2::theme(legend.position = "none")  

plot_right <- results_gdm_lep +
  ggplot2::theme(axis.title.y = ggplot2::element_blank())  

combined_plot <- cowplot::plot_grid(plot_left, plot_right, labels = c("A", "B"), ncol = 2, align = "hv")

# Add legend
combined_plot_gdm <- cowplot::plot_grid(
  combined_plot,
  cowplot::draw_plot_label(label = "Legend", size = 15, hjust = 0, vjust = 1),  # Positionner la légende à droite
  ncol = 2,
  rel_widths = c(1, 0.1),  
  align = "v"  
)

# Save
path_to_my_object = here::here("figures","figure4.png")
ggplot2::ggsave(filename = path_to_my_object, plot = combined_plot_gdm, device = "png")
