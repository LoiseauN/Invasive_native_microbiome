#======== PROJECT COM2LIFE ========
## This script is a continuation of the reads_percentage_species script 
## plot the ASVs shared between water and species based on chlorohyll a level of each make

# import data
data_import_2 <- function(name) {
  file_path <- here::here(paste("Data/", name, sep = "")) 
  data <- read.csv(file_path, sep = ";", header = TRUE)
}
param_table_july_2021 <- data_import_2("param_table_july2021.csv")

# create df
total_lakes_per <- tibble::tibble(
  lake = c("CRJ1", "CRJ2", "CHA", "CTL", "GDP", "VER"),
  percentage_total_shared = c(47.12, 50.77, 55.46, 19.3, 66.36, 35.14),
  species_type = "Perca fluviatilis"
)

total_lakes_lep <- tibble::tibble(
  lake = c("CRJ1", "CRJ2", "CHA", "CTL", "GDP", "VER"),
  percentage_total_shared = c(35.75, 37.57, 45.04, 11.88, 42.94, 27.92),
  species_type = "Lepomis gibbosus"
)

# combined two df
combined_df <- dplyr::bind_rows(total_lakes_per, total_lakes_lep)

#Merge dataframes on the "lake" column
combined_df <- merge(combined_df, param_table_july_2021[, c("lake", "Chla_median")], by="lake", all.x=TRUE)

# replace with the right name of lake 
combined_df <- combined_df %>%
  mutate(lake = case_when(
    lake == 'CHA' ~ 'CSM',
    lake == 'CRJ1' ~ 'CERL',
    lake == 'CRJ2' ~ 'CERS',
    lake == 'VER' ~ 'VSS',
    lake == 'CTL' ~ 'CRE',
    lake == 'GDP' ~ 'LGP'
  ))
# plot
percentage_shared_chla <- ggplot2::ggplot(combined_df, ggplot2::aes(x = Chla_median, y = percentage_total_shared, color = species_type, label = lake)) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_text(hjust = 0.5, vjust = -0.5) +  
  ggplot2::labs(
       x = "Concentration of chlorophyll a (µg.l-1)",
       y = "Total percentage of reads corresponding to ASVs shared with water",
       color = "Species") +
  ggplot2::scale_color_manual(values = c("Lepomis gibbosus" = "darkorange", "Perca fluviatilis" =  "#A90C38"),
                              labels = c(expression(italic("Lepomis gibbosus")), expression(italic("Perca fluviatilis")))) +  # Définir les couleurs pour chaque type d'espèce
  ggplot2::theme_minimal() +
  theme(legend.text = element_text(size = 10.5))


path_to_my_object = here::here("Figures", "percentage_shared_chla.png")
ggplot2::ggsave(filename = path_to_my_object, plot = percentage_shared_chla, device = "png")
