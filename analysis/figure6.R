#==== PROJECT COM2LIFE=====
# perform ecdf 

#load data
physeq_per <- readRDS(here::here("Data",
                                 "mon_objet_physeq_per.rds"))
physeq_lep <- readRDS(here::here("Data",
                                 "mon_objet_physeq_lep.rds"))
sample_ID_conv <- readRDS(here::here("Data",
                                     "sample_ID_conv.rds"))
all_metabo_filtered <- readRDS(here::here("Data",
                                          "all_metabo_filtered.rds"))
metadata <- readRDS(here::here("Data",
                               "metadata.rds"))
# otu_table
otu_table_per <- as.data.frame(physeq_per@otu_table)
otu_table_lep <- as.data.frame(physeq_lep@otu_table)

# function to rename samples by rownames
rename_metabolomic_samples <- function(all_metabo_filtered, sample_id_conv, id_metabolomics_col, id_sequencing_col) {
  all_metabo_filtered <- as.data.frame(all_metabo_filtered)
  sample_ids <- rownames(all_metabo_filtered)
  sample_id_filtered <- sample_id_conv %>%
    filter(!!sym(id_metabolomics_col) %in% sample_ids)
  rename_mapping <- setNames(sample_id_filtered[[id_sequencing_col]], sample_id_filtered[[id_metabolomics_col]])
  rownames(all_metabo_filtered) <- ifelse(sample_ids %in% names(rename_mapping), rename_mapping[sample_ids], sample_ids)
  return(all_metabo_filtered)
}

# rename samples
all_metabo_filtered<- rename_metabolomic_samples(
  all_metabo_filtered = all_metabo_filtered,
  sample_id_conv = sample_ID_conv,
  id_metabolomics_col = "ID_Metabolomics",
  id_sequencing_col = "ID_Sequencing"
)

rownames(metadata) <- metadata$samples
colnames(metadata)[colnames(metadata) == "lake"] <- "Lake"

# Perform the join directly by creating temporary 'rowname' columns
all_metabo_filtered <- all_metabo_filtered %>%
  tibble::rownames_to_column(var = "rowname") %>%
  left_join(metadata %>% tibble::rownames_to_column(var = "rowname") %>% dplyr::select(rowname, Chla_median), 
            by = "rowname") %>%
  tibble::column_to_rownames(var = "rowname")

#combine df
lepomis_data <- all_metabo_filtered %>%
  filter(Species == "Lepomis") %>%
  mutate(Species = "Lepomis gibbosus")

perca_data <- all_metabo_filtered %>%
  filter(Species == "Perca") %>%
  mutate(Species = "Perca fluviatilis")

# Combine data
combined_data <- bind_rows(lepomis_data, perca_data)

# ecdf
ecdf <- ggplot2::ggplot(combined_data, ggplot2::aes(x = Chla_median, color = Species)) +
  ggplot2::stat_ecdf(geom = "step") +
  ggplot2::scale_color_manual(values = c("Lepomis gibbosus" = "darkorange", "Perca fluviatilis" = "#A90C38")) +
  ggplot2::labs(x = "BMD",
                y = "ECDF",
                color = "Species") +
  ggplot2::theme_minimal() +
  ggplot2::theme(  
    axis.ticks = ggplot2::element_line(size = 0.5),  # Taille des ticks
    axis.line = ggplot2::element_line(size = 0.5),
    legend.text = ggplot2::element_text(face = "italic", size = 11),
    axis.text.x = ggplot2::element_text(size = 10),  # Agrandir la taille des graduations de l'axe x
    axis.text.y = ggplot2::element_text(size = 10)# Taille de la ligne de l'axe
  )

path_to_my_object = here::here("figures", "figure6.png")
ggplot2::ggsave(filename = path_to_my_object, plot = ecdf, device = "png")
