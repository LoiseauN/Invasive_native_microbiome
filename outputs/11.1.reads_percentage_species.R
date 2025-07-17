#======== PROJECT COM2LIFE ========
## This script is a continuation of the ASVs_by_lakes script 
## plot the ASVs shared between water and species and ASVs unique

#load physeq object
physeq_filtered <- readRDS(here::here("data",
                                      "mon_objet_physeq_filtered_all.rds"))
asv_table_filtered <- as.data.frame(physeq_filtered@otu_table)

# Define patterns for each site
sites <- list(
  cerl = "CERL",
  cers = "CERS",
  cre = "CRE",
  vss = "VSS",
  lgp = "LGP",
  csm = "CSM"
)

# Function to filter ASVs based on a given pattern
filter_asv <- function(asv_table, pattern) {
  asv_table <- asv_table[grepl(pattern, rownames(asv_table)), ]
  return(asv_table)
}

# Function to filter columns with a non-zero total of reads
filter_non_zero_columns <- function(asv_table) {
  col_totals <- colSums(asv_table)
  non_zero_cols <- col_totals != 0
  asv_table <- asv_table[, non_zero_cols]
  return(asv_table)
}

# Function to calculate percentages
calculate_percentages <- function(total_shared, total_unique) {
  total <- total_shared + total_unique
  percentage_shared <- (total_shared / total) * 100
  percentage_unique <- (total_unique / total) * 100
  return(c(total, total_shared, total_unique, percentage_shared, percentage_unique))
}

# Function to create a histogram
create_histogram <- function(pourcentages, title) {
  noms_colonnes <- c("All reads shared", "Reads unique")
  bar_centers <- barplot(pourcentages,
                         names.arg = noms_colonnes,
                         col = c("#AE123A", "#440154ff"),
                         main = title,
                         ylab = "Percentage of reads",
                         ylim = c(0, 100),
                         border = FALSE)
  text(x = bar_centers,
       y = pourcentages,
       labels = paste0(round(pourcentages, 2), "%"),
       pos = 3,
       cex = 0.8,
       col = "black",
       font = 2)
}

# ======= PERCA FLUVIATILIS =======
# Function to obtain ASVs shared by all (specific criteria)
get_shared_asv <- function(asv_table) {
  asv_table_shared <- asv_table %>%
    purrr::keep(~ sum(. > 0) != 1) %>%
    .[, -which(apply(., 2, function(col) {
      per_values <- col[grepl("PER", rownames(asv_table))]
      adn_values <- col[grepl("ADN", rownames(asv_table))]
      return(all(per_values == 0) | all(adn_values == 0))
    }))]
  return(asv_table_shared)
}

# Iterate over each site and apply transformations
for (site in names(sites)) {
  # Filter ASVs for the site
  asv_table <- filter_asv(asv_table_filtered, sites[[site]])
  asv_table <- asv_table[!grepl(paste0(site, "2"), rownames(asv_table)), ]
  
  # Filter ASVs for PER and ADN
  asv_table_per <- asv_table[grep("PER|ADN1", rownames(asv_table)), ]
  
  # Filter columns with non-zero totals
  asv_table_per <- filter_non_zero_columns(asv_table_per)
  
  # Create dataframes for shared and unique ASVs
  asv_shared <- get_shared_asv(asv_table_per)
  total_shared <- sum(as.matrix(asv_shared))
  
  # Check if columns contain positive values for PER and ADN rows
  contains_positive <- function(col, pattern) {
    any(col[grepl(pattern, rownames(asv_unique_per))] > 0)
  }
  
  columns_to_remove <- sapply(asv_unique_per, function(col) {
    contains_positive(col, "PER") & contains_positive(col, "ADN")
  })
  asv_unique_per <- asv_unique_per[, !columns_to_remove]
  
  total_unique <- sum(as.matrix(asv_unique_per))
  
  # Calculate percentages
  pourcentages <- calculate_percentages(total_shared, total_unique)[4:5]
  
  # Create a histogram
  create_histogram(pourcentages, paste("Percentage of reads of ASVs for", site))
}


# ======= LEPOMIS GIBBOSUS =======
# Function to obtain ASVs shared by all (specific criteria)
get_shared_asv <- function(asv_table) {
  asv_table_shared <- asv_table %>%
    purrr::keep(~ sum(. > 0) != 1) %>%
    .[, -which(apply(., 2, function(col) {
      lep_values <- col[grepl("GAR", rownames(asv_table))]
      adn_values <- col[grepl("ADN", rownames(asv_table))]
      return(all(lep_values == 0) | all(adn_values == 0))
    }))]
  return(asv_table_shared)
}

# Iterate over each site and apply transformations
for (site in names(sites)) {
  # Filter ASVs for the site
  asv_table <- filter_asv(asv_table_filtered, sites[[site]])
  asv_table <- asv_table[!grepl(paste0(site, "2"), rownames(asv_table)), ]
  
  # Filter ASVs for PER and ADN
  asv_table_lep <- asv_table[grep("GAR|ADN1", rownames(asv_table)), ]
  
  # Filter columns with non-zero totals
  asv_table_lep <- filter_non_zero_columns(asv_table_lep)
  
  # Create dataframes for shared and unique ASVs
  asv_shared <- get_shared_asv(asv_table_lep)
  total_shared <- sum(as.matrix(asv_shared))
  
  # Calculate unique ASVs
  water_asv <- get(paste0("water_asv_", site))
  lep_asv <- get(paste0("lep_asv_", site))
  water_diff <- setdiff(water_asv, lep_asv)
  lep_diff <- setdiff(lep_asv, water_asv)
  asv_unique <- union(lep_diff, water_diff)
  
  asv_unique_lep <- asv_table_lep[, colnames(asv_table_lep) %in% asv_unique]
  
  # Check if columns contain positive values for PER and ADN rows
  contains_positive <- function(col, pattern) {
    any(col[grepl(pattern, rownames(asv_unique_lep))] > 0)
  }
  
  columns_to_remove <- sapply(asv_unique_lep, function(col) {
    contains_positive(col, "GAR") & contains_positive(col, "ADN")
  })
  asv_unique_lep <- asv_unique_lep[, !columns_to_remove]
  
  total_unique <- sum(as.matrix(asv_unique_lep))
  
  # Calculate percentages
  pourcentages <- calculate_percentages(total_shared, total_unique)[4:5]
  
  # Create a histogram
  create_histogram(pourcentages, paste("Percentage of reads of ASVs for", site))
}

