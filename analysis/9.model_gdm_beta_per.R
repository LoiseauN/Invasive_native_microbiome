#======== PROJECT COM2LIFE ========
## Modelisation GDM on beta diversity estimates in "estimates_beta_diversity.R" 
## First load the libraries and object. 
##### 4 models GDM by species, here we modelise 5 environmental varibales for Perca fluviatilis

# load libraries
library(gdm)

# Load objects
physeq_per <- readRDS(here::here("Data",
                                      "mon_objet_physeq_per.rds"))

metadata <- readRDS(here::here("Data",
                               "metadata.rds"))

phylo_q0_beta_matrix_per <- readRDS(here::here("Data",
                                 "matrix_phylo_q0_per.rds"))

phylo_q1_beta_matrix_per <- readRDS(here::here("Data",
                                           "matrix_phylo_q1_per.rds"))

tax_q0_beta_matrix_per <- readRDS(here::here("Data",
                                           "matrix_taxo_q0_per.rds"))

tax_q1_beta_matrix_per <- readRDS(here::here("Data",
                                         "matrix_taxo_q1_per.rds"))

### GDM
# -------------- EDIT DATA -------------------
#keep only metadata with PER samples
metadata <- metadata[grepl("PER", metadata$origin),]

#delete unwanted columns
metadata_filtered <- dplyr::select(metadata, -month, -origin, -lake_origin, -replicate,
                                   -month_full, -TPN, -TPC, -TPC.TPN, -NH4, -PO4, -NO3_N02, 
                                   -PROK_1_median, -PROK_2_median, -LARGE_CELLS_median, 
                                   -VLP_TOT_median,-CELLS_TOT_median, -Secchi)

# generate ID for site
metadata_filtered$site <- seq_along(rownames(metadata_filtered))

# add GPD coordinates for lakes
lake_coordinates <- data.frame(
  lake = c("CSM", "CERL", "CERS", "CRE", "VSS", "LGP"),  # Exemple de noms de lac avec leurs coordonnées
  latitude = c(48.86, 49.02, 49.03, 48.77, 48.99, 48.37),  # Exemple de latitudes correspondantes
  longitude = c(-2.597, -2.046, -2.05, -2.45, -1.968, -2.899)  # Exemple de longitudes correspondantes
)

# Save orginal rownames
original_rownames <- metadata_filtered$samples

# merge metadata with lake_coordinates by lake
metadata_filtered <- merge(metadata_filtered, lake_coordinates, by = "lake", all.x = TRUE, sort = FALSE)

# add rownames 
rownames(metadata_filtered) <- original_rownames

#delete samples and lake
metadata_filtered <- dplyr::select(metadata_filtered, -samples, -lake)

# add site column to distance matrix
site <- metadata_filtered$site

# Add the 'site' column to the distance matrix
tax_q0_beta_matrix <- cbind(site, tax_q0_beta_matrix_per)

# replace "-" by "."
rownames(metadata_filtered) <- gsub("\\-", ".", rownames(metadata_filtered))

# -------------- TAXO Q0 -------------------
# edit data for gdm function
gdmTab.dis <- formatsitepair(bioData= tax_q0_beta_matrix,
                             bioFormat=3, #diss matrix
                             XColumn="longitude",
                             YColumn="latitude",
                             predData= metadata_filtered,
                             siteColumn="site")

# apply gdm function
gdm.1 <- gdm(data=gdmTab.dis, geo=TRUE)

length(gdm.1$predictors) # get ideal of number of panels
plot(gdm.1, plot.layout=c(3,3))

# -------------- TAXO Q1 -------------------
# Add the 'site' column to the distance matrix
tax_q1_beta_matrix <- cbind(site, tax_q1_beta_matrix_per)

# edit data for gdm function
gdmTab.dis <- formatsitepair(bioData= tax_q1_beta_matrix,
                             bioFormat=3, #diss matrix
                             XColumn="longitude",
                             YColumn="latitude",
                             predData= metadata_filtered,
                             siteColumn="site")

# apply gdm function
gdm.2 <- gdm(data=gdmTab.dis, geo=TRUE)

length(gdm.2$predictors) # get ideal of number of panels
plot(gdm.2, plot.layout=c(3,3))


# -------------- PHYLO Q0 -------------------
# Add the 'site' column to the distance matrix
phylo_q0_beta_matrix <- cbind(site, phylo_q0_beta_matrix_per)

# edit for gdm function
gdmTab.dis <- formatsitepair(bioData= phylo_q0_beta_matrix,
                             bioFormat=3, #diss matrix
                             XColumn="longitude",
                             YColumn="latitude",
                             predData= metadata_filtered,
                             siteColumn="site")

# apply gdm function
gdm.3 <- gdm(data=gdmTab.dis, geo=TRUE)

length(gdm.3$predictors) # get ideal of number of panels
plot(gdm.3, plot.layout=c(3,3))


# -------------- PHYLO Q1 -------------------
# Add the 'site' column to the distance matrix
phylo_q1_beta_matrix <- cbind(site, phylo_q1_beta_matrix_per)

# edit for gdm function
gdmTab.dis <- formatsitepair(bioData= phylo_q1_beta_matrix,
                             bioFormat=3, #diss matrix
                             XColumn="longitude",
                             YColumn="latitude",
                             predData= metadata_filtered,
                             siteColumn="site")

# apply gdm function
gdm.4 <- gdm(data=gdmTab.dis, geo=TRUE)

length(gdm.4$predictors) # get ideal of number of panels
plot(gdm.4, plot.layout=c(3,3))

