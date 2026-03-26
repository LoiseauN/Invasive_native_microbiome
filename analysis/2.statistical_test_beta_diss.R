#======== PROJECT COM2LIFE ========
## Statistical testing 

# load beta_diss
beta_diss <- readRDS(here::here("data",
                                "beta_diss.rds"))
### Statistics for intra species

# ======= LEPOMIS GIBBOSUS ========
## Taxo Q0
# Filter out data to keep only lepomis gibbosus
lep_intra_region_data <- beta_diss[(beta_diss$species_a == "LEP" & beta_diss$reg_lev == "intra_region") |
                                     (beta_diss$species_b == "LEP" & beta_diss$reg_lev == "intra_region"), ]

# Select taxo_q0 and region
taxo_q0_lep_intra <- c(lep_intra_region_data$taxo_q0, lep_intra_region_data$taxo_q0)
region_lep_intra <- c(lep_intra_region_data$region_a, lep_intra_region_data$region_b)

# Transform region into factor
region_factor_intra <- factor(region_lep_intra)

# Test Kruskal-Wallis to compare groups of region for LEp 
kruskal_test_result_lep_intra <- kruskal.test(taxo_q0_lep_intra ~ region_factor_intra)
print(kruskal_test_result_lep_intra)

dunn_test_result <- dunn.test(taxo_q0_lep_intra, region_factor_intra, method = "bonferroni")
print(dunn_test_result)

## Taxo Q1
# Select taxo_q1 and intra region
taxo_q1_lep_intra <- c(lep_intra_region_data$taxo_q1, lep_intra_region_data$taxo_q1)

# Test Kruskal-Wallis to compare groups of region for LEp 
kruskal_test_result_lep_intra <- kruskal.test(taxo_q1_lep_intra ~ region_factor_intra)
print(kruskal_test_result_lep_intra)

dunn_test_result <- dunn.test::dunn.test(taxo_q1_lep_intra, region_factor_intra, method = "bonferroni")
print(dunn_test_result)

## Phylo_q0
# Select phylo_q0 and intra region
phylo_q0_lep_intra <- c(lep_intra_region_data$phylo_q0, lep_intra_region_data$phylo_q0)

# Test Kruskal-Wallis to compare groups of region for LEp 
kruskal_test_result_lep_intra <- kruskal.test(phylo_q0_lep_intra ~ region_factor_intra)
print(kruskal_test_result_lep_intra)

dunn_test_result <- dunn.test::dunn.test(phylo_q0_lep_intra, region_factor_intra, method = "bonferroni")
print(dunn_test_result)

## Phylo_q1
# Select phylo_q1 and intra region
phylo_q1_lep_intra <- c(lep_intra_region_data$phylo_q1, lep_intra_region_data$phylo_q1)

# Test Kruskal-Wallis to compare groups of region for LEp 
kruskal_test_result_lep_intra <- kruskal.test(phylo_q1_lep_intra ~ region_factor_intra)
print(kruskal_test_result_lep_intra)

dunn_test_result <- dunn.test::dunn.test(phylo_q1_lep_intra, region_factor_intra, method = "bonferroni")
print(dunn_test_result)


# ========= PERCA FLUVIATILIS ========
## Taxo Q0
# Filter out data to keep only Perca fluviatilis
per_intra_region_data <- beta_diss[(beta_diss$species_a == "PER" & beta_diss$reg_lev == "intra_region") |
                                     (beta_diss$species_b == "PER" & beta_diss$reg_lev == "intra_region"), ]

# Select taxo_q0 and intra region
taxo_q0_per_intra <- c(per_intra_region_data$taxo_q0, per_intra_region_data$taxo_q0)
region_per_intra <- c(per_intra_region_data$region_a, per_intra_region_data$region_b)

# Trasnform region into factor
region_factor_intra <- factor(region_per_intra)

# Test Kruskal-Wallis to compare groups of region for LEp 
kruskal_test_result_per_intra <- kruskal.test(taxo_q0_per_intra ~ region_factor_intra)
print(kruskal_test_result_per_intra)

dunn_test_result <- dunn.test::dunn.test(taxo_q0_per_intra, region_factor_intra, method = "bonferroni")
print(dunn_test_result)

## Taxo Q1
# Select taxo_q1 and intra region
taxo_q1_per_intra <- c(per_intra_region_data$taxo_q1, per_intra_region_data$taxo_q1)

# Test Kruskal-Wallis to compare groups of region for LEp 
kruskal_test_result_per_intra <- kruskal.test(taxo_q1_per_intra ~ region_factor_intra)
print(kruskal_test_result_per_intra)

dunn_test_result <- dunn.test::dunn.test(taxo_q1_per_intra, region_factor_intra, method = "bonferroni")
print(dunn_test_result)

## Phylo_q0
# Select phylo_q0 and intra region
phylo_q0_per_intra <- c(per_intra_region_data$phylo_q0, per_intra_region_data$phylo_q0)

# Test Kruskal-Wallis to compare groups of region for LEp 
kruskal_test_result_per_intra <- kruskal.test(phylo_q0_per_intra ~ region_factor_intra)
print(kruskal_test_result_per_intra)

dunn_test_result <- dunn.test::dunn.test(phylo_q0_per_intra, region_factor_intra, method = "bonferroni")
print(dunn_test_result)

## Phylo_q1
# Select phylo_q1 and intra region
phylo_q1_per_intra <- c(per_intra_region_data$phylo_q1, per_intra_region_data$phylo_q1)

# Test Kruskal-Wallis to compare groups of region for LEp 
kruskal_test_result_per_intra <- kruskal.test(phylo_q1_per_intra ~ region_factor_intra)
print(kruskal_test_result_per_intra)

dunn_test_result <- dunn.test::dunn.test(phylo_q1_per_intra, region_factor_intra, method = "bonferroni")
print(dunn_test_result)