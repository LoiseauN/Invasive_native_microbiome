#======== PROJECT COM2LIFE ========
## Statistical testing 

#======== PROJECT COM2LIFE ========
## Statistical testing: intra-species beta-diversity among regions

beta_diss <- readRDS(here::here("data", "beta_diss.rds"))

species_vec <- c("LEP", "PER")
metrics     <- c("taxo_q0", "taxo_q1", "phylo_q0", "phylo_q1")

# Build the value / region vectors for one species and one metric
get_intra <- function(sp, metric) {
  d <- beta_diss[beta_diss$reg_lev == "intra_region" &
                   (beta_diss$species_a == sp | beta_diss$species_b == sp), ]
  data.frame(value  = c(d[[metric]], d[[metric]]),
             region = factor(c(d$region_a, d$region_b)))
}

# ---- 1. Kruskal-Wallis for every species x metric ----
kw_tab <- expand.grid(species = species_vec, metric = metrics,
                      stringsAsFactors = FALSE)

kw_tab[, c("chi2", "df", "p_raw")] <- t(mapply(function(sp, m) {
  d <- get_intra(sp, m)
  k <- kruskal.test(value ~ region, data = d)
  c(k$statistic, k$parameter, k$p.value)
}, kw_tab$species, kw_tab$metric))

# ---- 2. FDR correction ACROSS all KW tests ----
kw_tab$p_fdr <- p.adjust(kw_tab$p_raw, method = "BH")
print(kw_tab)

# ---- 3. Post hoc Dunn only for KW tests passing FDR < 0.05 ----
sig <- kw_tab[kw_tab$p_fdr < 0.05, ]

posthoc <- lapply(seq_len(nrow(sig)), function(i) {
  d <- get_intra(sig$species[i], sig$metric[i])
  res <- dunn.test::dunn.test(d$value, d$region, method = "bonferroni",
                              kw = FALSE, table = FALSE)
  data.frame(species = sig$species[i], metric = sig$metric[i],
             comparison = res$comparisons, Z = res$Z,
             p_adj_bonf = res$P.adjusted)
})
posthoc <- do.call(rbind, posthoc)
print(posthoc)



#OLD VERSION








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

dunn_test_result <- dunn.test::dunn.test(taxo_q0_lep_intra, region_factor_intra, method = "bonferroni")
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

# Test Kruskal-Wallis to compare groups of region for PER 
kruskal_test_result_per_intra <- kruskal.test(taxo_q0_per_intra ~ region_factor_intra)
print(kruskal_test_result_per_intra)

dunn_test_result <- dunn.test::dunn.test(taxo_q0_per_intra, region_factor_intra, method = "bonferroni")
print(dunn_test_result)

## Taxo Q1
# Select taxo_q1 and intra region
taxo_q1_per_intra <- c(per_intra_region_data$taxo_q1, per_intra_region_data$taxo_q1)

# Test Kruskal-Wallis to compare groups of region for PER 
kruskal_test_result_per_intra <- kruskal.test(taxo_q1_per_intra ~ region_factor_intra)
print(kruskal_test_result_per_intra)

dunn_test_result <- dunn.test::dunn.test(taxo_q1_per_intra, region_factor_intra, method = "bonferroni")
print(dunn_test_result)

## Phylo_q0
# Select phylo_q0 and intra region
phylo_q0_per_intra <- c(per_intra_region_data$phylo_q0, per_intra_region_data$phylo_q0)

# Test Kruskal-Wallis to compare groups of region for PER 
kruskal_test_result_per_intra <- kruskal.test(phylo_q0_per_intra ~ region_factor_intra)
print(kruskal_test_result_per_intra)

dunn_test_result <- dunn.test::dunn.test(phylo_q0_per_intra, region_factor_intra, method = "bonferroni")
print(dunn_test_result)

## Phylo_q1
# Select phylo_q1 and intra region
phylo_q1_per_intra <- c(per_intra_region_data$phylo_q1, per_intra_region_data$phylo_q1)

# Test Kruskal-Wallis to compare groups of region for PER 
kruskal_test_result_per_intra <- kruskal.test(phylo_q1_per_intra ~ region_factor_intra)
print(kruskal_test_result_per_intra)

dunn_test_result <- dunn.test::dunn.test(phylo_q1_per_intra, region_factor_intra, method = "bonferroni")
print(dunn_test_result)



