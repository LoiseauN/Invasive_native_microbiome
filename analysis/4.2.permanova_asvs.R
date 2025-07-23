#======== PROJECT COM2LIFE ========
## permanova on beta diversity matrix to see if there is lake effect or species effect

#library
library(dplyr)
library(spaa)
library(vegan)

# load data 
tax_q0_matrix_per <- readRDS(here::here("Data",
                                  "matrix_taxo_q0_per.rds"))
tax_q0_matrix_lep <- readRDS(here::here("Data",
                                        "matrix_taxo_q0_lep.rds"))
tax_q1_matrix_per <- readRDS(here::here("Data",
                                  "matrix_taxo_q1_per.rds"))
tax_q1_matrix_lep <- readRDS(here::here("Data",
                                        "matrix_taxo_q1_lep.rds"))
phylo_q0_matrix_per <- readRDS(here::here("Data",
                                        "matrix_phylo_q0_per.rds"))
phylo_q0_matrix_lep <- readRDS(here::here("Data",
                                          "matrix_phylo_q0_lep.rds"))
phylo_q1_matrix_per <- readRDS(here::here("Data",
                                          "matrix_phylo_q1_per.rds"))
phylo_q1_matrix_lep <- readRDS(here::here("Data",
                                          "matrix_phylo_q1_lep.rds"))
metadata <- readRDS(here::here("Data",
                               "metadata.rds"))

# ===== PERMANOVA - LAKE EFFECT - PERCA FLUVIATILIS =====
metadata_per <- metadata[grep("PER", metadata$origin), ]


permanova_tax_q0_per <- vegan::adonis2(tax_q0_matrix_per ~ lake,
                                        data = metadata_per,
                                        perm = 1000)
permanova_tax_q0_per

permanova_tax_q1_per <- vegan::adonis2(tax_q1_matrix_per ~ lake,
                                        data = metadata_per,
                                        perm = 1000)
permanova_tax_q1_per

permanova_phylo_q0_per <- vegan::adonis2(phylo_q0_matrix_per ~ lake,
                                          data = metadata_per,
                                          perm = 1000)
permanova_phylo_q0_per

permanova_phylo_q1_per <- vegan::adonis2(phylo_q1_matrix_per ~ lake,
                                          data = metadata_per,
                                          perm = 1000)
permanova_phylo_q1_per

## Betadisp()
# === TAXO Q0 ===
tax_q0_dist_per <- dist(tax_q0_matrix_per)

disp <- betadisper(tax_q0_dist_per, metadata_per$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

# === TAXO Q1 ===
tax_q1_dist_per <- dist(tax_q1_matrix_per)

disp <- betadisper(tax_q1_dist_per, metadata_per$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

# === PHYLO Q0 ===
phylo_q0_dist_per <- dist(phylo_q0_matrix_per)

disp <- betadisper(phylo_q0_dist_per, metadata_per$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

# === PHYLO Q1 ===
phylo_q1_dist_per <- dist(phylo_q1_matrix_per)

disp <- betadisper(phylo_q1_dist_per, metadata_per$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

# ===== PERMANOVA - LAKE EFFECT - LEPOMIS GIBBOSUS =====
metadata_lep <- metadata[grep("LEP", metadata$origin), ]


lepmanova_tax_q0_lep <- vegan::adonis2(tax_q0_matrix_lep ~ lake,
                                       data = metadata_lep,
                                       lepm = 1000)
lepmanova_tax_q0_lep

lepmanova_tax_q1_lep <- vegan::adonis2(tax_q1_matrix_lep ~ lake,
                                       data = metadata_lep,
                                       lepm = 1000)
lepmanova_tax_q1_lep

lepmanova_phylo_q0_lep <- vegan::adonis2(phylo_q0_matrix_lep ~ lake,
                                         data = metadata_lep,
                                         lepm = 1000)
lepmanova_phylo_q0_lep

lepmanova_phylo_q1_lep <- vegan::adonis2(phylo_q1_matrix_lep ~ lake,
                                         data = metadata_lep,
                                         lepm = 1000)
lepmanova_phylo_q1_lep

## Betadisp()
# === TAXO Q0 ===
tax_q0_dist_lep <- dist(tax_q0_matrix_lep)

disp <- betadisper(tax_q0_dist_lep, metadata_lep$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

# === TAXO Q1 ===
tax_q1_dist_lep <- dist(tax_q1_matrix_lep)

disp <- betadisper(tax_q1_dist_lep, metadata_lep$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

# === PHYLO Q0 ===
phylo_q0_dist_lep <- dist(phylo_q0_matrix_lep)

disp <- betadisper(phylo_q0_dist_lep, metadata_lep$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

# === PHYLO Q1 ===
phylo_q1_dist_lep <- dist(phylo_q1_matrix_lep)

disp <- betadisper(phylo_q1_dist_lep, metadata_lep$lake )
print(disp)

# Significativity of dispersion
perm_test <- permutest(disp, permutations = 999)
print(perm_test)

# Visualize
plot(disp)

## No parwise adonis here, as lake variance homogeneity shows no significant differences