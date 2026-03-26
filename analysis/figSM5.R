# === Data preparation ===

#==== PERCA FLUVIATILIS =====
# === Data preparation ===

# Convert sample names to vectors
otu_samples <- rownames(otu_table_per)
metabo_samples <- rownames(annotated_metabo_per)

# Find common samples
common_samples <- intersect(otu_samples, metabo_samples)

# Filter tables to keep only common samples
otu_table_filtered <- otu_table_per[common_samples, ]
metabo_table_filtered <- annotated_metabo_per[common_samples, ]

# Check if row names match
stopifnot(all(rownames(otu_table_filtered) == rownames(metabo_table_filtered)))

# Remove columns with zero standard deviation
otu_table_filtered <- otu_table_filtered[, apply(otu_table_filtered, 2, sd) != 0]
metabo_table_filtered <- metabo_table_filtered[, apply(metabo_table_filtered, 2, sd) != 0]

# === Taxonomic labels for OTUs ===
tax_per <- as.data.frame(physeq_per@tax_table)
tax_per$Label <- ifelse(is.na(tax_per$Genus), tax_per$Family, tax_per$Genus)

# Replace column names in filtered OTU table by taxonomic labels
colnames(otu_table_filtered) <- tax_per[colnames(otu_table_filtered), "Label"]

# === Required packages ===
library(igraph)
library(dplyr)
library(paletteer)

# === Calculate correlations between OTUs and metabolites ===
cor_matrix <- cor(otu_table_filtered, metabo_table_filtered)

# Convert to long-format dataframe for visualization
cor_df <- as.data.frame(as.table(cor_matrix))
colnames(cor_df) <- c("OTU", "Metabolite", "Correlation")

# Filter correlations above threshold (absolute value)
threshold <- 0.77
cor_df_filtered <- cor_df %>% filter(abs(Correlation) > threshold)

# Convert to character and remove NA rows
cor_df_filtered$OTU <- as.character(cor_df_filtered$OTU)
cor_df_filtered$Metabolite <- as.character(cor_df_filtered$Metabolite)
cor_df_filtered <- na.omit(cor_df_filtered)

cat("Number of edges after filtering:", nrow(cor_df_filtered), "\n")

# === Create bipartite graph ===
bipartite_edges <- cor_df_filtered[, c("OTU", "Metabolite", "Correlation")]
g <- graph_from_data_frame(bipartite_edges, directed = FALSE)

# Identify node types (TRUE = OTU, FALSE = Metabolite)
V(g)$type <- bipartite_mapping(g)$type

# === Node colors ===
# Metabolite colors from paletteer diverging palette
metabolite_colors <- paletteer_c("ggthemes::Red-Blue-White Diverging", 30)
# Grayscale gradient for OTUs
asv_colors <- gray.colors(sum(V(g)$type), start = 0.9, end = 0.3)

# Assign colors based on node type
vertex_colors <- ifelse(V(g)$type, asv_colors, metabolite_colors)

# === Node sizes ===
vertex_sizes <- ifelse(V(g)$type, 5, 7)

# === Bipartite layout ===
layout <- matrix(NA, nrow = vcount(g), ncol = 2)
layout[V(g)$type == TRUE, 1] <- 1   # OTUs on the left
layout[V(g)$type == FALSE, 1] <- 2  # Metabolites on the right

layout[V(g)$type == TRUE, 2] <- seq(-1, 1, length.out = sum(V(g)$type == TRUE))
layout[V(g)$type == FALSE, 2] <- seq(-1, 1, length.out = sum(V(g)$type == FALSE))

# === Label positioning ===
vertex_label_dist <- ifelse(V(g)$type, 8, 1.5)  # Closer labels for metabolites
vertex_label_degree <- ifelse(V(g)$type, pi, 0) # OTUs labels left, metabolites right

# === Edge colors ===
edge_colors <- rep("#006400", ecount(g))  # Dark green

# === Final plot ===
plot(g, 
     layout = layout,
     vertex.label.cex = 0.7,
     edge.width = abs(E(g)$Correlation) * 2,
     edge.color = edge_colors,
     vertex.color = vertex_colors,
     vertex.size = vertex_sizes,
     vertex.label.color = "black",
     vertex.label.dist = vertex_label_dist,
     vertex.label.degree = vertex_label_degree,
     main = "Bipartite Network Analysis of OTUs and Metabolites")


#==== LEPOMIS GIBBOSUS =====
# Convert sample names to vectors
otu_samples <- rownames(otu_table_lep)
metabo_samples <- rownames(annotated_metabo_lep)

# Find common samples
common_samples <- intersect(otu_samples, metabo_samples)

# Filter tables to keep only common samples
otu_table_filtered <- otu_table_lep[common_samples, ]
metabo_table_filtered <- annotated_metabo_lep[common_samples, ]

# Check if row names match
stopifnot(all(rownames(otu_table_filtered) == rownames(metabo_table_filtered)))

# Remove columns with zero standard deviation
otu_table_filtered <- otu_table_filtered[, apply(otu_table_filtered, 2, sd) != 0]
metabo_table_filtered <- metabo_table_filtered[, apply(metabo_table_filtered, 2, sd) != 0]

# === Taxonomic labels for OTUs ===
tax_lep <- as.data.frame(physeq_lep@tax_table)
tax_lep$Label <- ifelse(is.na(tax_lep$Genus), tax_lep$Family, tax_lep$Genus)

# Replace column names in filtered OTU table by taxonomic labels
colnames(otu_table_filtered) <- tax_lep[colnames(otu_table_filtered), "Label"]

# === Required packages ===
library(igraph)
library(dplyr)
library(paletteer)

# === Calculate correlations between OTUs and metabolites ===
cor_matrix <- cor(otu_table_filtered, metabo_table_filtered)

# Convert to long-format dataframe for visualization
cor_df <- as.data.frame(as.table(cor_matrix))
colnames(cor_df) <- c("OTU", "Metabolite", "Correlation")

# Filter correlations above threshold (absolute value)
threshold <- 0.77
cor_df_filtered <- cor_df %>% filter(abs(Correlation) > threshold)

# Convert to character and remove NA rows
cor_df_filtered$OTU <- as.character(cor_df_filtered$OTU)
cor_df_filtered$Metabolite <- as.character(cor_df_filtered$Metabolite)
cor_df_filtered <- na.omit(cor_df_filtered)

cat("Number of edges after filtering:", nrow(cor_df_filtered), "\n")

# === Create bipartite graph ===
bipartite_edges <- cor_df_filtered[, c("OTU", "Metabolite", "Correlation")]
g <- graph_from_data_frame(bipartite_edges, directed = FALSE)

# Identify node types (TRUE = OTU, FALSE = Metabolite)
V(g)$type <- bipartite_mapping(g)$type

# === Node colors ===
# Metabolite colors from paletteer diverging palette
metabolite_colors <- paletteer_c("ggthemes::Red-Blue-White Diverging", 30)
# Grayscale gradient for OTUs
asv_colors <- gray.colors(sum(V(g)$type), start = 0.9, end = 0.3)

# Assign colors based on node type
vertex_colors <- ifelse(V(g)$type, asv_colors, metabolite_colors)

# === Node sizes ===
vertex_sizes <- ifelse(V(g)$type, 5, 7)

# === Bipartite layout ===
layout <- matrix(NA, nrow = vcount(g), ncol = 2)
layout[V(g)$type == TRUE, 1] <- 1   # OTUs on the left
layout[V(g)$type == FALSE, 1] <- 2  # Metabolites on the right

layout[V(g)$type == TRUE, 2] <- seq(-1, 1, length.out = sum(V(g)$type == TRUE))
layout[V(g)$type == FALSE, 2] <- seq(-1, 1, length.out = sum(V(g)$type == FALSE))

# === Label positioning ===
vertex_label_dist <- ifelse(V(g)$type, 8, 1.5)  # Closer labels for metabolites
vertex_label_degree <- ifelse(V(g)$type, pi, 0) # OTUs labels left, metabolites right

# === Edge colors ===
edge_colors <- rep("#006400", ecount(g))  # Dark green

# === Final plot ===
plot(g, 
     layout = layout,
     vertex.label.cex = 0.7,
     edge.width = abs(E(g)$Correlation) * 2,
     edge.color = edge_colors,
     vertex.color = vertex_colors,
     vertex.size = vertex_sizes,
     vertex.label.color = "black",
     vertex.label.dist = vertex_label_dist,
     vertex.label.degree = vertex_label_degree,
     main = "Bipartite Network Analysis of OTUs and Metabolites")
