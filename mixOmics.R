## install BiocManager if not installed if (!requireNamespace("BiocManager", quietly = TRUE))     install.packages("BiocManager") 
## install mixOmics 
BiocManager::install('mixOmics')
## library 
library(mixOmics)

# load data
physeq_per <- readRDS(here::here("Data",
                                 "mon_objet_physeq_per.rds"))
physeq_lep <- readRDS(here::here("Data",
                                      "mon_objet_physeq_lep.rds"))
annotated_metabo_per <- readRDS(here::here("Data",
                                      "annotated_metabo_per.rds"))
annotated_metabo_lep <- readRDS(here::here("Data",
                                           "annotated_metabo_lep.rds"))
sample_ID_conv <- readRDS(here::here("Data",
                                           "sample_ID_conv.rds"))
# otu_table
otu_table_per <- as.data.frame(physeq_per@otu_table)
otu_table_lep <- as.data.frame(physeq_lep@otu_table)

# Définir la fonction pour renommer les échantillons basés sur les rownames
rename_metabolomic_samples <- function(annotated_metabo, sample_id_conv, id_metabolomics_col, id_sequencing_col) {
  # Convertir en data frame si ce n'est pas déjà le cas
  annotated_metabo <- as.data.frame(annotated_metabo)
  
  # Extraire les rownames actuels
  sample_ids <- rownames(annotated_metabo)
  
  # 1. Filtrer sample_ID_conv pour ne garder que les IDs présents dans les rownames de annotated_metabo
  sample_id_filtered <- sample_id_conv %>%
    filter(!!sym(id_metabolomics_col) %in% sample_ids)
  
  # 2. Créer un vecteur de correspondance pour le renommage
  rename_mapping <- setNames(sample_id_filtered[[id_sequencing_col]], sample_id_filtered[[id_metabolomics_col]])
  
  # 3. Renommer les rownames de annotated_metabo en utilisant le vecteur de correspondance
  rownames(annotated_metabo) <- ifelse(sample_ids %in% names(rename_mapping), rename_mapping[sample_ids], sample_ids)
  
  # 4. Retourner le dataframe renommé
  return(annotated_metabo)
}
  
# Renommer les échantillons
annotated_metabo_lep<- rename_metabolomic_samples(
  annotated_metabo = annotated_metabo_lep,
  sample_id_conv = sample_ID_conv,
  id_metabolomics_col = "ID_Metabolomics",
  id_sequencing_col = "ID_Sequencing"
)
annotated_metabo_lep <- subset(annotated_metabo_lep, !(rownames(annotated_metabo_lep) %in% c("LGP-B-GAR-GUT3", "LGP-B-GAR-GUT8")))
rownames(annotated_metabo_lep) <- gsub("\\-", ".", rownames(annotated_metabo_lep))

annotated_metabo_per<- rename_metabolomic_samples(
  annotated_metabo = annotated_metabo_per,
  sample_id_conv = sample_ID_conv,
  id_metabolomics_col = "ID_Metabolomics",
  id_sequencing_col = "ID_Sequencing"
)
annotated_metabo_per <- annotated_metabo_per %>%
  filter(!(rownames(annotated_metabo_per) == "CERL-B-PER-GUT1"))
rownames(annotated_metabo_per) <- gsub("\\-", ".", rownames(annotated_metabo_per))


# NETWORK ANALYSIS 
# ===== LEPOMIS GIBBOSUS =====
# Convertir les noms des échantillons en facteurs pour faciliter l'appariement
otu_samples <- rownames(otu_table_lep)
metabo_samples <- rownames(annotated_metabo_lep)

# Trouver les échantillons communs
common_samples <- intersect(otu_samples, metabo_samples)

# Filtrer les tables pour conserver uniquement les échantillons communs
otu_table_filtered <- otu_table_lep[common_samples, ]
metabo_table_filtered <- annotated_metabo_lep[common_samples, ]

# Vérifier que les noms des échantillons correspondent maintenant
all(rownames(otu_table_filtered) == rownames(metabo_table_filtered))

# Supprimer les colonnes avec une déviation standard nulle dans les deux tables
otu_table_filtered <- otu_table_filtered[, apply(otu_table_filtered, 2, sd) != 0]
metabo_table_filtered <- metabo_table_filtered[, apply(metabo_table_filtered, 2, sd) != 0]

# Vérifier les dimensions après la suppression
dim(otu_table_filtered)
dim(metabo_table_filtered)

tax_lep <- as.data.frame(physeq_lep@tax_table)

# Création d'une nouvelle colonne avec le nom du genre ou de la famille
tax_lep$Label <- ifelse(is.na(tax_lep$Genus), tax_lep$Family, tax_lep$Genus)

# Vérifier les 10 premières lignes pour s'assurer que la nouvelle colonne est correcte
head(tax_lep, 10)
# Remplacer les noms des colonnes dans la table filtrée par les labels de tax_lep
colnames(otu_table_filtered) <- tax_lep[colnames(otu_table_filtered), "Label"]

# Charger les packages nécessaires
library(igraph)
library(dplyr)

# Calculer les corrélations entre les OTUs et les métabolites
cor_matrix <- cor(otu_table_filtered, metabo_table_filtered)

# Transformer la matrice de corrélation en un dataframe pour la visualisation
cor_df <- as.data.frame(as.table(cor_matrix))
colnames(cor_df) <- c("OTU", "Metabolite", "Correlation")

# Définir un seuil pour la corrélation
threshold <- 0.7
cor_df_filtered <- cor_df %>% filter(abs(Correlation) > threshold)

# Vérifier combien de connexions restent après filtrage
cat("Nombre de connexions après filtrage:", nrow(cor_df_filtered), "\n")

# Créer un graphe bipartite à partir du dataframe de corrélations filtrées
bipartite_edges <- cor_df_filtered[, c("OTU", "Metabolite", "Correlation")]
g <- graph_from_data_frame(bipartite_edges, directed = FALSE)

# Ajouter un attribut pour distinguer les types de nœuds (OTU ou Metabolite)
V(g)$type <- bipartite.mapping(g)$type

# Définir des couleurs pour les types de nœuds
vertex_colors <- ifelse(V(g)$type, "lightblue", "salmon")

# Ajuster les tailles des nœuds pour améliorer la lisibilité
vertex_sizes <- ifelse(V(g)$type, 5, 7)

# Visualiser le réseau bipartite
network_lep <- plot(g, 
     vertex.label.cex = 0.7, 
     edge.width = abs(E(g)$Correlation) * 2,  # Échelle des largeurs des arêtes
     vertex.color = vertex_colors, 
     vertex.size = vertex_sizes,
     vertex.label.color = "black",
     main = "Bipartite Network Analysis of OTUs and Metabolites of Lepomis gibbosus",
     layout = layout_with_fr)  # Utiliser un layout de Fruchterman-Reingold pour une meilleure lisibilité

# ===== PERCA FLUVIATILIS =====
# Convertir les noms des échantillons en facteurs pour faciliter l'appariement
otu_samples <- rownames(otu_table_per)
metabo_samples <- rownames(annotated_metabo_per)

# Trouver les échantillons communs
common_samples <- intersect(otu_samples, metabo_samples)

# Filtrer les tables pour conserver uniquement les échantillons communs
otu_table_filtered <- otu_table_per[common_samples, ]
metabo_table_filtered <- annotated_metabo_per[common_samples, ]

# Vérifier que les noms des échantillons correspondent maintenant
all(rownames(otu_table_filtered) == rownames(metabo_table_filtered))

# Supprimer les colonnes avec une déviation standard nulle dans les deux tables
otu_table_filtered <- otu_table_filtered[, apply(otu_table_filtered, 2, sd) != 0]
metabo_table_filtered <- metabo_table_filtered[, apply(metabo_table_filtered, 2, sd) != 0]

# Vérifier les dimensions après la suppression
dim(otu_table_filtered)
dim(metabo_table_filtered)

tax_per <- as.data.frame(physeq_per@tax_table)

# Création d'une nouvelle colonne avec le nom du genre ou de la famille
tax_per$Label <- ifelse(is.na(tax_per$Genus), tax_per$Family, tax_per$Genus)

# Vérifier les 10 premières lignes pour s'assurer que la nouvelle colonne est correcte
head(tax_per, 10)
# Remplacer les noms des colonnes dans la table filtrée par les labels de tax_lep
colnames(otu_table_filtered) <- tax_per[colnames(otu_table_filtered), "Label"]

# Charger les packages nécessaires
library(igraph)
library(dplyr)

# Calculer les corrélations entre les OTUs et les métabolites
cor_matrix <- cor(otu_table_filtered, metabo_table_filtered)

# Transformer la matrice de corrélation en un dataframe pour la visualisation
cor_df <- as.data.frame(as.table(cor_matrix))
colnames(cor_df) <- c("OTU", "Metabolite", "Correlation")

# Définir un seuil pour la corrélation
threshold <- 0.7
cor_df_filtered <- cor_df %>% filter(abs(Correlation) > threshold)

# Vérifier combien de connexions restent après filtrage
cat("Nombre de connexions après filtrage:", nrow(cor_df_filtered), "\n")

# Créer un graphe bipartite à partir du dataframe de corrélations filtrées
bipartite_edges <- cor_df_filtered[, c("OTU", "Metabolite", "Correlation")]
g <- graph_from_data_frame(bipartite_edges, directed = FALSE)

# Ajouter un attribut pour distinguer les types de nœuds (OTU ou Metabolite)
V(g)$type <- bipartite.mapping(g)$type

# Définir des couleurs pour les types de nœuds
vertex_colors <- ifelse(V(g)$type, "lightblue", "salmon")

# Ajuster les tailles des nœuds pour améliorer la lisibilité
vertex_sizes <- ifelse(V(g)$type, 5, 7)

# Visualiser le réseau bipartite
network_per <- plot(g, 
     vertex.label.cex = 0.7, 
     edge.width = abs(E(g)$Correlation) * 2,  # Échelle des largeurs des arêtes
     vertex.color = vertex_colors, 
     vertex.size = vertex_sizes,
     vertex.label.color = "black",
     main = "Bipartite Network Analysis of OTUs and Metabolites of Perca fluviatilis",
     layout = layout_with_fr)  # Utiliser un layout de Fruchterman-Reingold pour une meilleure lisibilité


