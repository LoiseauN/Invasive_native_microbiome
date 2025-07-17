#### Script plotcrastination ####
## The aim of this script is to generate a lot of plots to discover the data 
## and to know it better, to see if there is some pattern, correlation, ... ##

# Install necessary packages
install.packages("mFD")
library("mFD")
install.packages("tidyverse")
library("tidyverse")
install.packages("dplyr")
library("dplyr")
install.packages("ggplot2")
library("ggplot2")
install.packages("reshape2")
library(reshape2)
install.packages("magrittr")
library(magrittr)
install.packages("stringr")
library(stringr)



# modifcation de la table d'ASV 
ASV_table_t <- ASV_table_july_2021 %>%
  `rownames<-`(.$X) %>%
  select(-1) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate(labels = rowname) %>%
  rownames_to_column(var = "NouvelleLignes") %>%
  column_to_rownames(var = "rowname") %>%
  select(-NouvelleLignes)
ASV_table_t$lac <- substr(ASV_table_t$labels, 1,3)
# Filtrer les data et spprimer les ASV uniquement présent dans moins de 3 échantillons 
ASV_table <- ASV_table_t[,apply(ASV_table_t > 0,2,sum) > 3]

##Plot about Environmental conditions ##
# Fonction pour créer un boxplot
create_boxplot <- function(data, title) {
  ggplot(data, aes(x = lake, y = value)) +
    geom_point() +
    facet_wrap(~ variable, scales = "free_y") +
    labs(title = title, x = "Lac", y = "Valeurs") +
    theme_minimal()
}

# Pipeline pour les conditions environnementales
param_env_long <- param_table_july2021 %>%
  select(lake, Temperature_median, Salinity_median, Oxygen_median, pH_median, Chla_median) %>%
  as.data.frame() %>%
  melt(id.vars = "lake")

create_boxplot(param_env_long, "Boxplot des Conditions Environnementales par Lac")

# Pipeline pour les paramètres chimiques
param_ch_long <- param_table_july2021 %>%
  select(lake, TPN, TPC, TPC.TPN, NH4, PO4) %>%
  melt(id.vars = "lake")

create_boxplot(param_ch_long, "Boxplot des Paramètres Chimiques par Lac")

## Plot of occurrence of ASV in the Water, Sediments and fishes ##
# Fonction pour créer dataframe occurrence en format long qu'on souhaite plot
process_ASV_table <- function(data_table, pattern) {
  result <- data_table %>%
    filter(str_detect(labels, pattern)) %>%
    group_by(labels) %>%
    pivot_longer(
      cols = starts_with("ASV_"),
      names_to = "ASV",
      values_to = "nb_occurrence"
    ) %>%
    ungroup() %>%
    filter(nb_occurrence > 3500)
  
  return(result)
}

# Fonction pour générer le graphique
plot_occurrence <- function(data, title) {
  data %>%
    ggplot(aes(x = ASV, y = nb_occurrence, fill = labels)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = title,
         x = "ASV",
         y = "Nombre d'Occurrence",
         fill = "Labels") +
    theme_minimal() +
    facet_wrap(~ labels, scales = "free_y", ncol = 2)
}

## BOUCLE LAPPLY POUR CREER LES DF ET GENERER LES PLOT D'OCCCURRENCE ET LES SAVE
process_and_plot <- function(echantillon) {
  df <- process_ASV_table(ASV_table, echantillon)
  plot <- plot_occurrence(df, paste("Comparaison occurrence -", echantillon))
  # Afficher le plot
  print(plot)
  # Sauvegarder le plot en tant que fichier image 
  ggsave(paste("plot_occurrence_", gsub("C:/Users/nicol/Documents/COM2LIFE/M2_COM2LIFE", "_", echantillon), ".png"), plot)
  # Retourner le dataframe
  return(df)
}

# Liste des échantillons
echantillons <- c("BOI.B.GAR.GUT", "BOI.B.W", "CHA.B.GAR.GUT", "CHA.B.PER.GUT", "CHA.B.W", 
                  "CRJ.B.GAR.GUT", "CRJ.B.PER.GUT", "CRJ1.B.W", "CRJ2.B.W", "CTL.B.GAR.GUT", 
                  "CTL.B.PER.GUT", "CTL.B.W", "GDP.B.GAR.GUT", "GDP.B.PER.GUT", "GDP.W.B", 
                  "JAB.B.GAR.GUT", "JAB.B.PER.GUT", "TRI.B.GAR.GUT", "TRI.B.PER.GUT", 
                  "VAI.B.GAR.GUT", "VAI.B.W", "VER.B.GAR.GUT", "VER.B.PER.GUT", "VER.B.W")

# Appliquer la fonction à chaque échantillon avec lapply
resultats <- lapply(echantillons, process_and_plot)


#### PLOT ABONDANCE OF ASV IN EACH COMPARTIMENTS ####
library(RColorBrewer)

## Les 20 ASV les plus abondants pour chaque lac 
# Trouver les indices des colonnes qui commencent par "ASV"
asv_columns_indices <- grep("^ASV", names(ASV_table))
# Boucle sur les noms des lacs
for (lac_name in c("BOI", "CHA", "CRJ", "TRI", "CTL", "GDP", "JAB", "VAI", "VER")) {
  # Filtrer le dataframe pour le lac actuel
  filtered_data <- ASV_table[ASV_table$lac == lac_name, ]
  
  # Sélectionner les colonnes d'ASV
  asv_columns <- filtered_data[, asv_columns_indices]
  
  # Imprimer les dimensions
  print(dim(asv_columns))
  
  # Trouver les 20 valeurs maximales
  top_20_max_values <- head(sort(as.matrix(asv_columns), decreasing = TRUE), 20)
  
  # Imprimer les valeurs maximales 
  print(top_20_max_values)
  
  # Extraire les noms des 20 premières colonnes d'ASV
  top_20_asv_indices <- head(order(-colMeans(asv_columns)), 20)
  top_20_asv_names <- names(asv_columns)[top_20_asv_indices]
  
  # Imprimer les noms 
  print(top_20_asv_names)
  
  # Répéter le nom du lac pour correspondre à la longueur des autres vecteurs
  lac_name_rep <- rep(lac_name, 20)
  
  # Créer un dataframe avec les informations nécessaires
  result_df <- data.frame(
    Lac = lac_name_rep,
    ASV_Noms = top_20_asv_names,
    Abondances = top_20_max_values
  )
  
  # Ajouter le résultat à la liste
  result_list[[lac_name]] <- result_df
}

print(result_list)

# Fonction pour transformer les données pour chaque lac
transform_data <- function(result_df) {
  # Calculer la somme des abondances pour normaliser en pourcentage
  total_abundance <- sum(result_df$Abondances)
  
  # Calculer les pourcentages pour chaque ASV
  result_df$Pourcentage <- (result_df$Abondances / total_abundance) * 100
  
  return(result_df)
}

# Transformer les données pour chaque lac dans la liste
transformed_list <- lapply(result_list, transform_data)

# Concaténer les dataframes pour chaque lac en un seul dataframe
combined_df <- do.call(rbind, transformed_list)

# Créer un graphique à barres empilées
ggplot(combined_df, aes(x = Lac, y = Pourcentage, fill = ASV_Noms)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Pourcentage des 20 ASV par Lac", y = "Pourcentage", x = "Lac") +
  theme_minimal() +
  theme(legend.position = "right", legend.direction = "vertical")

# Fonction pour obtenir les 20 ASV les plus abondants pour un type de lac donné
get_top_20_asv_for_lac <- function(lac_prefix, lac_names, asv_table) {
  result_list <- list()
  
  for (lac_name in lac_names) {
    # Filtrer les lignes où la colonne 'labels' contient le préfixe du lac
    filtered_data <- asv_table[grepl(lac_prefix, asv_table$labels) & asv_table$lac == lac_name, ]
    
    # Trouver les indices des colonnes qui commencent par "ASV"
    asv_columns_indices <- grep("^ASV", colnames(filtered_data))
    
    # Sélectionner les colonnes d'ASV
    asv_columns <- filtered_data[, asv_columns_indices]
    
    # Imprimer les dimensions
    print(dim(asv_columns))
    
    # Trouver les 20 valeurs maximales
    top_20_max_values <- head(sort(as.matrix(asv_columns), decreasing = TRUE), 20)
    
    # Imprimer les valeurs maximales 
    print(top_20_max_values)
    
    # Extraire les noms des 20 premières colonnes d'ASV
    top_20_asv_indices <- head(order(-colMeans(asv_columns)), 20)
    top_20_asv_names <- names(asv_columns)[top_20_asv_indices]
    
    # Imprimer les noms 
    print(top_20_asv_names)
    
    # Répéter le nom du lac pour correspondre à la longueur des autres vecteurs
    lac_name_rep <- rep(lac_name, 20)
    
    # Créer un dataframe avec les informations nécessaires
    result_df <- data.frame(
      Lac = lac_name_rep,
      ASV_Noms = top_20_asv_names,
      Abondances = top_20_max_values
    )
    
    # Ajouter le résultat à la liste
    result_list[[paste(lac_prefix, lac_name, sep = "_")]] <- result_df
  }
  
  return(result_list)
}

# Liste des préfixes de lac et des noms de lac correspondants
lac_prefixes <- c("GAR", "PER")
lac_names_list <- list(c("BOI", "CHA", "CRJ", "TRI", "CTL", "GDP", "JAB", "VAI", "VER"),
                       c("CHA", "CRJ", "TRI", "CTL", "GDP", "JAB", "VER"))

# Utiliser la fonction pour obtenir les résultats pour chaque type de lac
result_list_gar <- get_top_20_asv_for_lac(lac_prefixes[1], lac_names_list[[1]], ASV_table)
result_list_per <- get_top_20_asv_for_lac(lac_prefixes[2], lac_names_list[[2]], ASV_table)

# Fonction pour transformer les données pour chaque lac
transform_data <- function(result_df) {
  # Calculer la somme des abondances pour normaliser en pourcentage
  total_abundance <- sum(result_df$Abondances)
  
  # Calculer les pourcentages pour chaque ASV
  result_df$Pourcentage <- (result_df$Abondances / total_abundance) * 100
  
  return(result_df)
}

# Transformer les données pour chaque lac dans la liste (GAR)
transformed_list_gar <- lapply(result_list_gar, transform_data)

# Transformer les données pour chaque lac dans la liste (PER)
transformed_list_per <- lapply(result_list_per, transform_data)

# Concaténer les dataframes pour chaque lac en un seul dataframe (GAR)
combined_df_gar <- do.call(rbind, transformed_list_gar)

# Concaténer les dataframes pour chaque lac en un seul dataframe (PER)
combined_df_per <- do.call(rbind, transformed_list_per)

# Créer un graphique à barres empilées pour GAR
ggplot(combined_df_gar, aes(x = Lac, y = Pourcentage, fill = ASV_Noms)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ASV_Noms), position = position_stack(vjust = 0.5), color = "white", size = 3) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Pourcentage des 20 ASV par Lac pour LEP", y = "Pourcentage", x = "Lac") +
  theme_minimal() +
  theme(legend.position = "right", legend.direction = "vertical")

# Créer un graphique à barres empilées pour PER
ggplot(combined_df_per, aes(x = Lac, y = Pourcentage, fill = ASV_Noms)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = ASV_Noms), position = position_stack(vjust = 0.5), color = "white", size = 3) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Pourcentage des 20 ASV par Lac pour PER", y = "Pourcentage", x = "Lac") +
  theme_minimal() +
  theme(legend.position = "right", legend.direction = "vertical")

#### DISTANCE DE BRAY-CURTIS ET PCoA ####
install.packages("proxy")
library(proxy)

### PCoA sur matrice de Bray_Curtis de chaque lacs tous compartiments mélangés 
library(RColorBrewer)
library(viridis)
color_palette <- brewer.pal(n = 11, name = "Set3")  

ASV_table %>%
  select_if(is.numeric) %>%
  proxy::dist(method = "bray") %>%
  cmdscale(eig = TRUE) %>%
  `$`("points") %>%
  as.data.frame() %>%
  setNames(c("PCoA1", "PCoA2")) %>%
  rownames_to_column(var = "Sample") %>%
  mutate(Group = sub("\\..*", "", Sample),
         Group = factor(Group)) %>%
  ggplot(aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point() +
  scale_color_manual(values = color_palette) +
  labs(x = "PCoA Axis 1", y = "PCoA Axis 2", title = "PCoA Plot")

### Fonction pour créer les plots PCoA
plot_pcoa <- function(data, prefix) {
  ASV_subset <- data %>%
    rownames_to_column(var = "Sample") %>%
    filter(grepl(prefix, Sample))
  
  pcoa_result <- ASV_subset %>%
    select_if(is.numeric) %>%
    proxy::dist(method = "bray") %>%
    cmdscale(eig = TRUE) %>%
    `$`("points") %>%
    as.data.frame() %>%
    setNames(c("PCoA1", "PCoA2"))
  
  ASV_subset$Color <- substring(ASV_subset$Sample, 1, 2)
  
  color_palette <- brewer.pal(9, "Set1")
  
  pcoa_df <- cbind(ASV_subset[c("Sample", "Color")], pcoa_result)
  
  ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Color)) +
    geom_point() +
    scale_color_manual(values = setNames(color_palette, unique(ASV_subset$Color))) +
    labs(x = "PCoA Axis 1", y = "PCoA Axis 2", title = paste("PCoA Plot -", prefix))
}
# Utilisation de la fonction avec différents préfixes
data <- ASV_table

plot_pcoa(data, "LEP") #PCoA sur matrice de Bray-Curtis de la perche soleil dans chaque lac 
plot_pcoa(data, "PER") #PCoA sur matrice de Bray-curtis de la perche dans chaque lac 
plot_pcoa(data, "ADN") #PCoA sur matrice de Bray-Curtis de l'ADN de l'eau de chaque lac 
plot_pcoa(data, "SED") #PcoA sur matrice de Bray-Curtis de l'ADN du sédiment de chaque lac 

## Fonction pour créer un plot pcoa combiné des 2 espèces 
plot_pcoa_combined <- function(data, prefix) {
  ASV_subset <- data %>%
    rownames_to_column(var = "Sample") %>%
    filter(grepl(prefix, Sample))
  
  pcoa_result <- ASV_subset %>%
    select_if(is.numeric) %>%
    proxy::dist(method = "bray") %>%
    cmdscale(eig = TRUE) %>%
    `$`("points") %>%
    as.data.frame() %>%
    setNames(c("PCoA1", "PCoA2"))
  
  ASV_subset$Color <- substring(ASV_subset$Sample, 1, 2)
  
  pcoa_df <- cbind(ASV_subset[c("Sample", "Color")], pcoa_result, prefix = prefix)
  
  return(pcoa_df)
}

# Utilisation de la fonction avec différents préfixes
data <- ASV_table

per_df <- plot_pcoa_combined(data, "PER")
lep_df <- plot_pcoa_combined(data, "GAR")

# Combinaison des données des deux préfixes
combined_df <- rbind(per_df, lep_df)

# Création du graphique avec les points des deux préfixes
color_palette <- brewer.pal(9, "Set1")
ggplot(combined_df, aes(x = PCoA1, y = PCoA2, color = Color, shape = prefix)) +
  geom_point() +
  scale_color_manual(values = setNames(color_palette, unique(combined_df$Color))) +
  labs(x = "PCoA Axis 1", y = "PCoA Axis 2", title = "Combined PCoA Plot")

