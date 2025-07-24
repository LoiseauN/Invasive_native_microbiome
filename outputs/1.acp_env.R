#======== PROJECT COM2LIFE ========
## ACP on lakes and environmental variables



param_table_july2021 <- read.csv(here::here("data",
                                 "param_table_july2021.csv"), sep = ";")

# clean param_table
param_table_july2021 <- param_table_july2021[,-1]
rownames(param_table_july2021) <- param_table_july2021$lake

param_table_july2021 <- param_table_july2021 %>%
  filter(lake != "BOI" & lake != "JAB" & lake != "VAI")

# select columns
variables_env <- as.data.frame(param_table_july2021[, c("Temperature_median", "Salinity_median", "Oxygen_median", "pH_median", "Chla_median", "TPN", "TPC", "TPC.TPN", "NH4", "PO4", "NO3_N02")])
rownames(variables_env) <- rownames(param_table_july2021)
variables_env <- variables_env %>%
  dplyr::rename(
    Temperature = Temperature_median,
    Salinite = Salinity_median,
    O2 = Oxygen_median,
    pH = pH_median,
    Chla = Chla_median
  )
# Effectuez l'ACP
acp_result <- FactoMineR::PCA(variables_env, graph = FALSE)

# Affichez un graphique des valeurs propres
factoextra::fviz_eig(acp_result, addlabels = TRUE, ylim = c(0, 50))

# Affichez un graphique des individus
factoextra::fviz_pca_ind(acp_result, col.ind = "cos2", pointsize = 2)

# Affichez un graphique des variables
factoextra::fviz_pca_var(acp_result, col.var = "cos2", pointsize = 2)

biplot <- factoextra::fviz_pca_biplot(acp_result, pointsize = 2, col.var = "contrib")
biplot

path_to_my_object = here::here("figures", "acp_env.png")
ggplot2::ggsave(filename = path_to_my_object, plot = biplot, device = "png")