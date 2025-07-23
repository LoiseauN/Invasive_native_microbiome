#======== PROJECT COM2LIFE ========
## Modelisation GLMM on alpha diversity estimates in "estimates_alpha_diversity.R" 
## First load the libraries. Do a correlation matrix to chosse environmeental variables 
##### Here, we choose 5 variables 

#load library
library(gridExtra)

# Load metadata object
metadata_alpha <- readRDS(here::here("data",
                                     "metadata_alpha.rds"))

# ====== CORRELATION MATRIX =======
cor_metadata <- cor(metadata_alpha[,c("Temperature_median", "Oxygen_median", "pH_median", "Chla_median", "Salinity_median", "TPC", "TPN",
                                      "TPC.TPN", "NH4", "PO4", "NO3_N02")], method = "spearman")

#function to do test te collineraity between variables
cor_mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p_mat <- matrix(NA, n, n)
  diag(p_mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method = "spearman", ...)
      p_mat[i, j] <- p_mat[j, i] <- tmp$p.value
    }
  }
  colnames(p_mat) <- rownames(p_mat) <- colnames(mat)
  p_mat
}

# matrix of the p-value of the correlation
p_mat <- cor_mtest(metadata_alpha[,c("Temperature_median", "Oxygen_median", "pH_median", "Chla_median", "Salinity_median", "TPC", "TPN",
                                     "TPC.TPN", "NH4", "PO4", "NO3_N02")])
# Leave blank on no significant coefficient
corrplot::corrplot(cor_metadata,
                   type = "upper",
                   order = "hclust",
                   p.mat = p_mat,
                   sig.level = 0.05,
                   insig = "blank",
                   method = "number")

# VIF with all variables to see which variable to keep 
car::vif(glm( taxo_q0 ~ Temperature_median + Chla_median + Salinity_median + pH_median + Oxygen_median , 
              data = metadata_alpha, 
              family = gaussian))

# ====== MODEL GLMM ======
### on alpha diversity ####
# 4 models for each Hill indices (taxo_q, taxo_q, phylo_q and phylo_q)

# MODEL TAXO_Q0
model_glm_taxo_q0 <- glmmTMB::glmmTMB( taxo_q0 ~ Temperature_median + Chla_median + Salinity_median + pH_median + Oxygen_median + (1|lake) +(1|origin),
                                       data = metadata_alpha,
                                       family = gaussian )
# Residuals plots
plot_1 <- sjPlot::plot_model(model_glm_taxo_q0, show.values = TRUE, type = "pred", terms = c("Temperature_median"))
plot_2 <- sjPlot::plot_model(model_glm_taxo_q0, show.values = TRUE, type = "pred", terms = c("Chla_median"))
plot_3 <- sjPlot::plot_model(model_glm_taxo_q0, show.values = TRUE, type = "pred", terms = c("Salinity_median"))
plot_4 <-sjPlot::plot_model(model_glm_taxo_q0, show.values = TRUE, type = "pred", terms = c("pH_median"))
plot_5 <-sjPlot::plot_model(model_glm_taxo_q0, show.values = TRUE, type = "pred", terms = c("Oxygen_median"))
# put all the plot on same page 
grid.arrange(plot_1, plot_2, plot_3, plot_4, plot_5)
#results of model
sjPlot::tab_model(model_glm_taxo_q0, show.aic = TRUE)

# MODEL TAXO_Q1
model_glm_taxo_q1 <- glmmTMB::glmmTMB( taxo_q1 ~ Temperature_median + Chla_median + Salinity_median + pH_median+ Oxygen_median + (1|lake) +(1|origin),
                                   data = metadata_alpha,
                                   family = gaussian )
# Residuals plots
plot_1 <- sjPlot::plot_model(model_glm_taxo_q1, show.values = TRUE, type = "pred", terms = c("Temperature_median"))
plot_2 <- sjPlot::plot_model(model_glm_taxo_q1, show.values = TRUE, type = "pred", terms = c("Chla_median"))
plot_3 <- sjPlot::plot_model(model_glm_taxo_q1, show.values = TRUE, type = "pred", terms = c("Salinity_median"))
plot_4 <-sjPlot::plot_model(model_glm_taxo_q1, show.values = TRUE, type = "pred", terms = c("pH_median"))
plot_5 <-sjPlot::plot_model(model_glm_taxo_q1, show.values = TRUE, type = "pred", terms = c("Oxygen_median"))
# put all the plot on same page 
grid.arrange(plot_1, plot_2, plot_3, plot_4, plot_5)
#results of model
sjPlot::tab_model(model_glm_taxo_q1, show.aic = TRUE)

# MODEL PHYLO_Q0
model_glm_phylo_q0 <- glmmTMB::glmmTMB( phylo_q0 ~ Temperature_median + Chla_median + Salinity_median + pH_median + Oxygen_median + (1|lake) +(1|origin),
                                    data = metadata_alpha,
                                    family = gaussian )
# Residuals plots
plot_1 <- sjPlot::plot_model(model_glm_phylo_q0, show.values = TRUE, type = "pred", terms = c("Temperature_median"))
plot_2 <- sjPlot::plot_model(model_glm_phylo_q0, show.values = TRUE, type = "pred", terms = c("Chla_median"))
plot_3 <- sjPlot::plot_model(model_glm_phylo_q0, show.values = TRUE, type = "pred", terms = c("Salinity_median"))
plot_4 <-sjPlot::plot_model(model_glm_phylo_q0, show.values = TRUE, type = "pred", terms = c("pH_median"))
plot_5 <-sjPlot::plot_model(model_glm_phylo_q0, show.values = TRUE, type = "pred", terms = c("Oxygen_median"))
# put all the plot on same page 
grid.arrange(plot_1, plot_2, plot_3, plot_4, plot_5)
# results of model
sjPlot::tab_model(model_glm_phylo_q0, show.aic = TRUE)

# MODEL PHYLO_Q1
model_glm_phylo_q1 <- glmmTMB::glmmTMB( phylo_q1 ~ Temperature_median + Chla_median + Salinity_median + pH_median+ Oxygen_median + (1|lake) +(1|origin),
                                    data = metadata_alpha,
                                    family = gaussian )
# Residuals plots
plot_1 <- sjPlot::plot_model(model_glm_phylo_q1, show.values = TRUE, type = "pred", terms = c("Temperature_median"))
plot_2 <- sjPlot::plot_model(model_glm_phylo_q1, show.values = TRUE, type = "pred", terms = c("Chla_median"))
plot_3 <- sjPlot::plot_model(model_glm_phylo_q1, show.values = TRUE, type = "pred", terms = c("Salinity_median"))
plot_4 <-sjPlot::plot_model(model_glm_phylo_q1, show.values = TRUE, type = "pred", terms = c("pH_median"))
plot_5 <-sjPlot::plot_model(model_glm_phylo_q1, show.values = TRUE, type = "pred", terms = c("Oxygen_median"))
# put all the plot on same page 
grid.arrange(plot_1, plot_2, plot_3, plot_4, plot_5)
# results of model
sjPlot::tab_model(model_glm_phylo_q1, show.aic = TRUE)