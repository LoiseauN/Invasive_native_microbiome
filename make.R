
#' Run the Entire Project
#'
#' This script runs the entire project and produces all figures present in the
#' Navarro et al paper.
#'
#' @author Alice Navarro  \email{alice.navarro@@hotmail.fr},
#'         Nicolas Loiseau, \email{nicolas.loiseau@@cnrs.fr},
#'



## Parameters ----
rm(list=ls())
set.seed(666)

if (!("here" %in% installed.packages())) install.packages("here")


#------------------Running code------------------------

source(here::here("analyses", "1.compo_indiv_lake.R"))
source(here::here("analyses", "2.composition_water.R"))
source(here::here("analyses", "3.estimates_alpha_diversity.R"))
source(here::here("analyses", "4.0.estimates_beta_diversity.R"))
source(here::here("analyses", "4.1.estimates_beta_diversity_metabo.R"))
source(here::here("analyses", "4.2.permanova_asvs.R"))
source(here::here("analyses", "4.3.permanova_metabo.R"))
source(here::here("analyses", "5.pcoA_dbRDA_beta_diversity.R"))
source(here::here("analyses", "6.statistical_test_beta_diss.R"))
source(here::here("analyses", "7.core_microbiota.R"))
source(here::here("analyses", "8.model_glmm_alpha.R"))
source(here::here("analyses", "9.model_gdm_beta_per.R"))
source(here::here("analyses", "10.model_gdm_beta_lep.R"))


