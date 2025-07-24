#======== PROJECT COM2LIFE ========
## LEFse and LDA to see which ASVs is abundant in each species

# Load libraries 
library(microbiomeMarker)
library(ggplot2)

#### LEFse et LDA
# Load my physeq object
physeq_filtered <- readRDS(here::here("data",
                                      "mon_objet_physeq_filtered.rds"))

#LEFSE
mm_lefse <- microbiomeMarker::run_lefse(physeq_filtered, norm = "CPM",
                                        wilcoxon_cutoff = 0.01,
                                        group = "origin",
                                        taxa_rank = "none",
                                        kw_cutoff = 0.01,
                                        multigrp_strat = TRUE,
                                        lda_cutoff = 4)

mm_lefse_table <- data.frame(mm_lefse@marker_table)
mm_lefse_table

p_LDAsc <- microbiomeMarker::plot_ef_bar(mm_lefse)
p_LDAsc <- p_LDAsc + 
  ggplot2::scale_fill_manual(
    values = c("LEP" = "darkorange", "PER" = "#A90C38"),
    labels = c("LEP" = "Lepomis gibbosus", "PER" = "Perca fluviatilis")
  ) +
  ggplot2::scale_color_manual(
    values = c("LEP" = "darkorange", "PER" = "#A90C38"),
    labels = c("LEP" = "Lepomis gibbosus", "PER" = "Perca fluviatilis")
  )

# Retrieve graph data from p_LDAsc
y_labs <- ggplot2::ggplot_build(p_LDAsc)$layout$panel_params[[1]]$y$get_labels()

p_abd <- microbiomeMarker::plot_abundance(mm_lefse, group = "origin") +
  ggplot2::scale_y_discrete(limits = y_labs) +
  ggplot2::scale_fill_manual(
    values = c("LEP" = "darkorange", "PER" = "#A90C38"),
    labels = c("LEP" = "Lepomis gibbosus", "PER" = "Perca fluviatilis")
  ) +
  ggplot2::scale_color_manual(
    values = c("LEP" = "darkorange", "PER" = "#A90C38"),
    labels = c("LEP" = "Lepomis gibbosus", "PER" = "Perca fluviatilis")
  )

# Display graphics
LEFse_lda <- gridExtra::grid.arrange(p_LDAsc, p_abd, nrow = 1)

path_to_my_object = here::here("Figures", "hill", "lefse_lda.png")
ggplot2::ggsave(filename = path_to_my_object, plot = LEFse_lda, device = "png")
