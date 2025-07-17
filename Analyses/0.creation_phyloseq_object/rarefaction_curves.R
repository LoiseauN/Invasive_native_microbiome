ggrare <- function(physeq, step = 10, label = NULL, color = NULL,
                   plot = TRUE, parallel = FALSE, se = TRUE) {
  
  require("ggplot2")
  
  x <- as(phyloseq::otu_table(physeq), "matrix")
  if (phyloseq::taxa_are_rows(physeq)) x <- t(x)
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  
  df <- do.call(rbind, out)
  
  ## Get sample data 
  if (!is.null(phyloseq::sample_data(physeq, FALSE))) {
    sdf <- as(phyloseq::sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if (length(color) > 1) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if (length(label) > 1) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot2::ggplot(data = data,
              ggplot2::aes_string(x = "Size", y = ".S",
                         group = "Sample", color = color)) +
    ggplot2::labs(x = "Sample Size", y = "Species Richness")
  
  if (!is.null(label)) {
    p <- p + geom_text(data = labels,
                       aes_string(x = "x", y = "y",
                                  label = label, color = color),
                       size = 4, hjust = 0)
  }
  
  p <- p + ggplot2::geom_line()
  
  if (se) { ## add standard error if available
    p <- p +
      geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se",
                             color = NULL, fill = color),
                  alpha = 0.2)
  }
  
  if (plot) {
    plot(p)
  }
  
  invisible(p)
  
}

physeq <- readRDS(here::here("data",
                             "mon_objet_physeq_filtered.rds"))
#Make rarefaction curves & Add min sample size line
library(ggplot2)
library(phyloseq)

# Définir les couleurs associées aux noms des échantillons
colors <- c(
  "CHA" = "black",
  "CRJ" = "#CFEDD5",
  "TRI" = "#B3E0A6",
  "CTL" = "#40AA5F",
  "GDP" = "#004616",
  "VER" = "#27763D"
)

# Créer le graphique avec ggrare et appliquer les couleurs spécifiques et les personnalisations
rarefaction_curves <- ggrare(physeq, step = 10, color = "sample_color", se = FALSE) +
  geom_vline(xintercept = min(sample_sums(physeq)), color = "red") +
  scale_color_identity() +
  theme_bw() +
  theme(legend.position = "none")

path_to_my_object = here::here("Figures", "rarefaction_curves.png")
ggsave(filename = path_to_my_object, plot = rarefaction_curves, device = "png")
