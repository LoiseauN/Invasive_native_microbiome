#====== CREATE PHYSEQ OBJECT =========
#### CLEAN DATA ####
### NETTOYER LES DATA, FILTRER LES ASV ET PREPARER L'OBJET PHYLOSEQ_ALL

## Load phyloseq 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

library(phyloseq)
library(here)

##Importer les data 
asv_table <- read.table(here::here("Data", "asv_table",
                                   "asv_table.tsv"), header = TRUE)
metadata_july_2021 <- read.csv(here::here("Data", "metadata_july_2021.csv"), sep = ";")
taxonomy <- read.table(here::here("Data", "asv_table", "taxonomy.tsv"), sep="\t", header = TRUE)
param_table_july_2021 <- read.csv(here::here("Data", "param_table_july2021.csv"), sep = ";")

## passer la première colonne en nom de lignes
rownames(metadata_july_2021) <-metadata_july_2021$sample
metadata_july_2021 <- metadata_july_2021[,-1]
rownames(asv_table) <- asv_table$asv
asv_table <- asv_table[,-1]
rownames(taxonomy) <- taxonomy$asv
taxonomy <- taxonomy[,-1]
rownames(param_table_july_2021) <- param_table_july_2021$X
param_table_july_2021 <- param_table_july_2021[,-1]
# merge les 2 df metadata et param_table pour avoir qu'une grande table qui regroupe le contexte du sampling
context <- merge(metadata_july_2021, param_table_july_2021, by = c("lake"), all.x = TRUE)
rownames(context) <-context$sample

## remplacer les tirets par des points 
new_row <- gsub("-", ".", rownames(context))
rownames(context) <- new_row

## Exporter le DataFrame vers un fichier texte
write.table(context, here::here("context.txt"), sep = "\t", quote = FALSE)

# add the mock in df context
df_tmp <- data.frame(matrix(NA, nrow = 1, ncol = 24))
colnames(df_tmp) <- colnames(context)
row.names(df_tmp) <- "T.extract.2.ADN1"
context <- rbind(context, df_tmp)

#Check sample file 
setdiff(x = colnames(asv_table),
        y = rownames(context))

## Assemble the physeq object 
physeq <- phyloseq::phyloseq(
  phyloseq::otu_table(asv_table, taxa_are_rows = TRUE),
  phyloseq::tax_table(as.matrix(taxonomy)),
  phyloseq::sample_data(context),
  phyloseq::refseq(asv_seq)
)
path_to_my_object = here::here("Data","mon_objet_physeq_17_06.rds")
saveRDS(physeq, file = path_to_my_object)

# Load my physeq object 
physeq <- readRDS(here::here("data",
                             "mon_objet_physeq_17_06.rds"))

# Turn my matrix to have ASV in column and samples in rows
physeq@otu_table <- t(as.matrix(physeq@otu_table))

# Check that it is correct 
physeq@otu_table[1:10,1:6]

# Remove Eukaryota Archaea Chloroplast and Mitochondria in all the samples 
physeq_sub <- subset_taxa(physeq, (Kingdom!="Eukaryota")|is.na(Kingdom))
physeq_sub <- subset_taxa(physeq_sub, (Kingdom!="d__Archaea")|is.na(Kingdom))
physeq_sub <- subset_taxa(physeq_sub, (Order!="Chloroplast")|is.na(Order))
physeq_sub <- subset_taxa(physeq_sub, (Family!="Mitochondria")|is.na(Family))

##Remove mocks
# Create a new df without the row 133
new_physeq_otu_table <- subset(physeq_sub@otu_table, row.names(physeq_sub@otu_table) != "T.extract.2.ADN1")

# reassign the new df at physeq_sub@otu_table
physeq_sub@otu_table <- new_physeq_otu_table

# Create a new df without the row 133
new_physeq_sam_table <- subset(physeq_sub@sam_data, row.names(physeq_sub@sam_data) != "T.extract.2.ADN1")

# reassign the new df at physeq_sub@sam-data
physeq_sub@sam_data <- new_physeq_sam_table

physeq_sub@tax_table <- physeq_sub@tax_table[rownames(physeq_sub@tax_table) %in% colnames(physeq_sub@otu_table),]
physeq_sub@otu_table <- physeq_sub@otu_table[rownames(physeq_sub@otu_table) %in% colnames(physeq_sub@tax_table),]

# Sum otu_table rows 
relative_abundance <- physeq_sub@otu_table / rowSums(physeq_sub@otu_table)

# Set threshold (1%)
threshold <- 0.01

# Identify ASVs representing at least 1% of reads in at least one sample
asvs_to_keep <- colSums(relative_abundance >= threshold) > 0

# Filter these ASVs
physeq_sub@otu_table <- physeq_sub@otu_table[, asvs_to_keep]
physeq_sub@tax_table <- physeq_sub@tax_table[rownames(physeq_sub@tax_table) %in% colnames(physeq_sub@otu_table),]

# Rarefies the samples
# set the seed for random sampling
# it allows reproductibility
set.seed(10000)
#keep sample with more than 5500 reads
physeq_rar <- rarefy_even_depth(physeq_sub, sample.size = 5500, rngseed = TRUE)
rowSums(physeq_rar@otu_table@.Data) #how many reads per sample

# Make a tree 
# Align the sequences 
aln <- refseq(physeq_rar) |>
  DECIPHER::AlignSeqs(anchor = NA)

### Infering the phylogenetic tree

# We will infer a phylogenetic from our alignement using the library `phangorn`.
# First, let's convert our `DNAStringSet` alignment to the `phangorn` `phyDat` format.
phang_align <- as.matrix(aln) |> phangorn::phyDat(type = "DNA")

# Then, we compute pairwise distances of our aligned sequences using equal base frequencies (JC69 model used by default).
dm <- phangorn::dist.ml(phang_align, model = "JC69")

# Finally, we reconstruct a neighbour joining tree. 
treeNJ <- phangorn::NJ(dm)

# We need the tree to be rooted for future analysis.
# We can do that using the function `phangorn::midpoint()`
treeNJ <- phangorn::midpoint(tree = treeNJ)

# Once we have a rooted tree, we can add it to the phyloseq object.
physeq <- phyloseq::merge_phyloseq(physeq_rar,treeNJ)

# Keep only PER et LEP in physeq objet
physeq_filtered <- subset_samples(physeq, origin %in% c("LEP", "PER"))

# Save my phyloseq object 
path_to_my_object = here::here("data","mon_objet_physeq_17_06.rds")
saveRDS(physeq, file = path_to_my_object)

#== Rarefaction curves ===
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