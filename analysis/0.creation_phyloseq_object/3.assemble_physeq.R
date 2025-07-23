#### CLEAN DATA ####
### NETTOYER LES DATA, FILTRER LES ASV ET PREPARER L'OBJET PHYLOSEQ_ALL

## Load phyloseq 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

library(phyloseq)
library(here)

##Importer les data 
asv_table <- read.table(here::here("outputs", "dada2", "asv_table",
                   "asv_table.tsv"), header = TRUE)
metadata_july_2021 <- read.csv(here::here("Data", "metadata_july_2021.csv"), sep = ";")
taxonomy <- read.table(here::here("outputs", "dada2", "asv_table", "taxonomy.tsv"), sep="\t", header = TRUE)
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
