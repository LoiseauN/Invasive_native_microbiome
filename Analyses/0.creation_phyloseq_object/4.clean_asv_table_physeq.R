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