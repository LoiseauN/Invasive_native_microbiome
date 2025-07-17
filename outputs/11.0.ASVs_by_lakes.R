#======== PROJECT COM2LIFE ========
## This script subsample the ASVs of each compartiments

# Load my physeq object 
physeq <- readRDS(here::here("data",
                             "mon_objet_physeq_17_06.rds"))
# Filter physeq to keep only samples of fish
physeq_filtered_all <- subset_samples(physeq, !(grepl("BOI|JAB|VAI|SED", sample_names(physeq))))
rows_to_remove <- grepl("BOI|JAB|VAI|SED", rownames(sample_data(physeq_filtered_all)))

# keep the sam data
physeq_filtered_all@sam_data <- sample_data(physeq_filtered_all)[!rows_to_remove, ]

# Modify lakes's name
# Access samples name
sample_names_current <- sample_names(physeq_filtered_all)
# Modify names
sample_names_modified <- sample_names_current %>%
  gsub("CHA", "CSM", .) %>%
  gsub("CRJ", "CERL", .) %>% 
  gsub("TRI", "CERS", .) %>%
  gsub("CTL", "CRE", .) %>%
  gsub("GDP", "LGP", .) %>%
  gsub("VER", "VSS", .)
# Apply new names
sample_names(physeq_filtered_all) <- sample_names_modified
# Modify names for sam_data
sample_data <- sample_data(physeq_filtered_all)
# for lake 
if ("lake" %in% colnames(sample_data)) {
  sample_data[["lake"]] <- sample_data[["lake"]] %>%
    gsub("CHA", "CSM", .) %>%
    gsub("CRJ1", "CERL", .) %>% 
    gsub("CERS", "CERS", .) %>%
    gsub("CTL", "CRE", .) %>%
    gsub("GDP", "LGP", .) %>%
    gsub("VER", "VSS", .)
}
# for lake_origin
if ("lake_origin" %in% colnames(sample_data)) {
  sample_data[["lake_origin"]] <- sample_data[["lake_origin"]] %>%
    gsub("CHA", "CSM", .) %>%
    gsub("CRJ1", "CERL", .) %>% 
    gsub("CERS", "CERS", .) %>%
    gsub("CTL", "CRE", .) %>%
    gsub("GDP", "LGP", .) %>%
    gsub("VER", "VSS", .)
}
#for samples 
if ("samples" %in% colnames(sample_data)) {
  sample_data[["samples"]] <- sample_data[["samples"]] %>%
    gsub("CHA", "CSM", .) %>%
    gsub("CRJ", "CERL", .) %>% 
    gsub("TRI", "CERS", .) %>%
    gsub("CTL", "CRE", .) %>%
    gsub("GDP", "LGP", .) %>%
    gsub("VER", "VSS", .)
}
# Update sample_data in phyloseq object
sample_data(physeq_filtered_all) <- sample_data

# Save
path_to_my_object = here::here("Data", "mon_objet_physeq_filtered_all.rds")
saveRDS(physeq_filtered_all, file = path_to_my_object)

#Load the needed library
library(VennDiagram)
library(venn)
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)

# get the asv table from the phyloseq objecf
asv_table <- as.data.frame(physeq@otu_table)
asv_table_filtered <- asv_table[!grepl("BOI|JAB|VAI|SED", rownames(asv_table)),]

# subsampling the asv 
water_asv <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("ADN1", rownames(asv_table_filtered)),]) > 0]
lep_asv <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("GAR", rownames(asv_table_filtered)),]) > 0]
per_asv <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("PER", rownames(asv_table_filtered)),]) > 0]

## CERL
# subsampling the asv 
water_asv_CERL <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CERL.*ADN1|ADN1.*CERL", rownames(asv_table_filtered)), ]) != 0]
lep_asv_CERL <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CERL.*GAR|GAR.*CERL", rownames(asv_table_filtered)), ]) != 0]
per_asv_CERL <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CERL.*PER|PER.*CERL", rownames(asv_table_filtered)), ]) != 0]

## CERS
# subsampling the asv 
water_asv_CERS <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CERS.*ADN1|ADN1.*CERS", rownames(asv_table_filtered)), ]) != 0]
lep_asv_CERS <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CERS.*GAR|GAR.*CERS", rownames(asv_table_filtered)), ]) != 0]
per_asv_CERS <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CERS.*PER|PER.*CERS", rownames(asv_table_filtered)), ]) != 0]

## CRE
# subsampling the asv
water_asv_CRE <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CRE.*ADN1|ADN1.*CRE", rownames(asv_table_filtered)), ]) != 0]
lep_asv_CRE <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CRE.*GAR|GAR.*CRE", rownames(asv_table_filtered)), ]) != 0]
per_asv_CRE <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CRE.*PER|PER.*CRE", rownames(asv_table_filtered)), ]) != 0]

## LGP
# subsampling the asv
water_asv_LGP <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("LGP.*ADN1|ADN1.*LGP", rownames(asv_table_filtered)), ]) != 0]
lep_asv_LGP <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("LGP.*GAR|GAR.*LGP", rownames(asv_table_filtered)), ]) != 0]
per_asv_LGP <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("LGP.*PER|PER.*LGP", rownames(asv_table_filtered)), ]) != 0]

## VSS
# subsampling the asv
water_asv_VSS <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("VSS.*ADN1|ADN1.*VSS", rownames(asv_table_filtered)), ]) != 0]
lep_asv_VSS <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("VSS.*GAR|GAR.*VSS", rownames(asv_table_filtered)), ]) != 0]
per_asv_VSS <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("VSS.*PER|PER.*VSS", rownames(asv_table_filtered)), ]) != 0]

## CSM
# subsampling the asv
water_asv_CSM <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CSM.*ADN1|ADN1.*CSM", rownames(asv_table_filtered)), ]) != 0]
lep_asv_CSM <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CSM.*GAR|ADN1.*GAR", rownames(asv_table_filtered)), ]) !=0]
per_asv_CSM <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CSM.*PER|ADN1.*PER", rownames(asv_table_filtered)), ]) != 0]
