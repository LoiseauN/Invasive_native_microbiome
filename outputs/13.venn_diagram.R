#======== PROJECT COM2LIFE ========
## Venn diagramm ##


# Load my physeq object 
physeq <- readRDS(here::here("data",
                             "mon_objet_physeq_filtered_all.rds"))

# Load the needed library
library(VennDiagram)
library(venn)
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)

# get the asv table from the phyloseq object
asv_table <- as.data.frame(physeq@otu_table)

# subsampling the asv 
water_asv <- colnames(asv_table)[colSums(asv_table[grep("ADN1", rownames(asv_table)),]) > 0]
lep_asv <- colnames(asv_table)[colSums(asv_table[grep("GAR", rownames(asv_table)),]) > 0]
per_asv <- colnames(asv_table)[colSums(asv_table[grep("PER", rownames(asv_table)),]) > 0]

### Venn Diagram with all the compartments
#create a list with all subsampling asv
x1 <- list(water = water_asv, perca_fluviatilis = per_asv, lepomis_gibbosus = lep_asv)
#plot the venn diagram
venn.diagram(
  x1,
  filename = 'Figures/venn_fish_water.png',
  category.names = c("Perca fluviatilis", "Lepomis gibbosus","Water"),
  output = TRUE ,
  imagetype="png" ,
  col = c("#A90C38",'darkorange','#21908dff'),
  fill = c("#A90C38",'darkorange','#21908dff'),
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.7,
  cat.col = c("#A90C38",'darkorange','#21908dff')
)

## CRJ1
# subsampling the asv 
water_asv_crj1 <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CRJ1.*ADN1|ADN1.*CRJ1", rownames(asv_table_filtered)), ]) != 0]
lep_asv_crj1 <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CRJ.*GAR|GAR.*CRJ", rownames(asv_table_filtered)), ]) != 0]
per_asv_crj1 <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CRJ.*PER|PER.*CRJ", rownames(asv_table_filtered)), ]) != 0]
#create a list with the subsampling asv list needed
x3 <- list(water = water_asv_crj1, perca_fluviatilis = per_asv_crj1, lepomis_gibbosus = lep_asv_crj1)
# plot the venn diagram
venn.diagram(
  x3,
  filename = 'Figures/venn_fish_water_crj.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  category.names = c("Perca fluvatilis", "Lepomis gibbosus","Water"),
  col=c('darkorange', "#440154ff",'#21908dff'),
  fill = c( alpha('darkorange',0.3),alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c('darkorange', "#440154ff",'#21908dff'),
  rotation = 1,
  main = "Cergy 1",
  main.cex = 0.5,
)

## CRJ2
# subsampling the asv 
water_asv_crj2 <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CRJ2.*ADN1|ADN1.*CRJ2", rownames(asv_table_filtered)), ]) != 0]
lep_asv_crj2 <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("TRI.*GAR|GAR.*TRI", rownames(asv_table_filtered)), ]) != 0]
per_asv_crj2 <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("TRI.*PER|PER.*TRI", rownames(asv_table_filtered)), ]) != 0]
#create a list with the subsampling asv list needed
x4 <- list(water = water_asv_crj2, perca_fluviatilis = per_asv_crj2, lepomis_gibbosus = lep_asv_crj2)
# plot the venn diagram
venn.diagram(
  x4,
  filename = 'Figures/venn_fish_water_crj2.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  category.names = c("Perca fluvatilis", "Lepomis gibbosus","Water"),
  col=c('darkorange', "#440154ff",'#21908dff'),
  fill = c( alpha('darkorange',0.3),alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c('darkorange', "#440154ff",'#21908dff'),
  rotation = 1,
  main = "Cergy 2",
  main.cex = 0.5,
)

## CTL
# subsampling the asv
water_asv_ctl <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CTL.*ADN1|ADN1.*CTL", rownames(asv_table_filtered)), ]) != 0]
lep_asv_ctl <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CTL.*GAR|GAR.*CTL", rownames(asv_table_filtered)), ]) != 0]
per_asv_ctl <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CTL.*PER|PER.*CTL", rownames(asv_table_filtered)), ]) != 0]
#create a list with the subsampling asv list needed
x5 <- list(water = water_asv_ctl, perca_fluviatilis = per_asv_ctl, lepomis_gibbosus = lep_asv_ctl)
# plot the venn diagram
venn.diagram(
  x5,
  filename = 'Figures/venn_fish_water_ctl.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  category.names = c("Perca fluvatilis", "Lepomis gibbosus","Water"),
  col=c('darkorange', "#440154ff",'#21908dff'),
  fill = c( alpha('darkorange',0.3),alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c('darkorange', "#440154ff",'#21908dff'),
  rotation = 1,
  main = "Créteil",
  main.cex = 0.5,
)

## GDP
# subsampling the asv
water_asv_gdp <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("GDP.*ADN1|ADN1.*GDP", rownames(asv_table_filtered)), ]) != 0]
lep_asv_gdp <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("GDP.*GAR|GAR.*GDP", rownames(asv_table_filtered)), ]) != 0]
per_asv_gdp <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("GDP.*PER|PER.*GDP", rownames(asv_table_filtered)), ]) != 0]

#create a list with the subsampling asv list needed
x6 <- list(water = water_asv_gdp, perca_fluviatilis = per_asv_gdp, lepomis_gibbosus = lep_asv_gdp)
# plot the venn diagram
venn.diagram(
  x6,
  filename = 'Figures/venn_fish_water_gdp.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  category.names = c("Perca fluvatilis", "Lepomis gibbosus","Water"),
  col=c('darkorange', "#440154ff",'#21908dff'),
  fill = c( alpha('darkorange',0.3),alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c('darkorange', "#440154ff",'#21908dff'),
  rotation = 1,
  main = "Grande paroisse",
  main.cex = 0.5,
)

## VER
# subsampling the asv
water_asv_ver <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("VER.*ADN1|ADN1.*VER", rownames(asv_table_filtered)), ]) != 0]
lep_asv_ver <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("VER.*GAR|GAR.*VER", rownames(asv_table_filtered)), ]) != 0]
per_asv_ver <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("VER.*PER|PER.*VER", rownames(asv_table_filtered)), ]) != 0]
#create a list with the subsampling asv list needed
x7 <- list(water = water_asv_ver, perca_fluviatilis = per_asv_ver, lepomis_gibbosus = lep_asv_ver)
# plot the venn diagram
venn.diagram(
  x7,
  filename = 'Figures/venn_fish_water_ver.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  category.names = c("Perca fluvatilis", "Lepomis gibbosus","Water"),
  col=c('darkorange', "#440154ff",'#21908dff'),
  fill = c( alpha('darkorange',0.3),alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c('darkorange', "#440154ff",'#21908dff'),
  rotation = 1,
  main = "Verneuil",
  main.cex = 0.5,
)

## CHA
# subsampling the asv
water_asv_cha <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CHA.*ADN1|ADN1.*CHA", rownames(asv_table_filtered)), ]) != 0]
lep_asv_cha <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CHA.*GAR|ADN1.*GAR", rownames(asv_table_filtered)), ]) !=0]
per_asv_cha <- colnames(asv_table_filtered)[colSums(asv_table_filtered[grep("CHA.*PER|ADN1.*PER", rownames(asv_table_filtered)), ]) != 0]

#create a list with the subsampling asv list needed
x8 <- list(water = water_asv_cha, perca_fluviatilis = per_asv_cha, lepomis_gibbosus = lep_asv_cha)
# plot the venn diagram
venn.diagram(
  x8,
  filename = 'Figures/venn_fish_water_cha.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  category.names = c("Perca fluvatilis", "Lepomis gibbosus","Water"),
  col=c('darkorange', "#440154ff",'#21908dff'),
  fill = c( alpha('darkorange',0.3),alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c('darkorange', "#440154ff",'#21908dff'),
  rotation = 1,
  main = "Champs sur Marne",
  main.cex = 0.5,
)


library(gridExtra)
library(gridGraphics)
library(png)

# Charger les images des diagrammes de Venn
venn_crj1 <- readPNG("Figures/venn_fish_water_crj.png")
venn_crj2 <- readPNG("Figures/venn_fish_water_crj2.png")
venn_ctl <- readPNG("Figures/venn_fish_water_ctl.png")
venn_gdp <- readPNG("Figures/venn_fish_water_gdp.png")
venn_ver <- readPNG("Figures/venn_fish_water_ver.png")
venn_cha <- readPNG("Figures/venn_fish_water_cha.png")

# Afficher les diagrammes de Venn sur une même page
grid.arrange(
  rasterGrob(venn_crj1, interpolate=TRUE),
  rasterGrob(venn_crj2, interpolate=TRUE),
  rasterGrob(venn_ctl, interpolate=TRUE),
  rasterGrob(venn_gdp, interpolate=TRUE),
  rasterGrob(venn_ver, interpolate=TRUE),
  rasterGrob(venn_cha, interpolate=TRUE),
  nrow = 2, ncol = 3
)

