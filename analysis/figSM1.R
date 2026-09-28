#### Map without API or tiles: ggplot2 + sf + rnaturalearth + ggrepel ####

library(ggplot2)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(dplyr)
library(ggrepel)

# Data
points <- data.frame(
  nom = c("VER", "CRJ1", "GDP", "CHA", "CTL", "CRJ2"),
  lat = c(48.9922, 49.0308, 48.3733, 48.8633, 48.775, 49.0261),
  long = c(1.9664, 2.0556, 2.8992, 2.5967, 2.4497, 2.0492),
  trophie = c("hypereutrophe", "eutrophe", "eutrophe",
              "hypereutrophe", "eutrophe", "eutrophe")
)

couleurs <- c(
  "CHA"  = "#00332A",
  "CRJ1" = "#F2F26D",
  "CRJ2" = "#C5D163",
  "CTL"  = "#60A360",
  "GDP"  = "#035236",
  "VER"  = "#227548"
)

names_lakes <- c(
  "CHA"  = "Champs-sur-Marne (CSL)",
  "CRJ1" = "Cergy Large (CERL)",
  "CRJ2" = "Cergy Small (CERS)",
  "CTL"  = "Créteil (CTL)",
  "GDP"  = "La Grande-Paroisse (LGP)",
  "VER"  = "Verneuil-sur-Seine (VSS)"
)

# Lakes (sf) with labels
points_sf <- st_as_sf(points, coords = c("long", "lat"), crs = 4326) |>
  mutate(label = names_lakes[nom])

# Paris (city centre)
paris_sf <- st_as_sf(
  data.frame(label = "Paris", long = 2.3522, lat = 48.8566),
  coords = c("long", "lat"), crs = 4326
)

# All labels together so ggrepel avoids overlaps between them
labels_sf <- bind_rows(
  points_sf |> select(label),
  paris_sf  |> select(label)
) |>
  mutate(is_paris = label == "Paris")

# Department outlines (local data, no API)
departements <- ne_states(country = "France", returnclass = "sf")

idf <- departements |>
  filter(name %in% c("Paris", "Seine-et-Marne", "Yvelines", "Essonne",
                     "Hauts-de-Seine", "Seine-Saint-Denis",
                     "Val-de-Marne", "Val-d'Oise"))

# Map
map_idf <- ggplot() +
  geom_sf(data = departements, fill = "grey95", color = "grey80", linewidth = 0.2) +
  geom_sf(data = idf, fill = "white", color = "grey50", linewidth = 0.3) +
  geom_sf(data = points_sf, aes(color = nom), size = 4, alpha = 1) +
  geom_sf(data = paris_sf, shape = 8, size = 4, color = "black", stroke = 1.2) +
  geom_text_repel(
    data = labels_sf,
    aes(label = label, geometry = geometry,
        fontface = ifelse(is_paris, "bold", "plain")),
    stat = "sf_coordinates",
    size = 3.5,
    min.segment.length = 0,
    segment.color = "grey40",
    box.padding = 0.6,
    point.padding = 0.4,
    seed = 1
  ) +
  scale_color_manual(values = couleurs) +
  annotation_scale(location = "br", width_hint = 0.2) +
  coord_sf(xlim = c(1.4, 3.6), ylim = c(48.1, 49.3), expand = FALSE) +
  theme_void() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "#EAF1F5", color = NA))

map_idf

ggsave(here::here("figures", "figSM1.png"),
       plot = map_idf, width = 7, height = 6, dpi = 300)
