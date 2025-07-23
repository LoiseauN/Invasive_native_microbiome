#### Create a map with the package Leaflet ####


# Create the data frame
points <- data.frame(
  nom = c("VER_W1", "VER_W2", "VER_W3", "CRJ_W1", "CRJ_W2", "CRJ_W3", 
          "GDP_W1", "GDP_W2", "GDP_W3",
          "CHA_W1", "CHA_W2", "CHA_W3", "CTL_W1", "CTL_W2", "CTL_W3", "CRJ1_W1", "CRJ1_W2", "CRJ1_W3"),
  lat = c(48.9922, 48.9933, 48.9933, 49.0308, 49.0317, 49.0296, 48.3733, 48.3719, 48.3753,
          48.8633, 48.8631, 48.8628, 48.775, 48.7744, 48.7772, 49.0261, 49.0247, 49.0247),
  long = c(1.9664, 1.9683, 1.9706, 2.0556, 2.0528, 2.0533, 2.8992, 2.8994, 2.8983,
           2.5967, 2.5975, 2.6003, 2.4497, 2.45, 2.4514, 2.0492, 2.0461, 2.0444),
  trophie = c("hypereutrophe", "hypereutrophe", "hypereutrophe", "eutrophe", "eutrophe", "eutrophe", 
              "eutrophe", "eutrophe", "eutrophe", "hypereutrophe", "hypereutrophe", "hypereutrophe",
              "eutrophe", "eutrophe", "eutrophe",
              "eutrophe", "eutrophe", "eutrophe")
)

# Define colors for each point name
couleurs <- c(
  "CHA_W1" = "#608F3D", "CHA_W2" = "#608F3D", "CHA_W3" = "#608F3D",
  "CRJ_W1" = "#41AEBD","CRJ_W2" = "#41AEBD", "CRJ_W3" = "#41AEBD",
  "CRJ1_W1" = "#97E9D5", "CRJ1_W2" = "#97E9D5",  "CRJ1_W3" = "#97E9D5",
  "CTL_W1" = "#F4DE3A", "CTL_W2" = "#F4DE3A", "CTL_W3" = "#F4DE3A",
  "GDP_W1" = "#A2CF49", "GDP_W2" = "#A2CF49", "GDP_W3" = "#A2CF49",
  "VER_W1" = "#FCB11C", "VER_W2" = "#FCB11C", "VER_W3" = "#FCB11C"
)

# Create a color palette
palette_couleurs <- colorFactor(
  palette = couleurs,
  domain = names(couleurs)
)

# Create the map with a clean background
map_idf <- leaflet() %>%
  addTiles(
    urlTemplate = "https://cartodb-basemaps-{s}.global.ssl.fastly.net/light_all/{z}/{x}/{y}.png",
    options = providerTileOptions(noWrap = TRUE, minZoom = 9, maxZoom = 12, continuousWorld = TRUE),
    attribution = NULL
  ) %>%
  setView(lng = 2.35, lat = 48.85, zoom = 10) %>%  # Center the map on Paris
  addCircleMarkers(
    data = points, 
    lng = ~long, 
    lat = ~lat, 
    color = ~palette_couleurs(nom), 
    stroke = FALSE, 
    fillOpacity = 0.8
  ) %>%
  addScaleBar(
    position = "bottomright",   # Positioning the scale bar at the bottom right
    options = scaleBarOptions(
      maxWidth = 100,           # Maximum width of the scale bar
      metric = TRUE,            # Display in kilometers
      imperial = FALSE,         # Do not display in miles
      updateWhenIdle = TRUE     # Update the scale bar when the map view is idle
    )
  )# Add colored markers based on point names

# Display the map
map_idf
