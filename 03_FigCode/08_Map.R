####
#### µCeta project field map
#### 2024.12.12 Ushio
####

# Load library
library(tidyverse); packageVersion("tidyverse")
library(patchwork); packageVersion("patchwork")
library(ggmap); packageVersion("ggmap")


# ------------------------------------------ #
# Visualize sampling locations
# ------------------------------------------ #
# Register API (get your own API)
ggmap::register_stadiamaps("xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx")
# Draw Hong Kong Map
hk_lonlat <- c(left = 113.75, right = 114.5, bottom = 22.1, top = 22.6)
# Maptype "alidade_smooth", "outdoors", "stamen_terrain"
hk_map <- get_stadiamap(hk_lonlat, zoom = 11, maptype = "stamen_terrain")
#ggmap(hk_map)

# Point data
hk_p_df <- read.csv("data_table/data_natural_samples.csv")

# Sampling locations
g1 <- ggmap(hk_map) +
  geom_point(data = hk_p_df, aes(lat, lon), color = "red3", size = 1) +
  xlab(expression(paste("Longitude (", degree, "E)"))) +
  ylab(expression(paste("Latitude (", degree, "N)")))

# Save plot
ggsave("data_img/HK_map.jpg", g1, width = 6, height = 4)
