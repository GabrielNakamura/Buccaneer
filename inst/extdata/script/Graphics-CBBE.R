
df_assemblage_timeseries_richness <- df_assemblage_reg_richness$time_series_rich

df_assemblage_timeseries_richness |>
  ggplot(aes(x = time.slice, y = mean.rich.by.grid.slice)) +  # Line thickness for better visibility
  geom_point(aes(x = time.slice, y = mean.rich.by.grid.slice)) +
  geom_smooth(se = TRUE, method = "loess", size = 1) +
  labs(title = "",
       x = "Tempo",
       y = "Riqueza média") +
  theme_minimal()  +
  scale_x_continuous(breaks = seq(max(df_assemblage_timeseries_richness$time.slice), 0, by = -10)) +
  scale_x_reverse() +
  theme(legend.position = "none")# Use a clean theme



regional_richeness_clade <-
  regional_richness |>
  mutate(time.slice = as.numeric(time.slice)) |>
  ggplot(aes(x = time.slice, y = richness)) +
  geom_smooth(se = TRUE, method = "loess", size = 1) +  # Line thickness for better visibility
  geom_point(aes(x = time.slice, y = richness), shape = 21) +
  labs(title = "",
       x = "Tempo",
       y = "Riqueza") +
  scale_x_continuous(breaks = seq(max(as.numeric(regional_richness$time.slice)), 0, by = -10)) +
  scale_x_reverse() +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_line(linetype = "dashed"))# Use a clean theme


plot_mpd_regional1 <-
  regional_mpd |>
  ggplot(aes(x = time.slice, y = ses.mpd)) +
  geom_point(aes(x = time.slice, y = ses.mpd), shape = 21) +
  geom_smooth(se = TRUE, method = "loess", size = 1) +  # Line thickness for better visibility
  labs(title = "",
       x = "Tempo",
       y = "Distância média (ses mpd)") +
  scale_x_continuous(breaks = seq(40, 0, by = -10)) +
  xlim(40, 0) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_line(linetype = "dashed"))# Use a clean theme


plot_mpd_regional_all <-
  all_mpd_regional2 |>
  ggplot(aes(x = time.slice, y = ses.mpd, color = group_res)) +
  geom_point(aes(x = time.slice, y = ses.mpd), shape = 21) +
  geom_smooth(se = TRUE, method = "loess", size = 1) +  # Line thickness for better visibility
  labs(title = "",
       x = "Tempo",
       y = "Distância média (ses mpd)") +
  scale_x_continuous(breaks = seq(40, 0, by = -10)) +
  xlim(40, 0) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_line(linetype = "dashed"))# Use a clean theme


# individual species

plot_individual_species_all <-
  df_indiv_species_long |>
  filter(species == "Canis" | species == "Ursus" | species == "Ursavus") |>
  ggplot(aes(x = species, y = n.coexistence, fill = species)) +
  geom_violin(trim = FALSE) +  # Line thickness for better visibility
  geom_point(aes(x = species, y = n.coexistence, fill = species), shape = 21) +
  labs(title = "",
       x = "Espécies",
       y = "Coexistência média") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_line(linetype = "dashed"))# Use a clean theme


plot_individual_ncoexistence <-
  df_indiv_species_long |>
  filter(species == "Canis" | species == "Ursavus" | species == "Ursus" | species == "Panthera") |>
  ggplot(aes(x = time.slice, y = n.coexistence, color = species, group = species)) +  # Line thickness for better visibility
  geom_smooth(se = TRUE, method = "loess", size = 1) +
  geom_point(aes(x = time.slice, y = n.coexistence, color = species), shape = 21) +
  facet_wrap(~species) +
  labs(title = "",
       x = "Tempo",
       y = "Coexistência média por espécie") +
  theme_minimal()  +
  scale_x_continuous(breaks = seq(max(df_indiv_species_long$time.slice), 0, by = -10)) +
  scale_x_reverse() +
  theme(legend.position = "none")# Use a clean theme

df_indiv_species_mpd_long |>
  filter(species == "Canis" | species == "Ursavus" | species == "Ursus" | species == "Panthera") |>
  ggplot(aes(x = time.slice, y = mpd, color = species, group = species)) +  # Line thickness for better visibility
  geom_point(aes(x = time.slice, y = mpd, color = species), shape = 21) +
  geom_smooth(se = TRUE, method = "loess", size = 1) +
  facet_wrap(~species) +
  labs(title = "",
       x = "Tempo",
       y = "Distância média (mpd)") +
  theme_minimal()  +
  scale_x_continuous(breaks = seq(max(df_indiv_species_mpd_long$time.slice), 0, by = -10)) +
  scale_x_reverse() +
  theme(legend.position = "none",
        panel.grid.major = element_line(linetype = "dashed"))# Use a clean theme


# Assemblages

df_assemblage_timeseries_richness <- df_assemblage_reg_richness$time_series_rich

plot_assemblage_richness <-
  df_assemblage_timeseries_richness |>
  ggplot(aes(x = time.slice, y = mean.rich.by.grid.slice)) +  # Line thickness for better visibility
  geom_point(aes(x = time.slice, y = mean.rich.by.grid.slice), shape = 21) +
  geom_smooth(se = TRUE, method = "loess", size = 1) +
  labs(title = "",
       x = "Tempo",
       y = "Número médio de espécies/assembleia") +
  theme_minimal()  +
  scale_x_continuous(breaks = seq(max(df_assemblage_timeseries_richness$time.slice), 0, by = -10)) +
  scale_x_reverse() +
  theme(legend.position = "none",
        panel.grid.major = element_line(linetype = "dashed"))# Use a clean theme

plot_mpd_assemblage <-
  df_assemblage_reg_mpd2 |>
  ggplot(aes(x = time.slice, y = mean.mpd)) +  # Line thickness for better visibility
  geom_point(aes(x = time.slice, y = mean.mpd), shape = 21) +
  geom_smooth(se = TRUE, method = "loess", size = 1) +
  labs(title = "",
       x = "Tempo",
       y = "Distância média entre espécies/assembleia") +
  theme_minimal()  +
  scale_x_continuous(breaks = seq(max(df_assemblage_reg_richness$time_series_rich$time.slice), 0, by = -10)) +
  scale_x_reverse() +
  theme(legend.position = "none",
        panel.grid.major = element_line(linetype = "dashed"))# Use a clean theme


library(sf)
library(ggplot2)
library(rnaturalearth)
library(rcartocolor)

# Example data (replace with your actual data)
df_assemblage_grid_richness <- df_assemblage_reg_richness$grid_mean_age

# Load and transform continents
continents <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")
sf_mean_site <- sf::st_transform(df_assemblage_grid_richness, crs = "+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")
continents <- sf::st_transform(continents, crs = "+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")

# Define a bounding box for Eurasia in Mollweide projection coordinates
eurasia_bbox <- sf::st_as_sfc(sf::st_bbox(c(xmin = -180, xmax = 180, ymin = 0, ymax = 90), crs = sf::st_crs(sf_mean_site)))

# Crop the datasets
library(sf)

# Define Eurasia lat/lon bounding box
eurasia_bbox_lonlat <- st_as_sfc(st_bbox(c(
  xmin = -10,  # Western extent (longitude)
  xmax = 170,  # Eastern extent (longitude)
  ymin = 0,    # Southern extent (latitude)
  ymax = 80    # Northern extent (latitude)
), crs = 4326))  # WGS 84 CRS

# Transform to Mollweide projection
eurasia_bbox_moll <- st_transform(eurasia_bbox_lonlat, crs = "+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")



sf_mean_site_cropped <- sf::st_crop(sf_mean_site, eurasia_bbox_moll)
continents_cropped <- sf::st_crop(continents, eurasia_bbox_moll)

# Plot the map
map_rich_grid <-
  ggplot() +
  geom_sf(data = sf_mean_site_cropped, aes(geometry = geometry, fill = mean.age.grid),
          color = "white") +  # Your data
  geom_sf(data = continents_cropped, fill = NA, color = "black") +  # Continent boundaries
  rcartocolor::scale_fill_carto_c(name = "Riqueza", type = "quantitative", palette = "SunsetDark") +  # Optional: Use a color scale
  theme_bw() +
  labs(title = "",
       fill = "") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  )

# Variancia

map_var_grid <-
  ggplot() +
  geom_sf(data = sf_mean_site_cropped, aes(geometry = geometry, fill = var.age.grid),
          color = "white") +  # Your data
  geom_sf(data = continents_cropped, fill = NA, color = "black") +  # Continent boundaries
  rcartocolor::scale_fill_carto_c(name = "Variância da riqueza", type = "quantitative", palette = "SunsetDark") +  # Optional: Use a color scale
  theme_bw() +
  labs(title = "",
       fill = "") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  )


# plotting time slices in a map

df_assemblage_grid_richness_1 <-
  df_assemblage_reg_richness$time_series_rich |>
  filter(time.slice >= 0 & time.slice <= 10) |>
  group_by(grid_id) |>
  mutate(mean.rich.by.grid.slice.2 = mean(rich.by.grid)) |>
  distinct(grid_id, .keep_all = TRUE)
sf_mean_site_cropped1 <- sf::st_transform(sf::st_as_sf( df_assemblage_grid_richness_1), crs = "+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")

df_assemblage_grid_richness_2 <-
  df_assemblage_reg_richness$time_series_rich |>
  filter(time.slice >= 10.1 & time.slice <= 20) |>
  group_by(grid_id) |>
  mutate(mean.rich.by.grid.slice.2 = mean(rich.by.grid)) |>
  distinct(grid_id, .keep_all = TRUE)
sf_mean_site_cropped2 <- sf::st_transform(sf::st_as_sf(df_assemblage_grid_richness_2), crs = "+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")

df_assemblage_grid_richness_3 <-
  df_assemblage_reg_richness$time_series_rich |>
  filter(time.slice >= 20.1 & time.slice <= 30) |>
  group_by(grid_id) |>
  mutate(mean.rich.by.grid.slice.2 = mean(rich.by.grid)) |>
  distinct(grid_id, .keep_all = TRUE)
sf_mean_site_cropped3 <- sf::st_transform(sf::st_as_sf(df_assemblage_grid_richness_3), crs = "+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")


map_rich_grid1 <-
  ggplot() +
  geom_sf(data = sf_mean_site_cropped1, aes(geometry = geometry, fill = mean.rich.by.grid.slice.2),
          color = "white") +  # Your data
  geom_sf(data = continents_cropped, fill = NA, color = "black") +  # Continent boundaries
  rcartocolor::scale_fill_carto_c(name = "Riqueza", type = "quantitative", palette = "SunsetDark") +  # Optional: Use a color scale
  theme_bw() +
  labs(title = "",
       fill = "") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  )

map_rich_grid2 <-
  ggplot() +
  geom_sf(data = sf_mean_site_cropped2, aes(geometry = geometry, fill = mean.rich.by.grid.slice.2),
          color = "white") +  # Your data
  geom_sf(data = continents_cropped, fill = NA, color = "black") +  # Continent boundaries
  rcartocolor::scale_fill_carto_c(name = "Riqueza", type = "quantitative", palette = "SunsetDark") +  # Optional: Use a color scale
  theme_bw() +
  labs(title = "",
       fill = "") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  )

map_rich_grid3 <-
  ggplot() +
  geom_sf(data = sf_mean_site_cropped3, aes(geometry = geometry, fill = mean.rich.by.grid.slice.2),
          color = "white") +  # Your data
  geom_sf(data = continents_cropped, fill = NA, color = "black") +  # Continent boundaries
  rcartocolor::scale_fill_carto_c(name = "Riqueza", type = "quantitative", palette = "SunsetDark") +  # Optional: Use a color scale
  theme_bw() +
  labs(title = "",
       fill = "") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  )


# mpd in a map

df_assemblage_reg_mpd2 <-
  df_assemblage_reg_mpd$mean_mpd_grid
sf_mean_mpd_site <- sf::st_as_sf(df_assemblage_reg_mpd2)


continents <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
sf_mean_mpd_site2 <- sf::st_transform(sf_mean_mpd_site, crs = "+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")
continents <- sf::st_transform(continents, crs = "+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")
sf_mean_mpd_site_cropped <- sf::st_crop(sf_mean_mpd_site, eurasia_bbox_moll)
continents_cropped <- sf::st_crop(continents, eurasia_bbox_moll)


ggplot() +
  geom_sf(data = sf_mean_mpd_site2, aes(geometry = geometry, fill = grid.mean),
          color = "white") +  # Your data
  geom_sf(data = continents_cropped, fill = NA, color = "black") +  # Continent boundaries
  rcartocolor::scale_fill_carto_c(name = "MPD", type = "quantitative", palette = "SunsetDark") +  # Optional: Use a color scale
  theme_bw() +
  labs(title = "",
       fill = "Mean mpd by site") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom"
  )


# Example table

# Load necessary libraries
library(knitr)
library(kableExtra)

# Example dataset
data <- data.frame(
  Name = c("Alice", "Bob", "Charlie"),
  Age = c(25, 30, 35),
  Score = c(90, 80, 85)
)

# Generate the table
kable(df_assemblage_grid_richness[400:150, c(1, 2, 4)],
      col.names = c("Nível biológico", "Descritor biológico", "Coordenadas"),
      caption = "",
      format = "html") %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive"),
    full_width = FALSE,
    position = "center"
  ) %>%
  row_spec(0, bold = TRUE, background = "#D3D3D3") %>%  # Style the header
  row_spec(1:3, color = "black")  # Add row-specific styles


ggsave(here::here("output", "figures", "Fig4_models.png"),
       width = 6, height = 8,
       dpi = 600, plot = fig_model)

ggsave(here::here("output", "figures", "Fig4_models.png"),
       width = 6, height = 8,
       dpi = 600, plot = fig_model)
