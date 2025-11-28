# Building time series based on species coexistence

``` r
#library(buccaneer)
library(ggplot2)
library(dplyr)
data("df_TS_TE_faurby")
data("df_occ_faurby")
data("df_TS_TE_mass")
```

## Assemblage regional richness

``` r
df_assemblage_reg_richness <- 
  Assemblage_regional_richness(df.TS.TE = df_TS_TE_faurby,
                               df.occ = df_occ_faurby, 
                               time.slice = 0.1, 
                               grid.size = 5,
                               round.digits = 1,
                               species = "species",
                               TS = "TS",
                               TE = "TE",
                               Max.age = "Max.age",
                               Min.age = "Min.age", 
                               lat = "lat", 
                               lon = "lng", 
                               crs = 4326)
```

plotting time series based on mean assemblage regional richness

``` r

df_assemblage_timeseries_richness <- df_assemblage_reg_richness$time_series_rich

df_assemblage_timeseries_richness |> 
  ggplot(aes(x = time.slice, y = mean.rich.by.grid.slice)) +  # Line thickness for better visibility
  geom_point(aes(x = time.slice, y = mean.rich.by.grid.slice)) +
  geom_smooth(se = TRUE, method = "loess", size = 1) +
  labs(title = "",
       x = "Time",
       y = "Mean species richness assemblage") +
  theme_minimal()  +
  scale_x_continuous(breaks = seq(max(df_assemblage_timeseries_richness$time.slice), 0, by = -10)) +
  scale_x_reverse() +
  theme(legend.position = "none")# Use a clean theme
```

We can plot grid values at different time slices or a summary metric
considering all time slices

``` r
df_assemblage_grid_richness <- df_assemblage_reg_richness$grid_mean_age

continents <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")
sf_mean_site <- sf::st_transform(df_assemblage_grid_richness, crs = "+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")
continents <- sf::st_transform(continents, crs = "+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")

map_rich_grid <- 
  ggplot() +
  geom_sf(data = sf_mean_site, aes(geometry = geometry, fill = mean.age.grid),
          color = "white") +  # Your data
  geom_sf(data = continents, fill = NA, color = "black") +  # Continent boundaries
  rcartocolor::scale_fill_carto_c(name = "Richness", type = "quantitative", palette = "SunsetDark") +  # Optional: Use a color scale
  theme_bw() +
  labs(title = "",
       fill = "Mean mpd by site") +
  theme(
    axis.text = element_blank(),  
    axis.ticks = element_blank(),
    axis.line = element_blank()
  ) 
```

We can also plot the variance in richness through time

``` r
map_var_grid <- 
  ggplot() +
  geom_sf(data = sf_mean_site, aes(geometry = geometry, fill = var.age.grid),
          color = "white") +  # Your data
  geom_sf(data = continents, fill = NA, color = "black") +  # Continent boundaries
  rcartocolor::scale_fill_carto_c(name = "Richness variance", type = "quantitative", palette = "SunsetDark") +  # Optional: Use a color scale
  theme_bw() +
  labs(title = "",
       fill = "richness variance") +
  theme(
    axis.text = element_blank(),  
    axis.ticks = element_blank(),
    axis.line = element_blank()
  ) 
```

Plotting together

``` r
library(patchwork)
map_grid_rich_var <- 
  map_rich_grid + map_var_grid +
  plot_layout(nrow = 2)
map_grid_rich_var
```

## Assemblage regional mpd

``` r

df_assemblage_reg_mpd <- 
  Assemblage_regional_mpd(df.TS.TE = df_TS_TE_mass,
                        df.occ = df_occ_faurby, 
                        time.slice = 0.1, 
                        grid.size = 5, 
                        trait = "mean.size", 
                        round.digits = 1,
                        species = "species", 
                        TS = "TS", 
                        TE = "TE",
                        Max.age = "Max.age",
                        Min.age = "Min.age", 
                        lat = "lat",
                        lon = "lng",
                        crs = 4326)
```

plotting mpd assemblage

``` r
df_assemblage_reg_mpd2 <- 
  df_assemblage_reg_mpd$mean_mpd_timeslice

df_assemblage_reg_mpd2 |> 
  ggplot(aes(x = time.slice, y = mean.mpd)) +  # Line thickness for better visibility
  geom_point(aes(x = time.slice, y = mean.mpd)) +
  geom_smooth(se = TRUE, method = "loess", size = 1) +
  labs(title = "",
       x = "Time",
       y = "Mean mpd assemblage") +
  theme_minimal()  +
  scale_x_continuous(breaks = seq(max(df_assemblage_reg_richness$time_series_rich$time.slice), 0, by = -10)) +
  scale_x_reverse() +
  theme(legend.position = "none") # Use a clean theme
```

## Mapping diversity in space

Here we will map richness and mpd in space through time

``` r


# plotting in a map

sf_mean_site <- sf::st_as_sf(df_assemblage_reg_mpd$mean_mpd_grid) 

continents <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
sf_mean_site2 <- sf::st_transform(sf_mean_site, crs = "+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")
continents <- sf::st_transform(continents, crs = "+proj=moll +x_0=0 +y_0=0 +lat_0=0 +lon_0=0")


ggplot() +
  geom_sf(data = sf_mean_site2, aes(geometry = geometry, fill = grid.mean),
          color = "white") +  # Your data
  geom_sf(data = continents, fill = NA, color = "black") +  # Continent boundaries
  rcartocolor::scale_fill_carto_c(name = "MPD", type = "quantitative", palette = "SunsetDark") +  # Optional: Use a color scale
  theme_bw() +
  labs(title = "",
       fill = "Mean mpd by site") +
  theme(
    axis.text = element_blank(),  
    axis.ticks = element_blank(),
    axis.line = element_blank()
  ) 
```

## Clade
