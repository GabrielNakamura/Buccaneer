# Time series metrics

``` r
library(ggplot2)
library(dplyr)
```

## Background

In this article we will show how to calculate different metrics of
interspecific competition (or sometimes coexistence metrics) at
different taxonomic levels in deep time using {`Buccaneer`} R package.

For all metrics we will use data from [Graciotti et al.
(2025)](https://academic.oup.com/evolut/article-abstract/79/9/1835/8140865?redirectedFrom=fulltext).
This comprises fossil record of North American canids, a group that has
been well caracterized ecologically and regarding fossil record. The
data is part of the {`Buccaneer`} package and is comprised by:

- Species longevity data: A data frame describing the Origination Time
  (TS) and Extinction Time (TE) obtained with Bayesian framework called
  [PyRate](https://gabrielnakamura.github.io/Buccaneer/articles/)

- Occurrence recorda data:

- Species traits:

## Loading data and packages

``` r
data("df_longevities_canidae") # longevities for one replicate
#data("occurrences_canidae") # occurrence data
data("traits_canidae") # trait data
data("df_occ_canidae") # occurrence data
```

### Clade level analysis

Calculating time series at clade level

``` r
res_clade_regional <-
  clade_regional_coexistence(df.TS.TE = df_longevities_canidae,
                             time.slice = 0.1,
                             round.digits = 10,
                             species = "species",
                             TS = "TS",
                             TE = "TE")
```

Plotting the results

``` r
res_clade_regional |>
  mutate(time.slice = as.numeric(time.slice)) |>
  ggplot(aes(x = time.slice, y = coexistence)) +
  geom_smooth(
    se = TRUE,
    method = "loess",
    linewidth = 0.6
  ) +
  geom_line(
    linewidth = 0.4
  ) +
  labs(
    x = "Time (Ma)",
    y = "Coexistence"
  ) +
  scale_x_reverse(
    breaks = seq(
      from = 0,
      to   = max(res_clade_regional$time.slice, na.rm = TRUE),
      by   = 5
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.title = element_text(size = 9),
    axis.text  = element_text(size = 7),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    panel.grid.minor = element_blank()
  )
```

#### regional clade trait distances

Here we calculate the mean distances using a general function that
allows to set up the group comparison in which the mean will be
calculate. It can be `between` to calculate mean distances between a
focal and a target group or `within` to calculate the mean trait
distance only within the focal group

This function is also a general solution to calculate the mean trait
distance using different number of species/genus that are close to each
other. This is set up in the argument `nearest.taxon`. 1 is equivalent
of mnnd.

Calculating comparing mpd between and within groups using a distance
matrix and a distance threshold (1 = mnnd)

``` r
dist_body_mass <- dist(traits_canidae$LD1)

regional_mnd <- 
  clade_regional_distance(df.TS.TE = df_longevities_canidae, 
                          time.slice = 0.1,
                          dist.trait = dist_body_mass,
                          nearest.taxon = 1, 
                          round.digits = 1, 
                          species = "species", 
                          TS = "TS", 
                          TE = "TE")

# distances between groups using as focal group Caniformia
regional_mnd_between <- 
  clade_regional_distance(df.TS.TE = df_longevities_canidae, 
                          time.slice = 0.1,
                          dist.trait = dist_body_mass,
                          nearest.taxon = 1, 
                          round.digits = 1, 
                          species = "species", 
                          TS = "TS", 
                          TE = "TE", 
                          group = "group", 
                          group.focal.compare = c("Caniformia", "Feliformia"), 
                          type.comparison = "between")

# distance within the focal group, in this case is Caniformia
regional_mnd_within <- 
  clade_regional_distance(df.TS.TE = df_TS_TE_mass, 
                          time.slice = 0.1,
                          dist.trait = dist_body_mass,
                          nearest.taxon = 1, 
                          round.digits = 1, 
                          species = "species", 
                          TS = "TS", 
                          TE = "TE", 
                          group = "group", 
                          group.focal.compare = c("Caniformia", "Feliformia"), 
                          type.comparison = "within")
```

Joining all results to plot in one graphic

``` r
all_mnd_regional <- rbind(regional_mnd_between, regional_mnd_within)
all_mnd_regional2 <- 
  data.frame(all_mnd_regional, 
             group_res = rep(c("Between", "Within"), 
                             each = nrow(regional_mnd_within))
  )
```

Plotting mpd results

``` r
all_mnd_regional2 |> 
  mutate(time.slice.numeric = as.numeric(unlist(lapply(strsplit(time.slice, "_"), function(x) x[2])))) |> 
  ggplot(aes(x = time.slice.numeric, y = mean.distance, color = group_res)) +
  geom_point(aes(x = time.slice.numeric, y = mean.distance)) +
  geom_smooth(se = TRUE, method = "loess", size = 1) +  # Line thickness for better visibility
  labs(title = "",
       x = "Time",
       y = "Mean Distance") +
  scale_x_continuous(breaks = seq(40, 0, by = -10)) +
  xlim(40, 0) +
  theme_minimal() +
  theme(legend.position = "bottom")# Use a clean theme
```

#### Clade site richness

``` r
res_clade_site <- 
  clade_site_coexistence(df.TS.TE = df_longevities_canidae, 
                         df.occ = df_occ_canidae, 
                         time.slice = 0.1, 
                         round.digits = 1,
                         species = "species",
                         TS = "TS",
                         TE = "TE", 
                         Max.age = "max_T",
                         Min.age = "min_T",
                         site = "site.char")
```

Plotting the results

``` r
res_clade_site |> 
  ggplot(aes(x = as.numeric(time.slice), y = mean.coexistence)) +
  geom_line(aes(x = as.numeric(time.slice), y = mean.coexistence)) +
  geom_smooth(se = TRUE, method = "loess", size = 1) +  # Line thickness for better visibility
  labs(title = "",
       x = "Time (Ma)",
       y = "Mean site coexistence") +
  scale_x_continuous(breaks = seq(max(as.numeric(res_clade_site$time.slice)), 0, by = -10)) +
  scale_x_reverse() +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.title = element_text(size = 9),
    axis.text  = element_text(size = 7),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    panel.grid.minor = element_blank()
  )
```

### Individual species analysis

In individual species analyses the focus is on the “perception” of each
individual species with other co-occurring species. There are two main
components in this module. The first computes the number of species in
which a given species co-occur. The other is the mean distance between a
focal individual and its co-occurring species.

As the other modules there are four spatial levels in which we compute
the metrics. The regional level, that calculates the number of
co-occurrence by each focal species and all other species with
occurrence data in the specified time-slice. The site level, that
compute the taxonomic and trait metrics considering as co-existing
species only those that present occurrence in the same site (usually
defined accordingly to the PBDB definition of a site, or another
definition provided by the user). The reach level, that considers a
proxy of species dispersal capacity to infer the species that
potentially co-occurr. Finally, the assemblage level, that considers as
species co-occurring only those with an occurrence record in the same
assemblage, this last being defined as a grid of customized size.

#### Individual species regional coexistence - taxonomic

``` r
res_indiv_species_coex <- 
  IndivSpecies_regional_coex(df.TS.TE = df_longevities_canidae, 
                             time.slice = 0.1,
                             round.digits = 1,
                             species = "species",
                             TS = "TS",
                             TE = "TE")


res_indiv_species_coex |> 
  group_by(species) |> 
  mutate(media = mean(n.coexistence)) |> View()
```

We can also plot the time series for each species

``` r
res_indiv_species_coex |>
  filter(
    species %in% c(
      "Cynarctus_marylandica",
      "Psalidocyon_marianae",
      "Aelurodon_taxoides",
      "Leptocyon_douglassi"
    )
  ) |>
  ggplot(
    aes(
      x = time.slice,
      y = n.coexistence,
      color = species,
      group = species
    )
  ) +
  geom_smooth(
    se = TRUE,
    method = "loess",
    linewidth = 0.8
  ) +
  geom_line(
    linewidth = 0.5
  ) +
  facet_wrap(~ species, nrow = 2) +
  labs(
    x = "Time",
    y = "Number of coexistences by species"
  ) +
  scale_x_reverse(
    breaks = seq(
      from = 0,
      to   = max(res_indiv_species_coex$time.slice, na.rm = TRUE),
      by   = 10
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10) +
    theme(
    legend.position = "none",
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.4),
    axis.ticks.x       = element_line(color = "black", linewidth = 0.4),
    axis.text.x        = element_text(size = 8),
    axis.line.y.left   = element_line(color = "black", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 20, l = 5)
  )
```

#### Individual species regional mpd

We can also plot individual species mpd for all species coexisting in
regional context. For this we will also use data on body mass for all
species. There are different ways in which mpd can be calculated. We can
calculate it using a single trait or a combination of traits. The later
require a distance trait matrix (passed to `dist` argument) containing
the pairwise distances among species.

Here we show how to calculate mpd using a single trait, in this case
body mass.

``` r

res_indiv_regional_species_mpd <-
  IndivSpec_regional_distance(df.TS.TE = df_longevities_canidae,
                            time.slice = 0.1, 
                            dist.trait = dist_body_mass,
                            round.digits = 1, 
                            species = "species",
                            TS = "TS",
                            TE = "TE", 
                            nearest.taxon = "all")
```

Plotting the results for individual species mpd

``` r
res_indiv_regional_species_mpd |>
  filter(
    species %in% c(
      "Cynarctus_marylandica",
      "Psalidocyon_marianae",
      "Aelurodon_taxoides",
      "Leptocyon_douglassi"
    )
  ) |>
  ggplot(
    aes(
      x = time.slice,
      y = mean_dist_to_cooccur,
      color = species,
      group = species
    )
  ) +
  geom_smooth(
    se = TRUE,
    method = "loess",
    linewidth = 0.8
  ) +
  geom_line(
    linewidth = 0.5
  ) +
  facet_wrap(~ species, nrow = 2) +
  labs(
    x = "Time (Ma)",
    y = "Mean pairwise distance"
  ) +
  scale_x_reverse(
    breaks = seq(
      from = 0,
      to   = max(res_indiv_regional_species_mpd$time.slice, na.rm = TRUE),
      by   = 10
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10) +
    theme(
    legend.position = "none",
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.4),
    axis.ticks.x       = element_line(color = "black", linewidth = 0.4),
    axis.text.x        = element_text(size = 8),
    axis.line.y.left   = element_line(color = "black", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 20, l = 5)
  )
```

##### Individual species regional mpd with varying distances and multiple traits

In order to calculate the mpd for multiple traits the user must inform
an object of class dist containing all pairwise distances among species.
The function can handle distances in different ways. It can use all
pairwise distances to obtain the mpd, or the user can inform throught
the argument `nearest.taxon` the number of nearest species that will be
used to calculate the mpd. The default is `all` which is the same as
calculating the mpd considering all species co-occurring in a time
slice.

``` r
res_indiv_regional_species_mnd <- 
  IndivSpec_regional_distance(df.TS.TE = df_longevities_canidae, 
                            time.slice = 0.1, 
                            dist.trait = dist_body_mass, 
                            nearest.taxon = 1,
                            round.digits = 1, 
                            species = "species",
                            TS = "TS",
                            TE = "TE")
```

#### Individual species site coexistence

Here we will perform the same calculations but using site criteria for
individual species coexistence

``` r
res_indiv_species_site_coex <- 
  IndivSpec_site_coexistence(df.TS.TE = df_longevities_canidae,
                             df.occ = df_occ_canidae, 
                             time.slice = 0.1,
                             round.digits = 5,
                             species = "species",
                             TS = "TS",
                             TE = "TE",
                             Max.age = "max_T",
                             Min.age = "min_T",
                             site = "site.char")
```

plotting species coexistence considering site coexistence

``` r
res_indiv_species_site_coex |>
  filter(
    species %in% c(
      "Cynarctus_marylandica",
      "Psalidocyon_marianae",
      "Aelurodon_taxoides",
      "Leptocyon_douglassi"
    )
  ) |>
  ggplot(
    aes(
      x = time.slice,
      y = n.coexistence,
      color = species,
      group = species
    )
  ) +
  geom_smooth(
    se = TRUE,
    method = "loess",
    linewidth = 0.8
  ) +
  geom_line(
    linewidth = 0.5
  ) +
  facet_wrap(~ species, nrow = 2) +
  labs(
    x = "Time",
    y = "Number of coexistences by species"
  ) +
  scale_x_reverse(
    breaks = seq(
      from = 0,
      to   = max(res_indiv_species_coex$time.slice, na.rm = TRUE),
      by   = 10
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10) +
    theme(
    legend.position = "none",
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.4),
    axis.ticks.x       = element_line(color = "black", linewidth = 0.4),
    axis.text.x        = element_text(size = 8),
    axis.line.y.left   = element_line(color = "black", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 20, l = 5)
  )
```

#### Individual species site mpd

We can calculate mean pairwise distances of individual species
coexistence metrics based on site co-occurrence

``` r

df_TS_TE_mass3 <- 
  df_TS_TE_mass2 |> 
  rename(species = "Genus")


indiv_species_site_mpd <- 
  IndivSpec_site_mpd(df.TS.TE = df_TS_TE_mass3, 
                     df.occ = df_occ_faurby, 
                     time.slice = 0.1, 
                     trait = "mean.size",
                     round.digits = 1,
                     species = "species",
                     TS = "TS",
                     TE = "TE",
                     Max.age = "Max.age",
                     Min.age = "Min.age",
                     site = "site")
```

#### Individual reach coexistence

We will calculate individual mean coexistence using reach criteria

``` r

indiv_species_reach_coex <- 
  IndivSpec_reach_richness(df.TS.TE = df_TS_TE_faurby,
                           df.occ = df_occ_faurby,
                           time.slice = 0.1, 
                           round.digits = 1,
                           species = "species",
                           TS = "TS",
                           TE = "TE",
                           Max.age = "Max.age",
                           Min.age = "Min.age",
                           lat = "lat",
                           lon = "lng")
```

plotting individual species coexistence with reach criteria

``` r

indiv_species_reach_coex2 <- 
  indiv_species_reach_coex |> 
  mutate(time.slice = as.numeric(time.slice))

indiv_species_reach_coex2 |> 
  filter(species == "Ursus" | species == "Simocyon" | species == "Plionarctos" | species == "Vormela" | species == "Hesperocyon") |> 
  ggplot(aes(x = time.slice, y = n.coexistence, color = species, group = species)) +  # Line thickness for better visibility
  geom_smooth(se = TRUE, method = "loess", size = 1) +
  facet_wrap(~species) +
  labs(title = "",
       x = "Time",
       y = "mean coexistence") +
  theme_minimal()  +
  scale_x_continuous(breaks = seq(max(indiv_species_site_coex2$time.slice), 0, by = -10)) +
  scale_x_reverse() +
  theme(legend.position = "none")
```

#### Individual reach coexistence mpd
