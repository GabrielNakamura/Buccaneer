# Calculate Regional Mean Trait Distances Across Time Slices

This function computes mean pairwise trait distances between species at
the regional scale across different time slices. For each time slice, it
calculates the mean distance among all species whose temporal ranges
overlap with that interval, supporting both mean nearest neighbor
distance (MNND) and mean pairwise distance (MPD) metrics. The function
can perform comparisons between groups (e.g., clades, families) or
within a single group.

## Usage

``` r
clade_regional_distance(
  df.TS.TE,
  time.slice,
  dist.trait,
  nearest.taxon,
  trait = NULL,
  round.digits = 1,
  species = "species",
  TS = "TS",
  TE = "TE",
  group = NULL,
  group.focal.compare = NULL,
  type.comparison = NULL
)
```

## Arguments

- df.TS.TE:

  A data frame containing species temporal and trait data with at least
  three columns: species names, origination times (TS), extinction times
  (TE). Trait values are optional in this object. Additional columns may
  include group assignments.

- time.slice:

  Numeric. The time interval (in the same units as TS and TE) between
  consecutive time slices for temporal binning.

- dist.trait:

  A distance matrix object (class `dist` or `matrix`) containing
  pairwise trait distances between species. Row and column names must
  match species names in `df.TS.TE`. If NULL, distances will be computed
  from the trait column using Euclidean distance.

- nearest.taxon:

  Numeric or character. The number of nearest neighbors to consider when
  calculating mean distances. Use `1` for mean nearest neighbor distance
  (MNND), or `"all"` for mean pairwise distance (MPD).

- trait:

  Character. The name of the column in `df.TS.TE` containing trait
  values. If NULL (default), `dist.trait` must be provided.

- round.digits:

  Integer. The number of decimal places to round time slice values.
  Default is 1. This affects temporal binning precision.

- species:

  Character. The name of the column in `df.TS.TE` containing species
  identifiers. Default is "species".

- TS:

  Character. The name of the column in `df.TS.TE` containing origination
  (first appearance) times for each species. Default is "TS".

- TE:

  Character. The name of the column in `df.TS.TE` containing extinction
  (last appearance) times for each species. Default is "TE".

- group:

  Character. The name of the column in `df.TS.TE` containing group
  assignments for species (e.g., clade, family). Required if using
  `group.focal.compare`. Default is NULL.

- group.focal.compare:

  Character vector of length 2. The first element specifies the focal
  group and the second specifies the comparison group. If NULL
  (default), distances are calculated across all species regardless of
  group membership.

- type.comparison:

  Character. Specifies the type of distance comparison:

  - `"between"`: Calculate distances only between species from the focal
    and comparison groups.

  - `"within"`: Calculate distances only among species within the focal
    group.

  - NULL (default): Calculate distances among all species together.

## Value

A data frame with three columns:

- mean.distance:

  Numeric. The mean trait distance among species present in each time
  slice. Returns NA when insufficient data is available (e.g., only one
  species present, or no representatives from required groups).

- var.distance:

  Numeric. The variance of trait distances among species in each time
  slice. Returns NA when insufficient data is available.

- time.slice:

  Numeric. The time point representing each slice, typically the upper
  (older) boundary of the time bin.

## Details

The function performs the following steps:

1.  Creates a sequence of time slices from maximum TS to minimum TE

2.  Generates regional co-occurrence matrices using
    [`aux_matrix_regional_coex()`](https://gabrielnakamura.github.io/Buccaneer/reference/aux_matrix_regional_coex.md)

3.  For each time slice, identifies species whose temporal ranges
    overlap

4.  Computes pairwise trait distances among overlapping species

5.  Optionally filters comparisons by group membership

6.  Calculates mean and variance of distances based on `nearest.taxon`

Distance calculation options:

- **MPD (nearest.taxon = "all")**: Calculates mean pairwise distance
  considering all pairwise comparisons among species

- **MNND (nearest.taxon = 1)**: Calculates mean nearest neighbor
  distance using only the closest species for each focal species

- **Threshold (nearest.taxon = n)**: Uses the `n` nearest neighbors for
  each focal species

Missing values (NA) are returned for time slices where:

- Only one species is present

- No species from required groups are present (when using group
  comparisons)

- Insufficient data for distance calculations

This function calculates distances at the regional (pool) scale,
considering all species present during each time slice regardless of
their geographic distribution. For site-level distance calculations, see
[`IndivSpec_site_distance()`](https://gabrielnakamura.github.io/Buccaneer/reference/IndivSpec_site_distance.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Create example fossil data with traits
df_temporal <- data.frame(
  species = c("sp1", "sp2", "sp3", "sp4"),
  TS = c(100, 95, 90, 85),
  TE = c(50, 45, 40, 35),
  trait = c(1.2, 2.5, 3.1, 4.0),
  group = c("A", "A", "B", "B")
)

# Calculate regional MPD through time
result_mpd <- clade_regional_distance(
  df.TS.TE = df_temporal,
  time.slice = 10,
  dist.trait = NULL,
  nearest.taxon = "all",
  trait = "trait"
)

# View results
head(result_mpd)

# Plot mean distance through time
plot(result_mpd$time.slice,
     result_mpd$mean.distance,
     type = "l",
     xlab = "Time (Ma)",
     ylab = "Mean Trait Distance",
     main = "Regional Trait Distance Through Time")

# Calculate MNND between groups
result_between <- clade_regional_distance(
  df.TS.TE = df_temporal,
  time.slice = 10,
  dist.trait = NULL,
  nearest.taxon = 1,
  trait = "trait",
  group = "group",
  group.focal.compare = c("A", "B"),
  type.comparison = "between"
)

# Calculate distances within a single group
result_within <- clade_regional_distance(
  df.TS.TE = df_temporal,
  time.slice = 10,
  dist.trait = NULL,
  nearest.taxon = "all",
  trait = "trait",
  group = "group",
  group.focal.compare = c("A", "B"),
  type.comparison = "within"
)

# Using a pre-computed distance matrix
dist_matrix <- dist(df_temporal$trait)
result_custom_dist <- clade_regional_distance(
  df.TS.TE = df_temporal,
  time.slice = 10,
  dist.trait = dist_matrix,
  nearest.taxon = "all"
)
} # }
```
