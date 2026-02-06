# Calculate Mean Trait Distance for Individual Species Co-occurring at Local Sites

This function computes mean pairwise trait distances between species
that co-occur at specific fossil sites across different time slices. For
each species present at a site, it calculates the mean distance to its
co-occurring species, supporting both mean nearest neighbor distance
(MNND) and mean pairwise distance (MPD) metrics. The function can
perform comparisons between groups or within a single group.

## Usage

``` r
IndivSpec_site_distance(
  df.TS.TE,
  df.occ,
  time.slice,
  dist.trait,
  nearest.taxon,
  group = NULL,
  group.focal.compare = NULL,
  type.comparison = NULL,
  trait = NULL,
  round.digits = 1,
  species = "species",
  TS = "TS",
  TE = "TE",
  Max.age = "Max.age",
  Min.age = "Min.age",
  site = "site"
)
```

## Arguments

- df.TS.TE:

  A data frame containing species temporal and, optionally, trait data
  with at least three columns: species names, origination times (TS),
  extinction times (TE), and, optionally, trait values. Additional
  columns may include group assignments.

- df.occ:

  A data frame containing fossil occurrence records with at least four
  columns: species names, minimum age, maximum age, and site location
  ID. Each row represents a single occurrence record at a specific site.

- time.slice:

  Numeric. The time interval (in the same units as TS and TE) between
  consecutive time slices for binning occurrences.

- dist.trait:

  A distance matrix object (class `dist` or `matrix`) containing
  pairwise trait distances between species. Row and column names must
  match species names in `df.TS.TE`. If NULL, distances will be computed
  from the trait column in `df.TS.TE` using Euclidean distance.

- nearest.taxon:

  Numeric or character. The number of nearest neighbors to consider when
  calculating mean distances. Use `1` for mean nearest neighbor distance
  (MNND), or `"all"` for mean pairwise distance (MPD).

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

- trait:

  Character. The name of the column in `df.TS.TE` containing trait
  values. If NULL (default), `dist.trait` must be provided.

- round.digits:

  Integer. The number of decimal places to round time slice values.
  Default is 1. This affects temporal binning precision.

- species:

  Character. The name of the column in `df.TS.TE` and `df.occ`
  containing species identifiers. Default is "species".

- TS:

  Character. The name of the column in `df.TS.TE` containing origination
  times for each species. Default is "TS".

- TE:

  Character. The name of the column in `df.TS.TE` containing extinction
  times for each species. Default is "TE".

- Max.age:

  Character. The name of the column in `df.occ` containing the maximum
  time estimate for each occurrence record. Default is "Max.age".

- Min.age:

  Character. The name of the column in `df.occ` containing the minimum
  age estimate for each occurrence record. Default is "Min.age".

- site:

  Character. The name of the column in `df.occ` containing site location
  identifiers. Default is "site".

## Value

A data frame with four columns:

- species:

  Character. The name of each species.

- time.slice:

  Numeric. The time slice identifier.

- mean_dist_to_cooccur:

  Numeric. The mean trait distance from each species to its co-occurring
  species at each site and time slice. Returns 0 for singleton species
  (species with no co-occurring taxa), and NA when data is insufficient.

- is_singleton:

  Logical. TRUE if the species occurred alone at the site during that
  time slice, FALSE otherwise.

## Details

The function performs the following steps:

1.  Divides the temporal range into discrete time slices based on
    `time.slice`

2.  Determines which species co-occur at each site within each time
    slice

3.  For each species at each site, calculates mean trait distances to
    co-occurring species

4.  Optionally filters comparisons by group membership (focal vs.
    comparison groups)

Species-level results include:

- **Singleton species**: Species occurring alone at a site receive a
  `mean_dist_to_cooccur` value of 0 and `is_singleton = TRUE`

- **Missing data**: Time slices or sites with insufficient data return
  NA

- **Group filtering**: When using `group.focal.compare`, only distances
  between specified groups are computed

## Examples

``` r
# Create example fossil data
df_longevities <- data.frame(
  species = c("sp1", "sp2", "sp3", "sp4"),
  TS = c(100, 98, 98, 85),
  TE = c(50, 45, 40, 35),
  trait = c(1.2, 2.5, 3.1, 4.0),
  group = c("A", "A", "B", "B")
)

df_occurrences <- data.frame(
  species = c("sp1", "sp1", "sp2", "sp3", "sp4"),
  Max.age = c(90, 95, 95, 95, 85),
  Min.age = c(60, 65, 50, 80, 80),
  site = c("site1", "site2", "site2", "site1", "site2")
)

# Calculate MPD for all species at each site
result <- IndivSpec_site_distance(
  df.TS.TE = df_longevities,
  df.occ = df_occurrences,
  trait = "trait",
  time.slice = 5,
  dist.trait = NULL,
  nearest.taxon = "all"
)

# Calculate MNND between groups at each site
result_between <- IndivSpec_site_distance(
  df.TS.TE = df_longevities,
  df.occ = df_occurrences,
  time.slice = 5,
  trait = "trait",
  dist.trait = NULL,
  nearest.taxon = "all",
  group = "group",
  group.focal.compare = c("A", "B"),
  type.comparison = "between"
)
#> Warning: argument is not numeric or logical: returning NA
#> Warning: argument is not numeric or logical: returning NA
#> Warning: argument is not numeric or logical: returning NA
#> Warning: argument is not numeric or logical: returning NA
#> Warning: argument is not numeric or logical: returning NA
#> Warning: argument is not numeric or logical: returning NA
#> Warning: NAs introduced by coercion
#> Warning: NAs introduced by coercion
#> Warning: NAs introduced by coercion
#> Warning: NAs introduced by coercion
#> Warning: NAs introduced by coercion
#> Warning: NAs introduced by coercion


result_within <- IndivSpec_site_distance(
  df.TS.TE = df_longevities,
  df.occ = df_occurrences,
  time.slice = 5,
  trait = "trait",
  dist.trait = NULL,
  nearest.taxon = "all",
  group = "group",
  group.focal.compare = c("A", "B"),
  type.comparison = "within"
)
```
