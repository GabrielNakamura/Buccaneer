# Calculate Mean Site-Level Species Coexistence Across Time Slices

This function computes the mean number of co-occurring species at
individual sites across different time slices. For each time slice, it
determines which species co-occur at each site based on occurrence
records and temporal ranges, then calculates the mean and variance of
coexistence across all species for time slice.

## Usage

``` r
clade_site_coexistence(
  df.TS.TE,
  df.occ,
  time.slice,
  round.digits = 1,
  species = "species",
  TS = "TS",
  TE = "TE",
  Max.age = "Max.age",
  Min.age = "Min.age",
  site = "site",
  remove.singletons = TRUE,
  group = NULL,
  group.focal.compare = NULL,
  type.comparison = NULL
)
```

## Arguments

- df.TS.TE:

  A data frame containing species temporal data with at least three
  columns: species names, origination times (TS), and extinction times
  (TE). Additional columns may include group assignments.

- df.occ:

  A data frame containing fossil occurrence records with at least four
  columns: species names, minimum age, maximum age, and site location
  ID. Each row represents a single occurrence record at a specific site.

- time.slice:

  Numeric. The time interval (in the same units as TS and TE) between
  consecutive time slices for temporal binning.

- round.digits:

  Integer. The number of decimal places to round time slice values.
  Default is 1. This affects temporal binning precision.

- species:

  Character. The name of the column in `df.TS.TE` and `df.occ`
  containing species identifiers. Default is "species".

- TS:

  Character. The name of the column in `df.TS.TE` containing origination
  (first appearance) times for each species. Default is "TS".

- TE:

  Character. The name of the column in `df.TS.TE` containing extinction
  (last appearance) times for each species. Default is "TE".

- Max.age:

  Character. The name of the column in `df.occ` containing the maximum
  (oldest) age estimate for each occurrence record. Default is
  "Max.age".

- Min.age:

  Character. The name of the column in `df.occ` containing the minimum
  (youngest) age estimate for each occurrence record. Default is
  "Min.age".

- site:

  Character. The name of the column in `df.occ` containing site location
  identifiers. Default is "site".

- remove.singletons:

  Logical. Should singleton species (species occurring alone at a site
  with no co-occurring species) be excluded from mean and variance
  calculations? Default is TRUE. When TRUE, singletons are treated as
  NA; when FALSE, they contribute 0 to the mean.

- group:

  Character. The name of the column in `df.TS.TE` containing group
  assignments for species (e.g., clade, family). Required if using
  `group.focal.compare`. Default is NULL.

- group.focal.compare:

  Character vector of length 2. The first element specifies the focal
  group and the second specifies the comparison group. If NULL
  (default), coexistence is calculated across all species regardless of
  group membership.

- type.comparison:

  Character. Specifies the type of coexistence comparison:

  - `"between"`: Count only co-occurrences between species from the
    focal and comparison groups.

  - `"within"`: Count only co-occurrences among species within the focal
    group.

  - NULL (default): Count all co-occurrences regardless of group.

## Value

A data frame with three columns:

- mean.coexistence:

  Numeric. The mean number of co-occurring species in each time slice.
  This represents average local coexistence (excluding or including
  singletons based on `remove.singletons`).

- var.distance:

  Numeric. The variance in the number of co-occurring species in each
  time slice.

- time.slice:

  Numeric. The time point representing each slice, typically the upper
  (older) boundary of the time bin.

## Details

The function performs the following steps:

1.  Creates time slices from maximum TS to minimum TE

2.  Generates regional co-occurrence matrices using
    [`aux_matrix_regional_coex()`](https://gabrielnakamura.github.io/Buccaneer/reference/aux_matrix_regional_coex.md)

3.  Determines which species occur at which sites in each time slice

4.  Creates site-based co-occurrence matrices using
    [`comp_site_cooccurr()`](https://gabrielnakamura.github.io/Buccaneer/reference/comp_site_cooccurr.md)

5.  For each species at each site, counts the number of co-occurring
    species

6.  Calculates mean and variance of coexistence across species with
    cooccurrence in sites

7.  Optionally filters by group membership

Coexistence calculation details:

- **Self-coexistence**: Removed from counts (a species does not
  "co-occur" with itself)

- **Singleton species**: Species with 0 co-occurring taxa

  - If `remove.singletons = TRUE`: Treated as NA (excluded from mean)

  - If `remove.singletons = FALSE`: Contribute 0 to the mean

- **Group comparisons**: When using `group.focal.compare`, only counts
  co-occurrences between/within specified groups

This function calculates site-level (local) coexistence patterns, which
differ from regional-scale richness. For regional richness, see
[`clade_regional_coexistence()`](https://gabrielnakamura.github.io/Buccaneer/reference/clade_regional_coexistence.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Create example fossil data
df_temporal <- data.frame(
  species = c("sp1", "sp2", "sp3", "sp4"),
  TS = c(100, 95, 90, 85),
  TE = c(50, 45, 40, 35),
  group = c("A", "A", "B", "B")
)

df_occurrences <- data.frame(
  species = c("sp1", "sp1", "sp2", "sp3", "sp4", "sp4"),
  Max.age = c(100, 95, 95, 90, 85, 85),
  Min.age = c(90, 85, 85, 80, 75, 75),
  site = c("site1", "site2", "site1", "site1", "site2", "site3")
)

# Calculate mean site coexistence through time
result <- clade_site_coexistence(
  df.TS.TE = df_temporal,
  df.occ = df_occurrences,
  time.slice = 10
)

# View results
head(result)

# Plot mean coexistence through time
plot(result$time.slice,
     result$mean.coexistence,
     type = "l",
     xlab = "Time (Ma)",
     ylab = "Mean Site Richness",
     main = "Mean Species Coexistence at Sites Through Time")

# Add variance as error bars
arrows(result$time.slice,
       result$mean.coexistence - sqrt(result$var.distance),
       result$time.slice,
       result$mean.coexistence + sqrt(result$var.distance),
       length = 0.05, angle = 90, code = 3)

# Calculate coexistence including singletons
result_with_singletons <- clade_site_coexistence(
  df.TS.TE = df_temporal,
  df.occ = df_occurrences,
  time.slice = 10,
  remove.singletons = FALSE
)

# Calculate coexistence between groups
result_between <- clade_site_coexistence(
  df.TS.TE = df_temporal,
  df.occ = df_occurrences,
  time.slice = 10,
  group = "group",
  group.focal.compare = c("A", "B"),
  type.comparison = "between"
)

# Calculate coexistence within a single group
result_within <- clade_site_coexistence(
  df.TS.TE = df_temporal,
  df.occ = df_occurrences,
  time.slice = 10,
  group = "group",
  group.focal.compare = c("A", "B"),
  type.comparison = "within"
)
} # }
```
