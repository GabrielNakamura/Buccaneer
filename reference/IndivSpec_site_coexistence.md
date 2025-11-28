# Calculate Species-Level Co-occurrence at Sites Across Time Slices

This function computes the number of co-occurring species for each
individual species at fossil sites across different time slices. For
each species present at sites during a time slice, it counts how many
other species it co-occurs with, providing species-specific
co-occurrence patterns through time. The function can perform
comparisons between taxonomic groups or within a single group.

## Usage

``` r
IndivSpec_site_coexistence(
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

  Logical. Should singleton species (species occurring alone at sites
  with no co-occurring species) be excluded from the output? Default is
  TRUE. When TRUE, species with `n.coexistence = 0` are removed.

- group:

  Character. The name of the column in `df.TS.TE` containing group
  assignments for species (e.g., clade, family). Required if using
  `group.focal.compare`. Default is NULL.

- group.focal.compare:

  Character vector of length 2. The first element specifies the focal
  group and the second specifies the comparison group. If NULL
  (default), co-occurrence is calculated across all species regardless
  of group membership.

- type.comparison:

  Character. Specifies the type of co-occurrence comparison:

  - `"between"`: Count only co-occurrences between species from the
    focal and comparison groups.

  - `"within"`: Count only co-occurrences among species within the focal
    group.

  - NULL (default): Count all co-occurrences regardless of group.

## Value

A data frame with three columns:

- time.slice:

  Character. The time slice identifier (e.g., "100" for the time slice
  at 100 Ma).

- species:

  Character. The name of each species.

- n.coexistence:

  Integer. The number of other species that co-occur with the focal
  species at sites during that time slice. Self-coexistence is excluded
  (a species does not count as co-occurring with itself).

## Details

The function performs the following steps:

1.  Creates time slices from maximum TS to minimum TE

2.  Generates regional co-occurrence matrices using
    [`aux_matrix_regional_coex()`](https://gabrielnakamura.github.io/Buccaneer/reference/aux_matrix_regional_coex.md)

3.  Determines which species occur at which sites in each time slice

4.  Creates site-based co-occurrence matrices using
    [`comp_site_cooccurr()`](https://gabrielnakamura.github.io/Buccaneer/reference/comp_site_cooccurr.md)

5.  For each species, counts the number of co-occurring species across
    all sites

6.  Optionally removes singleton species (those with zero
    co-occurrences)

7.  Optionally filters by group membership

Co-occurrence calculation details:

- **Self-coexistence**: Automatically excluded from counts

- **Singleton species**: Species with `n.coexistence = 0`

  - If `remove.singletons = TRUE`: Excluded from output

  - If `remove.singletons = FALSE`: Included with value of 0

- **Site aggregation**: Co-occurrence counts are summed across all sites
  where a species occurs in each time slice

- **Group comparisons**: When using `group.focal.compare`, only counts
  co-occurrences between/within specified groups

This function differs from
[`clade_site_coexistence()`](https://gabrielnakamura.github.io/Buccaneer/reference/clade_site_coexistence.md)
by returning species-level results rather than time slice-level
aggregated means. It is useful for analyzing individual species'
ecological associations through time or identifying species with
consistently high or low co-occurrence patterns.

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

# Calculate species-level co-occurrence through time
result <- IndivSpec_site_coexistence(
  df.TS.TE = df_temporal,
  df.occ = df_occurrences,
  time.slice = 10
)

# View results
head(result)

# Plot co-occurrence patterns for a specific species
library(ggplot2)
sp1_data <- subset(result, species == "sp1")
ggplot(sp1_data, aes(x = time.slice, y = n.coexistence)) +
  geom_col() +
  labs(x = "Time Slice", y = "Number of Co-occurring Species",
       title = "Co-occurrence Pattern for sp1") +
  theme_minimal()

# Compare co-occurrence between species
ggplot(result, aes(x = time.slice, y = n.coexistence, fill = species)) +
  geom_col(position = "dodge") +
  labs(x = "Time Slice", y = "Number of Co-occurring Species") +
  theme_minimal()

# Include singleton species in results
result_with_singletons <- IndivSpec_site_coexistence(
  df.TS.TE = df_temporal,
  df.occ = df_occurrences,
  time.slice = 10,
  remove.singletons = FALSE
)

# Identify singleton species
singletons <- subset(result_with_singletons, n.coexistence == 0)
singletons

# Calculate co-occurrence between groups
result_between <- IndivSpec_site_coexistence(
  df.TS.TE = df_temporal,
  df.occ = df_occurrences,
  time.slice = 10,
  group = "group",
  group.focal.compare = c("A", "B"),
  type.comparison = "between"
)

# Calculate co-occurrence within a single group
result_within <- IndivSpec_site_coexistence(
  df.TS.TE = df_temporal,
  df.occ = df_occurrences,
  time.slice = 10,
  group = "group",
  group.focal.compare = c("A", "B"),
  type.comparison = "within"
)

# Summary statistics: mean co-occurrence per species across time
library(dplyr)
species_summary <- result %>%
  group_by(species) %>%
  summarise(
    mean_coexistence = mean(n.coexistence),
    sd_coexistence = sd(n.coexistence),
    max_coexistence = max(n.coexistence)
  )
} # }
```
