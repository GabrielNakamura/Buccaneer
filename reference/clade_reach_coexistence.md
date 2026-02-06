# Compute species coexistence for each time slice based on reach metric

Compute species coexistence for each time slice based on reach metric

## Usage

``` r
clade_reach_coexistence(
  df.TS.TE,
  df.occ,
  time.slice,
  round.digits,
  species = "species",
  TS = "TS",
  TE = "TE",
  lat = "lat",
  lon = "lon",
  Max.age = "Max.age",
  Min.age = "Min.age",
  crs = 4326,
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

- lat:

  Numeric. The latitude coordinate of the occurrence record in `df.occ`.

- lon:

  Numeric. The longitude coordinate of the occurrence records in
  `df.occ`.

- Max.age:

  Character. The name of the column in `df.occ` containing the maximum
  (oldest) age estimate for each occurrence record. Default is
  "Max.age".

- Min.age:

  Character. The name of the column in `df.occ` containing the minimum
  (youngest) age estimate for each occurrence record. Default is
  "Min.age".

- crs:

  Numeric. The code indicating the coordinate reference system to be
  used for latitude and longitude of occurrence records in `df.occ`

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

- var.coexistence:

  Numeric. The variance in the number of co-occurring species in each
  time slice.

- time.slice:

  Numeric. The time point representing each slice, typically the upper
  (older) boundary of the time bin.

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

# Calculate mean site reach through time
result <- clade_reach_coexistence(
  df.TS.TE = df_temporal,
  df.occ = df_occurrences,
  time.slice = 10
)

# View results
head(result)

# Calculate reach between groups
result_between <- clade_reach_coexistence(
  df.TS.TE = df_temporal,
  df.occ = df_occurrences,
  time.slice = 10,
  group = "group",
  group.focal.compare = c("A", "B"),
  type.comparison = "between"
)

# Calculate reach within a single group
result_within <- clade_reach_coexistence(
  df.TS.TE = df_temporal,
  df.occ = df_occurrences,
  time.slice = 10,
  group = "group",
  group.focal.compare = c("A", "B"),
  type.comparison = "within"
)
} # }
```
