# Compute Individual Species Coexistence based on reach criteria

Compute Individual Species Coexistence based on reach criteria

## Usage

``` r
IndivSpecies_reach_coexistence(
  df.TS.TE,
  df.occ,
  time.slice,
  round.digits,
  species,
  TS,
  TE,
  lat,
  lon,
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

- time.slice:

  Numeric. The time point representing each slice, typically the upper
  (older) boundary of the time bin.

- species:

  Character. The name of all species with at least one cooccurrence with
  another species.

- n.coexistence:

  Numeric. The number of co-occurring species in each time slice per
  species.
