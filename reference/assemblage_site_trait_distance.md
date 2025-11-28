# Compute mean distance for species co occurring in sites

Compute mean distance for species co occurring in sites

## Usage

``` r
assemblage_site_trait_distance(
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

  Data frame object containing at least four columns. Species names,
  origination time, extinction time and a trait value for each species.

- df.occ:

  a data frame object containing the occurrence records for each
  species. This must have at least a column indicating the name of
  species, its minimum and maximum age estimate, and its site location
  ID.

- time.slice:

  Scalar indicating the time interval between consecutive time slices.

- dist.trait:

  A dist object containing the clade pairwise distance matrices. The
  name of the clades in this object must be equal to the name of the
  clades in df.TS.TE data frame.

- nearest.taxon:

  A scalar indicating the number of nearest species/genus that will be
  used. 1 computes mnnd metric and the option "all" computes mpd.

- group:

  Character indicating the name of the column that contain the groups
  that will be used in comparison.

- group.focal.compare:

  Character vector indicating the focal (first element) and comparison
  (second element) groups used in the calculation. If NULL, the default,
  the metrics will be calculated using all clades.

- type.comparison:

  Character. It can be "between" to compute distances only between
  species/genus of two groups or "within" to calculate distance only
  inside the focal group. If null the distance is computed considering
  all clades together

- trait:

  Character indicating the name of the column containing values of the
  traits for each species. If NULL, the default, the user must provide a
  distance matrix.

- round.digits:

  Scalar indicating the precision of time slices.

- species:

  Character indicating the name of the column of the data frame
  containing the species name information.

- TS:

  Character indicating the name of the columns of the data frame
  containing the information on origination time.

- TE:

  Character indicating the name of the column of the data frame
  containing the information on extinction time.

- Max.age:

  Character indicating the name of the column containing the upper age
  limit for occurrence record.

- Min.age:

  Character indicating the name of the column containing the lower age
  limit for occurrence record.

- site:

  Character indicating the name of the column containing the information
  on site location.
