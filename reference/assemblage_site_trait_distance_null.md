# Compute Mean Pairwise Distances Between Species Cooccurring in Local Assemblages

Compute Mean Pairwise Distances Between Species Cooccurring in Local
Assemblages

## Usage

``` r
assemblage_site_trait_distance_null(
  df.TS.TE,
  df.occ,
  time.slice,
  dist_matrix_trait,
  round.digits = 10,
  species = "species",
  TS = "TS",
  TE = "TE",
  Max.age = "Max.age",
  Min.age = "Min.age",
  site = "site",
  nearest.taxon = c("mpd", "mnnd"),
  nperm = 1000,
  fixedmar = "rows",
  mtype = "prab"
)
```

## Arguments

- df.TS.TE:

  Data frame object containing at least four columns. Species names,
  origination time, extinction time and a trait value for each species.

- df.occ:

  A data frame containing fossil occurrence records with at least four
  columns: species names, minimum age, maximum age, and site location
  ID. Each row represents a single occurrence record at a specific site.

- time.slice:

  Numeric. The time interval (in the same units as TS and TE) between
  consecutive time slices for temporal binning.

- dist_matrix_trait:

  A distance matrix object (class `dist` or `matrix`) containing
  pairwise trait distances between species. Row and column names must
  match species names in `df.TS.TE`. If NULL, distances will be computed
  from the trait column using Euclidean distance.

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

- nearest.taxon:

  Character, indicating the mean distance to be used. If `"mpd"` mean
  pairwise distance will be calculated, if `mnnd` mean nearest neighbour
  distance will be computed, considering only the distant to the closest
  species in the morphospace.

- nperm:

  A scalar, indicating the number of permutations to be used in the null
  model

- fixedmar:

  Character, stating which row/column sums should be preserved in the
  assemblage composition matrix (sites in the rows and species in the
  columns). Options are "none", "rows", "columns" or "both".

- mtype:

  Character, indicanting the type of null model to be performed.
  `"prab"` is the option for presence-absence (detection/non-detection)
  matrix and `"count"` for matrix with count data.
