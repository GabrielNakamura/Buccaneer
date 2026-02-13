# Calculate Null Model for Trait Distances with Fixed Site Richness

This function generates null distributions of trait distances by
randomly permuting species occurrences across sites while maintaining
the observed number of species per site (fixed row sums) and total
occurrences per species (fixed column sums). It computes mean trait
distances for each permutation to create null expectations for comparing
against observed patterns.

## Usage

``` r
calc_null_model(
  list_occurrence,
  dist_matrix_trait,
  nearest.taxon,
  nperm = 1000
)
```

## Arguments

- list_occurrence:

  A list where each element represents a time slice and contains a data
  frame of species occurrences across sites or assemblages. Each data
  frame should have sites as rows and species as columns (plus one
  column for site identifiers). Presence is indicated by 1, absence by
  0.

- dist_matrix_trait:

  A distance matrix (class `matrix` or `dist`) containing pairwise trait
  distances between species. Row and column names must match species
  names in `list_occurrence`.

- nearest.taxon:

  Numeric or character. The number of nearest neighbors to consider when
  calculating mean distances. Use a numeric value (e.g., `1`) for mean
  nearest neighbor distance (MNND), or `"all"` for mean pairwise
  distance (MPD).

- nperm:

  Integer. The number of permutations to generate for the null model.
  Default is 1000. Higher values provide more robust null distributions
  but increase computation time.

## Value

A matrix with `nperm` rows and columns equal to the length of
`list_occurrence`. Each row represents one permutation, and each column
represents one time slice. Cell values contain mean trait distances
calculated from the randomized occurrence data.

## Details

The function performs the following steps for each time slice:

1.  Generates `nperm` random permutations of the occurrence matrix using
    [`vegan::permatfull()`](https://vegandevs.github.io/vegan/reference/permatfull.html)
    with fixed row and column margins

2.  For each permutation, creates a species co-occurrence matrix

3.  Filters the trait distance matrix to include only co-occurring
    species

4.  Calculates mean distances based on `nearest.taxon`:

    - If numeric: computes mean of the `nearest.taxon` nearest neighbors

    - If "all": computes mean pairwise distance (MPD) among all
      co-occurring species

The null model maintains:

- Species richness at each site (row sums)

- Total number of occurrences for each species (column sums)

- Presence/absence structure (binary data)

This approach tests whether observed trait distances deviate from
expectations under random community assembly constrained by species
prevalence and site richness.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create example occurrence data for two time slices
occ_slice1 <- data.frame(
  site = c("site1", "site2", "site3"),
  sp1 = c(1, 1, 0),
  sp2 = c(1, 0, 1),
  sp3 = c(0, 1, 1)
)

occ_slice2 <- data.frame(
  site = c("site1", "site2", "site3"),
  sp1 = c(1, 0, 1),
  sp2 = c(1, 1, 0),
  sp4 = c(0, 1, 1)
)

list_occ <- list(occ_slice1, occ_slice2)

# Create trait distance matrix
traits <- c(1.2, 2.5, 3.1, 4.0)
names(traits) <- c("sp1", "sp2", "sp3", "sp4")
dist_trait <- as.matrix(dist(traits))

# Calculate null model with 100 permutations for MNND
null_results_mnnd <- calc_null_model(
  list_occurrence = list_occ,
  dist_matrix_trait = dist_trait,
  nearest.taxon = 1,
  nperm = 100
)

# Calculate null model for MPD
null_results_mpd <- calc_null_model(
  list_occurrence = list_occ,
  dist_matrix_trait = dist_trait,
  nearest.taxon = "all",
  nperm = 100
)

# Calculate observed values (hypothetical)
observed_values <- c(2.3, 2.7)  # for two time slices

# Calculate standardized effect sizes (SES)
ses <- (observed_values - colMeans(null_results_mpd)) /
       apply(null_results_mpd, 2, sd)

# Test significance (two-tailed)
p_values <- colMeans(abs(null_results_mpd -
             rep(observed_values, each = nrow(null_results_mpd))) >=
             abs(observed_values - colMeans(null_results_mpd)))
} # }
```
