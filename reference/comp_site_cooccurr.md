# Calculate Species Co-occurrence Matrices Based on Site Overlap Across Time Slices

This auxiliary function computes species-by-species co-occurrence
matrices for each time slice based on shared site occupancy. Two species
are considered to co-occur if they are found at the same site during the
same time slice. The resulting matrices quantify how many sites each
pair of species shares, which can be used to calculate trait distances
among co-occurring species.

## Usage

``` r
comp_site_cooccurr(
  spp_slice,
  df.occ,
  species = "species",
  Max.age = "Max.age",
  Min.age = "Min.age",
  site = "site"
)
```

## Arguments

- spp_slice:

  A named list where each element contains a character vector of species
  names present in a given time slice. Names of list elements should be
  time slice identifiers (typically numeric values). Usually generated
  by `calc_spp_slice()`.

- df.occ:

  A data frame containing fossil occurrence records with at least four
  columns: species names, minimum age, maximum age, and site location
  ID. Each row represents a single occurrence record at a specific site.
  Column names should match the parameters below.

- species:

  Character. The name of the column in `df.occ` containing species
  identifiers. Default is "species".

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

## Value

A list with length equal to the number of time slices in `spp_slice`.
Each element is a square matrix (species × species) where:

- Row and column names are species names

- Diagonal elements represent the number of sites where each species
  occurs

- Off-diagonal elements represent the number of shared sites between
  pairs of species

- Values of 0 indicate no site overlap between species pairs

## Details

The function performs the following steps for each time slice:

1.  Filters `df.occ` to include only species present in the time slice

2.  Further filters occurrences where the age range overlaps with the
    time slice (Max.age \>= time slice AND Min.age \<= time slice)

3.  Removes duplicate species records within sites

4.  Creates a site-by-species presence/absence matrix

5.  Computes the cross-product (t(matrix) \\ species-by-species
    co-occurrence matrix

The resulting co-occurrence matrix quantifies spatial overlap:

- Diagonal values = total number of sites occupied by each species

- Off-diagonal values = number of sites shared by species pairs

- This differs from simple presence/absence (0/1) co-occurrence

This function is typically used internally by
[`IndivSpec_site_distance()`](https://gabrielnakamura.github.io/Buccaneer/reference/IndivSpec_site_distance.md)
and related functions to determine which species co-occur at sites
before calculating trait-based community metrics.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create example species composition per time slice
spp_slice <- list(
  "100" = c("sp1", "sp2", "sp3"),
  "90" = c("sp1", "sp2", "sp4"),
  "80" = c("sp2", "sp3", "sp4")
)

# Create example occurrence data
df_occurrences <- data.frame(
  species = c("sp1", "sp1", "sp1", "sp2", "sp2", "sp3", "sp3", "sp4"),
  Max.age = c(100, 100, 95, 100, 90, 95, 95, 85),
  Min.age = c(90, 90, 85, 80, 75, 80, 85, 70),
  site = c("site1", "site2", "site3", "site1", "site2", "site1", "site3", "site2")
)

# Generate co-occurrence matrices for each time slice
cooccur_matrices <- comp_site_cooccurr(
  spp_slice = spp_slice,
  df.occ = df_occurrences
)

# View first time slice co-occurrence matrix
cooccur_matrices[[1]]

# Convert to binary co-occurrence (present/absent)
binary_cooccur <- ifelse(cooccur_matrices[[1]] >= 1, 1, 0)

# Calculate co-occurrence frequency for a species pair
shared_sites_sp1_sp2 <- cooccur_matrices[[1]]["sp1", "sp2"]

} # }
```
