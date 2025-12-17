# Calculate Species Co-occurrence Matrices Based on Time Slice Overlap

This auxiliary function computes species-by-species co-occurrence
matrices for each time slice based on shared site occupancy. Two species
are considered to co-occur if they are found at the same site during the
same time slice. The resulting matrices quantify how many sites each
pair of species shares, which can be used to calculate trait distances
among co-occurring species.

## Usage

``` r
comp_slice_occ_cooccurr(
  spp_slice,
  df.occ,
  species = "species",
  Max.age = "Max.age",
  Min.age = "Min.age"
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

## Value

A list with length equal to the number of time slices in `spp_slice`.
Each element is a square matrix (species × species) where:

- Row and column names are species names

- Diagonal elements represent the number of sites where each species
  occurs

- Off-diagonal elements represent the number of shared sites between
  pairs of species

- Values of 0 indicate no site overlap between species pairs
