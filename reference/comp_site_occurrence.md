# Calculate Site-Based Species Occurrence Across Time Slices

This auxiliary function determines which species occur at which sites
within each time slice, creating presence/absence matrices for fossil
assemblages. It filters occurrence records based on temporal overlap
with time slices and generates site-by-species matrices that can be used
for downstream community ecology analyses.

## Usage

``` r
comp_site_occurrence(
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
Each element is a data frame in wide format representing a
site-by-species presence/absence matrix, with the following structure:

- First column: `site` - Site identifiers

- Subsequent columns: Species names as column headers, with values of 1
  (present) or 0 (absent) indicating occurrence at each site

## Details

The function performs the following steps for each time slice:

1.  Filters `df.occ` to include only species present in the time slice
    (based on `spp_slice`)

2.  Further filters occurrences where the age range overlaps with the
    time slice (Max.age \>= time slice AND Min.age \<= time slice)

3.  Groups occurrences by site and removes duplicate species records

4.  Creates a presence/absence matrix with sites as rows and species as
    columns

5.  Returns the matrix in wide format suitable for community ecology
    analyses

This function is typically used internally to prepare occurrence data
for calculating co-occurrence patterns, trait distances, or other
community-level metrics at the site scale.

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
  species = c("sp1", "sp1", "sp2", "sp2", "sp3", "sp4"),
  Max.age = c(100, 95, 100, 90, 95, 85),
  Min.age = c(90, 85, 80, 75, 80, 70),
  site = c("site1", "site2", "site1", "site2", "site1", "site2")
)

# Generate site occurrence matrices for each time slice
site_matrices <- comp_site_occurrence(
  spp_slice = spp_slice,
  df.occ = df_occurrences
)

# View first time slice matrix
site_matrices[[1]]

# Convert to traditional matrix format for analysis
mat_slice1 <- as.matrix(site_matrices[[1]][, -1])
rownames(mat_slice1) <- site_matrices[[1]]$site

# Calculate species richness per site
richness_per_site <- rowSums(mat_slice1)

# Calculate beta diversity between sites
library(vegan)
beta_div <- vegdist(mat_slice1, method = "jaccard")
} # }
```
