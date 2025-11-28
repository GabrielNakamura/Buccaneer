# Get Summary Information for Species Across Time Slices

This function processes fossil occurrence data to generate summaries
across discrete time slices. It calculates species temporal ranges
(origination and extinction times) from occurrence records, identifies
which species are present in each time slice, counts occurrence records
per species, and creates species-by-species co-occurrence matrices. The
function also flags occurrences with large age uncertainties based on a
user-defined threshold.

## Usage

``` r
get_summary_slices(
  df.TS.TE,
  timeframe,
  method.ages = c("midpoint", "upper", "lower"),
  thresh.age.range = 10,
  species = "species",
  Max.age = "Maximum_Age",
  Min.age = "Minimum_Age",
  TS = NULL,
  TE = NULL,
  lat = NULL,
  lng = NULL,
  site = NULL,
  group = NULL,
  trait = NULL,
  ...
)
```

## Arguments

- df.TS.TE:

  A data frame containing fossil occurrence records with at least three
  columns: species names, maximum (oldest) age estimates, and minimum
  (youngest) age estimates. Additional columns may include spatial
  coordinates, site information, group assignments, and trait values.

- timeframe:

  Numeric. The time interval (in millions of years or appropriate time
  units) between consecutive time slices. Negative values will create
  backwards intervals from the oldest to youngest occurrences.

- method.ages:

  Character. The method used to estimate species ages from occurrence
  records. Options include:

  - `"midpoint"` (default): Use the midpoint between max and min ages

  - `"upper"`: Use the maximum (oldest) age

  - `"lower"`: Use the minimum (youngest) age

- thresh.age.range:

  Numeric. The threshold for flagging occurrence records with large age
  uncertainties. Records with age ranges (Max.age - Min.age) greater
  than or equal to this value are flagged. Default is 10.

- species:

  Character. The name of the column in `df.TS.TE` containing species
  identifiers. Default is "species".

- Max.age:

  Character. The name of the column in `df.TS.TE` containing the maximum
  (oldest) age estimate for each occurrence record. Default is
  "Maximum_Age".

- Min.age:

  Character. The name of the column in `df.TS.TE` containing the minimum
  (youngest) age estimate for each occurrence record. Default is
  "Minimum_Age".

- TS:

  Character. The name of the column containing origination (first
  appearance) times. Default is NULL. If NULL, TS is calculated as the
  maximum midpoint age for each species.

- TE:

  Character. The name of the column containing extinction (last
  appearance) times. Default is NULL. If NULL, TE is calculated as the
  minimum midpoint age for each species.

- lat:

  Character. The name of the column containing latitude coordinates.
  Default is NULL. If provided, latitude information is retained in
  output.

- lng:

  Character. The name of the column containing longitude coordinates.
  Default is NULL. If provided, longitude information is retained in
  output.

- site:

  Character. The name of the column containing site location
  identifiers. Default is NULL. If provided, site information is
  retained in output.

- group:

  Character. The name of the column containing group assignments for
  species (e.g., clade, family). Default is NULL. If provided, group
  information is retained in output.

- trait:

  Character. The name of the column containing trait values for species.
  Default is NULL. If provided, trait information is retained in output.

- ...:

  Additional arguments (currently not used but reserved for future
  extensions).

## Value

A list containing two elements:

- df_sub_bin:

  A named list of data frames, one for each time slice. Each data frame
  contains all occurrence records present in that slice, with added
  columns:

  - `midpoint`: Numeric. The midpoint age of each occurrence

  - `age.range`: Numeric. The age uncertainty (Max.age - Min.age)

  - `flag.age.range`: Character. "TRUE" if age.range \>=
    thresh.age.range, "FALSE" otherwise

  - `TS`: Numeric. Species origination time (maximum midpoint)

  - `TE`: Numeric. Species extinction time (minimum midpoint)

  - `count_records`: Integer. Number of occurrence records for each
    species in that time slice

  List names follow the format "slice_X" where X is the time slice
  value.

- coex_matrix:

  A named list of binary co-occurrence matrices, one for each time
  slice. Each matrix:

  - Has dimensions: n × n, where n = total number of unique species in
    the dataset

  - Contains 1 if a species is present in the time slice, 0 if absent

  - Has row and column names as species identifiers

  List names follow the format "slice_X".

## Details

The function performs the following steps:

1.  Validates input data structure (must be a data frame)

2.  Calculates midpoint ages for each occurrence: (Max.age + Min.age) /
    2

3.  Calculates age uncertainties: Max.age - Min.age

4.  Flags occurrences with age uncertainties \>= thresh.age.range

5.  Computes species temporal ranges (TS and TE) from occurrence
    midpoints

6.  Creates time slices from oldest to youngest occurrence

7.  For each time slice, identifies species present (TS \>= slice AND TE
    \<= slice)

8.  Counts occurrence records per species per time slice

9.  Generates binary co-occurrence matrices for each time slice

Age uncertainty handling:

- **flag.age.range**: Identifies potentially problematic occurrences
  with large temporal uncertainties

- **Threshold**: Adjustable via `thresh.age.range` parameter

- **Use case**: Helps assess data quality and decide whether to filter
  or weight occurrences differently

Missing value handling:

- NA values in TS or TE columns trigger a warning and are removed

- Non-numeric TS or TE values trigger an error

This function is useful for preparing fossil occurrence data for
downstream analyses such as diversity dynamics, trait evolution, or
biogeographic patterns through time.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create example fossil occurrence data
df_fossils <- data.frame(
  species = c("sp1", "sp2", "sp3", "sp1", "sp2", "sp4"),
  Maximum_Age = c(100, 95, 90, 88, 85, 80),
  Minimum_Age = c(95, 90, 85, 83, 80, 75),
  lat = c(10, 15, 20, 12, 18, 25),
  lng = c(-50, -55, -60, -52, -57, -65),
  site = c("A", "B", "C", "A", "B", "D"),
  group = c("G1", "G1", "G2", "G1", "G1", "G2"),
  trait = c(1.2, 2.5, 3.1, 1.2, 2.5, 4.0)
)

# Get summary for 5 Ma time slices
results <- get_summary_slices(
  df.TS.TE = df_fossils,
  timeframe = 5,
  thresh.age.range = 10
)

# View occurrences in the first time slice
head(results$df_sub_bin[[1]])

# View co-occurrence matrix for the first time slice
results$coex_matrix[[1]]

# Count species richness per time slice
richness_per_slice <- sapply(results$df_sub_bin, function(df) {
  length(unique(df$species))
})
richness_per_slice

# Count total occurrences per time slice
records_per_slice <- sapply(results$df_sub_bin, nrow)
records_per_slice

# Identify occurrences with large age uncertainties
flagged_occurrences <- do.call(rbind, results$df_sub_bin) %>%
  filter(flag.age.range == "TRUE")
head(flagged_occurrences)

# Calculate sampling intensity through time
sampling_intensity <- data.frame(
  slice = names(results$df_sub_bin),
  n_records = records_per_slice,
  n_species = richness_per_slice,
  records_per_species = records_per_slice / richness_per_slice
)

# Plot richness through time
library(ggplot2)
ggplot(sampling_intensity, aes(x = slice, y = n_species)) +
  geom_col() +
  labs(x = "Time Slice", y = "Species Richness") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Use custom age method and threshold
results_upper <- get_summary_slices(
  df.TS.TE = df_fossils,
  timeframe = 5,
  method.ages = "upper",
  thresh.age.range = 5
)
} # }
```
