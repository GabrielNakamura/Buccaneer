# Check Occurrence Records Crossing Time Bin Boundaries

This function analyzes how fossil occurrence records align with
predefined time bins by classifying each record based on whether its age
uncertainty (minimum and maximum ages) falls completely within a bin,
crosses one boundary, or spans multiple bins. This is useful for
assessing data quality and understanding the temporal resolution of
fossil occurrences relative to time bins used in analyses.

## Usage

``` r
check_crossing_boundaries(
  df.occ.fossil,
  interval,
  species = "species",
  Max.age = "Maximum_Age",
  Min.age = "Minimum_Age",
  remove.sub.species = TRUE
)
```

## Arguments

- df.occ.fossil:

  A data frame containing fossil occurrence records with at least three
  columns: species names, maximum (oldest) age estimates, and minimum
  (youngest) age estimates. Each row represents a single occurrence
  record with age uncertainty.

- interval:

  Numeric vector. A sequence of time bin boundaries in descending order
  (from oldest to youngest). For example, `c(100, 90, 80, 70)` creates
  three time bins: 100-90, 90-80, and 80-70. Must have at least two
  elements.

- species:

  Character. The name of the column in `df.occ.fossil` containing
  species identifiers. Default is "species".

- Max.age:

  Character. The name of the column in `df.occ.fossil` containing the
  maximum (oldest) age estimate for each occurrence record. Default is
  "Maximum_Age".

- Min.age:

  Character. The name of the column in `df.occ.fossil` containing the
  minimum (youngest) age estimate for each occurrence record. Default is
  "Minimum_Age".

- remove.sub.species:

  Logical. Should subspecies-level identifications be removed, keeping
  only species-level records? Default is TRUE. This parameter is passed
  to
  [`clean_occ_fossil()`](https://gabrielnakamura.github.io/Buccaneer/reference/clean_occ_fossil.md).

## Value

A list containing two elements:

- summary_table:

  A named numeric vector showing the percentage of occurrence records in
  each category:

  - `within`: Records completely contained within a time bin

  - `upper_only`: Records crossing only the upper (older) boundary

  - `lower_only`: Records crossing only the lower (younger) boundary

  - `both`: Records spanning an entire time bin (crossing both
    boundaries)

- df_records_interval:

  A data frame containing all occurrence records with added columns:

  - `boundaries`: Classification of boundary crossing (one of the four
    categories above)

  - `time.interval`: The time bin identifier (e.g., "100-90")

  - All original columns from the cleaned occurrence data

## Details

The function classifies each occurrence record into one of four
categories:

1.  **within**: Both Max.age and Min.age fall within the time bin
    boundaries (Max.age \<= upper boundary AND Min.age \>= lower
    boundary)

2.  **upper_only**: The record crosses the upper (older) boundary
    (Max.age \> upper boundary BUT Min.age is within the bin)

3.  **lower_only**: The record crosses the lower (younger) boundary
    (Min.age \< lower boundary BUT Max.age is within the bin)

4.  **both**: The record spans the entire time bin, crossing both
    boundaries (Max.age \>= upper boundary AND Min.age \<= lower
    boundary)

Before classification, the function calls
[`clean_occ_fossil()`](https://gabrielnakamura.github.io/Buccaneer/reference/clean_occ_fossil.md)
to:

- Calculate midpoint ages

- Optionally remove subspecies-level identifications

- Standardize column names

Records crossing boundaries indicate age uncertainty that exceeds the
temporal resolution of the time bins. High percentages of crossing
records may suggest that finer time bins or different binning strategies
should be considered for the analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create example fossil occurrence data with age uncertainties
df_fossils <- data.frame(
  species = c("sp1", "sp2", "sp3", "sp4", "sp5"),
  Maximum_Age = c(95, 105, 88, 92, 110),
  Minimum_Age = c(85, 95, 82, 78, 70)
)

# Define time bins (100-90, 90-80, 80-70)
time_bins <- c(100, 90, 80, 70)

# Check boundary crossing patterns
results <- check_crossing_boundaries(
  df.occ.fossil = df_fossils,
  interval = time_bins
)

# View summary of boundary crossing
results$summary_table

# View detailed classification for all records
head(results$df_records_interval)

# Filter records that span entire bins
spanning_records <- subset(
  results$df_records_interval,
  boundaries == "both"
)

# Visualize boundary crossing patterns
library(ggplot2)
ggplot(results$df_records_interval,
       aes(x = time.interval, fill = boundaries)) +
  geom_bar(position = "fill") +
  labs(y = "Proportion", x = "Time Interval") +
  theme_minimal()

# Calculate statistics by time interval
library(dplyr)
interval_summary <- results$df_records_interval %>%
  group_by(time.interval, boundaries) %>%
  summarise(count = n(), .groups = "drop")
} # }
```
