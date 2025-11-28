# Clean and Flag Fossil Occurrence Records Based on Multiple Criteria

This function processes and standardizes fossil occurrence records by
calculating age estimates (midpoint, upper, or lower bounds), computing
age uncertainties, flagging records with large age ranges, optionally
computing species temporal ranges (origination and extinction times),
and identifying subspecies-level identifications. It provides a
comprehensive data cleaning pipeline for paleobiological analyses.

## Usage

``` r
clean_occ_fossil(
  df.occ.fossil,
  method.ages = c("midpoint", "upper", "lower"),
  thresh.age.range = 10,
  species = "species",
  Max.age = "Maximum_Age",
  Min.age = "Minimum_Age",
  remove.sub.species = TRUE,
  comp.TS.TE = TRUE,
  lat = NULL,
  lng = NULL,
  site = NULL,
  group = NULL,
  trait = NULL
)
```

## Arguments

- df.occ.fossil:

  A data frame containing fossil occurrence records with at least three
  columns: species names, maximum (oldest) age estimates, and minimum
  (youngest) age estimates. Additional columns may include spatial
  coordinates, site information, group assignments, and trait values.

- method.ages:

  Character. The method used to estimate representative ages from
  occurrence records. Options include:

  - `"midpoint"` (default): Calculate the midpoint between max and min
    ages

  - `"upper"`: Use the maximum (oldest) age

  - `"lower"`: Use the minimum (youngest) age

  Note: Currently, the function always calculates midpoint regardless of
  this parameter.

- thresh.age.range:

  Numeric. The threshold for flagging occurrence records with large age
  uncertainties (in millions of years or appropriate time units).
  Records with age ranges (Max.age - Min.age) greater than or equal to
  this value are flagged as `flag.age.range = "TRUE"`. Default is 10.

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

  Logical. Should subspecies-level identifications be flagged? Default
  is TRUE. When TRUE, adds a `subspecies` column identifying records
  with three or more words in the species name (e.g., "Genus species
  subspecies").

- comp.TS.TE:

  Logical. Should species-level temporal ranges (origination and
  extinction times) be computed? Default is TRUE. When TRUE, calculates
  TS (maximum of Max.age) and TE (minimum of Min.age) for each species.

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

## Value

A data frame containing all original columns plus the following added
columns:

- midpoint:

  Numeric. The midpoint age calculated as (Max.age + Min.age) / 2.

- age.range:

  Numeric. The age uncertainty calculated as Max.age - Min.age.

- flag.age.range:

  Character. "TRUE" if age.range \>= thresh.age.range, "FALSE"
  otherwise. Flags occurrences with potentially problematic age
  uncertainties.

- TS:

  Numeric (optional, if comp.TS.TE = TRUE). Species origination time,
  calculated as the maximum Max.age across all occurrences of that
  species.

- TE:

  Numeric (optional, if comp.TS.TE = TRUE). Species extinction time,
  calculated as the minimum Min.age across all occurrences of that
  species.

- subspecies:

  Character (optional, if remove.sub.species = TRUE). "subspecies" if
  the species name has three or more words, "species" otherwise.

## Details

The function performs the following data cleaning and flagging steps:

1.  Subsets and standardizes column names

2.  Converts Max.age and Min.age to numeric (removing non-numeric
    values)

3.  Calculates midpoint ages: (Max.age + Min.age) / 2

4.  Calculates age uncertainties: Max.age - Min.age

5.  Flags records with age.range \>= thresh.age.range

6.  Optionally computes species temporal ranges (TS and TE)

7.  Optionally identifies subspecies-level records

8.  Removes records with NA values in TS or TE (with warning)

Age uncertainty considerations:

- **Purpose**: Large age ranges indicate temporal imprecision

- **Threshold**: User-defined via `thresh.age.range`

- **Action**: Records are flagged but not removed, allowing users to
  decide how to handle them in downstream analyses

- **Typical thresholds**: 5-10 Ma for Cenozoic studies, 10-20 Ma for
  Mesozoic studies

Subspecies identification:

- Uses word count in species names (binomial = 2 words, trinomial \>= 3
  words)

- May require manual verification as some valid species names have
  multiple words (e.g., "Genus species complex")

- Useful for standardizing taxonomic resolution in analyses

Missing value handling:

- Non-numeric age values are converted to NA

- Records with NA in TS or TE are removed (when comp.TS.TE = TRUE)

- A warning is issued when NAs are found and removed

## Examples

``` r
if (FALSE) { # \dontrun{
# Create example fossil occurrence data
df_fossils <- data.frame(
  species = c("Genus species", "Genus species subspecies", "Genus other"),
  Maximum_Age = c(100, 95, 90),
  Minimum_Age = c(85, 80, 88),
  lat = c(10, 15, 20),
  lng = c(-50, -55, -60),
  site = c("A", "B", "C"),
  group = c("G1", "G1", "G2"),
  trait = c(1.2, 2.5, 3.1)
)

# Basic cleaning with default settings
cleaned_data <- clean_occ_fossil(df.occ.fossil = df_fossils)
head(cleaned_data)

# View flagged records with large age uncertainties
flagged <- subset(cleaned_data, flag.age.range == "TRUE")
flagged

# Identify subspecies-level records
subspecies_records <- subset(cleaned_data, subspecies == "subspecies")
subspecies_records

# Use a more stringent age uncertainty threshold
cleaned_strict <- clean_occ_fossil(
  df.occ.fossil = df_fossils,
  thresh.age.range = 5
)

# Clean without computing TS and TE
cleaned_no_ranges <- clean_occ_fossil(
  df.occ.fossil = df_fossils,
  comp.TS.TE = FALSE
)

# Clean without flagging subspecies
cleaned_all_taxa <- clean_occ_fossil(
  df.occ.fossil = df_fossils,
  remove.sub.species = FALSE
)

# Use custom column names
df_custom <- data.frame(
  taxon = c("sp1", "sp2", "sp3"),
  oldest = c(100, 95, 90),
  youngest = c(85, 80, 75)
)

cleaned_custom <- clean_occ_fossil(
  df.occ.fossil = df_custom,
  species = "taxon",
  Max.age = "oldest",
  Min.age = "youngest"
)

# Filter out flagged records for downstream analysis
library(dplyr)
reliable_records <- cleaned_data %>%
  filter(flag.age.range == "FALSE")

# Summary of data quality
quality_summary <- cleaned_data %>%
  group_by(flag.age.range) %>%
  summarise(
    n_records = n(),
    mean_age_range = mean(age.range),
    n_species = n_distinct(species)
  )
quality_summary
} # }
```
