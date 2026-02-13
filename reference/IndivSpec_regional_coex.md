# Calculate Regional Individual Species Coexistence Across Time Slices

This function computes the number of coexisting species at each time
slice in a regional scale analysis. It divides the temporal range into
discrete intervals and counts how many species overlap in each interval
based on their origination (TS) and extinction (TE) times.

## Usage

``` r
IndivSpec_regional_coex(
  df.TS.TE,
  time.slice,
  round.digits = 1,
  species = "species",
  TS = "TS",
  TE = "TE"
)
```

## Arguments

- df.TS.TE:

  A data frame containing species temporal data with at least three
  columns: species names, origination times, and extinction times.

- time.slice:

  Numeric. The time interval (in the same units as TS and TE) between
  consecutive time slices.

- round.digits:

  Integer. The number of decimal places to round time slice values.
  Default is 1. This affects the precision of temporal binning.

- species:

  Character. The name of the column in `df.TS.TE` containing species
  identifiers. Default is "species".

- TS:

  Character. The name of the column in `df.TS.TE` containing origination
  (start) times for each species. Default is "TS".

- TE:

  Character. The name of the column in `df.TS.TE` containing extinction
  (end) times for each species. Default is "TE".

## Value

A data frame with three columns:

- time.slice:

  Numeric. The time point representing each slice.

- species:

  Character. The name of each species.

- n.coexistence:

  Numeric. The number of species coexisting in each time slice for all
  species.

## Examples

``` r
# Example with fossil data
df <- data.frame(
  species = c("sp1", "sp2", "sp3"),
  TS = c(100, 90, 80),
  TE = c(50, 40, 30)
)
IndivSpec_regional_coex(df, time.slice = 10)
#>   time.slice species n.coexistence
#> 1        100     sp1             0
#> 2         90     sp1             1
#> 3         90     sp2             1
#> 4         40     sp2             1
#> 5         40     sp3             1
#> 6         30     sp3             0
```
