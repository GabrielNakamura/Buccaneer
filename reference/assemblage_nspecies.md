# Calculate the number of species in each grid/assemblage/site

Calculate the number of species in each grid/assemblage/site

## Usage

``` r
assemblage_nspecies(
  df.TS.TE,
  df.occ,
  time.slice,
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

  Data frame object containing the occurrence records for each species.
  This must have at least a column indicating the name of species, its
  minimum and maximum age estimate, and its site location ID.

- time.slice:

  Scalar indicating the time interval between consecutive time slices.

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
  on assemblage location.

## Value

A data frame with three columns:

- `time.slice`: time slice identifier

- `sites`: site or assemblage identifier

- `n.species`: number of species present at that site and time slice

## Examples

``` r
# Example species longevities
df_longevities <- data.frame(
  species = c("sp1", "sp2", "sp3", "sp4"),
  TS = c(100, 95, 90, 85),
  TE = c(60, 55, 50, 45)
)

# Example occurrence data
df_occurrences <- data.frame(
  species = c("sp1", "sp2", "sp3", "sp1", "sp4"),
  Max.age = c(90, 90, 90, 80, 80),
  Min.age = c(70, 70, 70, 60, 60),
  site = c("site1", "site1", "site2", "site2", "site1")
)

# Calculate number of species per site and time slice
assemblage_nspecies(
  df.TS.TE = df_longevities,
  df.occ = df_occurrences,
  time.slice = 10
)
#>    time.slice sites n.species
#> 1        <NA>  <NA>        NA
#> 2          90 site1         2
#> 3          90 site2         1
#> 4          80 site1         3
#> 5          80 site2         2
#> 6          70 site1         3
#> 7          70 site2         2
#> 8          60 site2         1
#> 9          60 site1         1
#> 10       <NA>  <NA>        NA
```
