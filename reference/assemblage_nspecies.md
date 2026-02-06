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
