# Compute individual species coexistence based on reach criteria

Compute individual species coexistence based on reach criteria

## Usage

``` r
IndivSpec_reach_richness(
  df.TS.TE,
  df.occ,
  time.slice,
  round.digits = 1,
  species = "species",
  TS = "TS",
  TE = "TE",
  Max.age = "Max.age",
  Min.age = "Min.age",
  lat = "lat",
  lon = "lng"
)
```

## Arguments

- df.TS.TE:

  a data frame object containing at least three columns. Species names,
  origination time and extinction time for each species.

- df.occ:

  a data frame object containing the occurrence records for each
  species. This must have at least a column indicating the name of
  species, its minimum and maximum age estimate, and its site location
  ID.

- time.slice:

  scalar indicating the time interval between consecutive time slices.

- round.digits:

  scalar indicating the number of digits for time of origination and
  time for extinction.

- species:

  character indicating the name of the column of the data frame
  containing the species name information.

- TS:

  character indicating the name of the columns of the data frame
  containing the information on origination time.

- TE:

  character indicating the name of the column of the data frame
  containing the information on extinction time.

- Max.age:

  character indicating the name of the column containing the upper age
  limit for occurrence record.

- Min.age:

  character indicating the name of the column containing the lower age
  limit for occurrence record.

- lat:

  character indicating the name of the column containing the latitude of
  occurrence record in df.occ.

- lon:

  character indicating the name of the column containing the longitude
  of occurrence record in df.occ.

## Value

data frame containing the name of species, its mean coexistence value
calculated by each time slice considering reach criteria of
co-occurrence
