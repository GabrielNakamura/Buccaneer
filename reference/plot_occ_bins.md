# Plot Fossil Occurrence Records and Species Age Ranges

This function creates a visualization of fossil occurrence records
showing the temporal distribution of species across geological time.
Each species is represented by horizontal line segments spanning from
maximum (oldest) to minimum (youngest) ages for each occurrence record.
The plot can optionally color-code species by number of occurrences and
overlay geological time scale boundaries.

## Usage

``` r
plot_occ_bins(
  df.occ.fossil,
  species = "accepted_name",
  max.age = "max_ma",
  min.age = "min_ma",
  occ.gradient = TRUE,
  log = TRUE,
  name.scale.occ = "viridis",
  x.axis.name = "age range",
  y.axis.name = "species.name",
  show.species.name = FALSE,
  age.scheme = TRUE
)
```

## Arguments

- df.occ.fossil:

  A data frame containing fossil occurrence records with at least three
  columns: species names, maximum (oldest) age estimates, and minimum
  (youngest) age estimates. Each row represents a single occurrence
  record with age uncertainty.

- species:

  Character. The name of the column in `df.occ.fossil` containing
  species identifiers. Default is "accepted_name".

- max.age:

  Character. The name of the column in `df.occ.fossil` containing the
  maximum (oldest) age estimate for each occurrence record. Default is
  "max_ma".

- min.age:

  Character. The name of the column in `df.occ.fossil` containing the
  minimum (youngest) age estimate for each occurrence record. Default is
  "min_ma".

- occ.gradient:

  Logical. Should line segments be colored by the number of occurrences
  for each species? If TRUE (default), species with more occurrence
  records will have different colors based on the chosen color scale.

- log:

  Logical. Should the number of occurrences be log-transformed for color
  scaling? Default is TRUE. This can help visualize species with highly
  variable occurrence counts.

- name.scale.occ:

  Character. The name of the viridis color scale to use for coloring by
  occurrence count. Options include "viridis" (default), "magma",
  "plasma", "inferno", "cividis", "mako", "rocket", and "turbo".

- x.axis.name:

  Character. Label for the x-axis. Default is "age range".

- y.axis.name:

  Character. Label for the y-axis. Default is "species.name".

- show.species.name:

  Logical. Should species names be displayed on the y-axis? Default is
  FALSE. When FALSE, y-axis labels are hidden to reduce clutter when
  plotting many species.

- age.scheme:

  Logical. Should geological time period boundaries be overlaid as
  vertical dashed lines? Default is TRUE. Currently uses North American
  Land Mammal Ages (NALMAs) boundaries.

## Value

A ggplot2 object displaying the occurrence records. Species are arranged
vertically (y-axis) sorted by their maximum age, with time (age) on the
reversed x-axis (oldest to the left, youngest to the right). Each
horizontal segment represents one occurrence record spanning its age
uncertainty.

## Details

The function performs the following visualization steps:

1.  Orders species by their oldest occurrence (maximum age)

2.  Counts the number of occurrence records per species

3.  Optionally log-transforms occurrence counts for better visualization

4.  Creates horizontal segments for each occurrence spanning min to max
    age

5.  Colors segments by occurrence count using viridis color scales

6.  Optionally overlays geological time boundaries

Visualization features:

- **X-axis**: Time in millions of years (Ma), reversed so older ages
  appear on the left

- **Y-axis**: Species names (ordered by first appearance)

- **Color gradient**: Number of occurrences per species (with or without
  log transformation)

- **Vertical lines**: Geological time period boundaries (if enabled)

The default geological time scheme uses North American Land Mammal Ages
(NALMAs) boundaries. To use custom time boundaries, set
`age.scheme = FALSE` and modify the returned ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create example fossil occurrence data
df_fossils <- data.frame(
  accepted_name = rep(c("sp1", "sp2", "sp3"), each = 3),
  max_ma = c(100, 95, 90, 85, 80, 75, 70, 65, 60),
  min_ma = c(95, 90, 85, 80, 75, 70, 65, 60, 55)
)

# Basic plot with default settings
plot_occ_bins(df.occ.fossil = df_fossils)

# Plot with species names visible
plot_occ_bins(
  df.occ.fossil = df_fossils,
  show.species.name = TRUE
)

# Plot without log transformation of occurrence counts
plot_occ_bins(
  df.occ.fossil = df_fossils,
  log = FALSE
)

# Plot without geological time boundaries
plot_occ_bins(
  df.occ.fossil = df_fossils,
  age.scheme = FALSE
)

# Use different color scale
plot_occ_bins(
  df.occ.fossil = df_fossils,
  name.scale.occ = "plasma"
)

# Customize the plot further
p <- plot_occ_bins(df.occ.fossil = df_fossils)
p +
  ggplot2::labs(title = "Fossil Occurrences Through Time") +
  ggplot2::theme_minimal()

# Use custom column names
df_custom <- data.frame(
  taxon = rep(c("sp1", "sp2"), each = 2),
  oldest = c(100, 95, 90, 85),
  youngest = c(95, 90, 85, 80)
)

plot_occ_bins(
  df.occ.fossil = df_custom,
  species = "taxon",
  max.age = "oldest",
  min.age = "youngest"
)
} # }
```
