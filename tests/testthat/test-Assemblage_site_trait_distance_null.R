test_that("Assemblage_site_trait_distance_null calculates null trait distance", {
  # Load test data
  data("df_longevities_canidae")
  data("df_occurrence_canidae")

  # Prepare data
  df_ts_te <- df_longevities_canidae
  df_occ <- df_occurrence_canidae

  # Rename columns
  colnames(df_occ)[colnames(df_occ) == "MaxT"] <- "Max.age"
  colnames(df_occ)[colnames(df_occ) == "MinT"] <- "Min.age"
  colnames(df_occ)[colnames(df_occ) == "Site"] <- "site"

  # Create a simple distance matrix
  dist_mat <- dist(matrix(rnorm(nrow(df_ts_te)), ncol = 1))

  # Run the function with reduced permutations for speed
  result <- assemblage_site_trait_distance_null(
    df.TS.TE = df_ts_te,
    df.occ = df_occ,
    time.slice = 5,
    dist_matrix_trait = dist_mat,
    nperm = 10,
    nearest.taxon = "mpd"
  )

  # Test that result is a data frame or list
  expect_true(is.data.frame(result))

  # Test if the result is numeric
  expect_true(is.numeric(result$dist.obs) &
                is.numeric(result$dist.obs.z) &
                is.numeric(result$p.value) &
                is.character(result$time.slice)
              )

  # Test if the result is returning non null values
  expect_false(all(is.na(result$dist.obs)) &
                 all(is.na(result$dist.obs.z)) &
                 all(is.na(result$p.value))
               )

  # Test if all p values are returning as probabilities
  expect_all_true(result$p.value[which(!is.na(result$p.value))] <= 1)
})
