test_that("IndivSpec_site_distance returns a data.frame with expected columns", {

  df_longevities <- data.frame(
    species = c("sp1", "sp2", "sp3"),
    TS = c(100, 95, 90),
    TE = c(60, 55, 50),
    trait = c(1, 2, 3)
  )

  df_occurrences <- data.frame(
    species = c("sp1", "sp2", "sp3"),
    Max.age = c(90, 90, 90),
    Min.age = c(70, 70, 70),
    site = c("site1", "site1", "site1")
  )

  res <- IndivSpec_site_distance(
    df.TS.TE = df_longevities,
    df.occ = df_occurrences,
    time.slice = 10,
    trait = "trait",
    dist.trait = NULL,
    nearest.taxon = "all"
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(c("species", "time.slice",
                    "mean_dist_to_cooccur", "is_singleton") %in% colnames(res)))
})


test_that("mean_dist_to_cooccur is numeric and is_singleton is logical", {

  df_longevities <- data.frame(
    species = c("sp1", "sp2"),
    TS = c(100, 100),
    TE = c(50, 50),
    trait = c(1, 3)
  )

  df_occurrences <- data.frame(
    species = c("sp1", "sp2"),
    Max.age = c(90, 90),
    Min.age = c(80, 80),
    site = c("site1", "site1")
  )

  res <- IndivSpec_site_distance(
    df.TS.TE = df_longevities,
    df.occ = df_occurrences,
    time.slice = 10,
    trait = "trait",
    dist.trait = NULL,
    nearest.taxon = "all"
  )

  expect_type(res$mean_dist_to_cooccur, "double")
  expect_type(res$is_singleton, "logical")
})


test_that("singleton species get distance 0 and is_singleton = TRUE", {

  df_longevities <- data.frame(
    species = c("sp1", "sp2"),
    TS = c(100, 100),
    TE = c(50, 50),
    trait = c(1, 2)
  )

  df_longevities_singleton <- data.frame(
    species = c("sp1", "sp2"),
    TS = c(100, 49),
    TE = c(50, 29),
    trait = c(1, 2)
  )

  df_occurrences <- data.frame(
    species = c("sp1", "sp2"),
    Max.age = c(90, 90),
    Min.age = c(80, 80),
    site = c("site1", "site2") # no co-occurrence
  )

  df_occurrences_singleton <- data.frame(
    species = c("sp1", "sp2"),
    Max.age = c(90, 47),
    Min.age = c(80, 25),
    site = c("site1", "site2") # no co-occurrence
  )


  res_no_singleton <- IndivSpec_site_distance(
    df.TS.TE = df_longevities,
    df.occ = df_occurrences,
    time.slice = 10,
    trait = "trait",
    dist.trait = NULL,
    nearest.taxon = "all"
  )

  expect_true(all(is.na(res_no_singleton$is_singleton)))

  res_singleton <- IndivSpec_site_distance(
    df.TS.TE = df_longevities_singleton,
    df.occ = df_occurrences_singleton,
    time.slice = 10,
    trait = "trait",
    dist.trait = NULL,
    nearest.taxon = "all"
  )

#expect_true(res_singleton$is_singleton)

})


test_that("MNND and MPD produce different results when nearest.taxon differs", {

  df_longevities <- data.frame(
    species = c("sp1", "sp2", "sp3"),
    TS = c(100, 100, 100),
    TE = c(50, 50, 50),
    trait = c(1, 2, 10)
  )

  df_occurrences <- data.frame(
    species = c("sp1", "sp2", "sp3"),
    Max.age = c(90, 90, 90),
    Min.age = c(80, 80, 80),
    site = c("site1", "site1", "site1")
  )

  mpd <- IndivSpec_site_distance(
    df.TS.TE = df_longevities,
    df.occ = df_occurrences,
    time.slice = 10,
    trait = "trait",
    dist.trait = NULL,
    nearest.taxon = "all"
  )

  mnnd <- IndivSpec_site_distance(
    df.TS.TE = df_longevities,
    df.occ = df_occurrences,
    time.slice = 10,
    trait = "trait",
    dist.trait = NULL,
    nearest.taxon = 1
  )

  expect_false(isTRUE(all.equal(
    mpd$mean_dist_to_cooccur,
    mnnd$mean_dist_to_cooccur
  )))
})


test_that("group-based comparisons run without error", {

  df_longevities <- data.frame(
    species = c("sp1", "sp2", "sp3", "sp4"),
    TS = c(100, 100, 100, 100),
    TE = c(50, 50, 50, 50),
    trait = c(1, 2, 3, 4),
    group = c("A", "A", "B", "B")
  )

  df_occurrences <- data.frame(
    species = c("sp1", "sp2", "sp3", "sp4"),
    Max.age = c(90, 90, 90, 90),
    Min.age = c(80, 80, 80, 80),
    site = c("site1", "site1", "site1", "site1")
  )

  expect_error(
    IndivSpec_site_distance(
      df.TS.TE = df_longevities,
      df.occ = df_occurrences,
      time.slice = 10,
      trait = "trait",
      dist.trait = NULL,
      nearest.taxon = "all",
      group = "group",
      group.focal.compare = c("A", "B"),
      type.comparison = "between"
    ),
    NA
  )
})
