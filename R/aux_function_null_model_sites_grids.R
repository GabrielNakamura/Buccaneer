#' Calculate Null Model for Trait Distances with Fixed Site Richness
#'
#' This function generates null distributions of trait distances by randomly
#' permuting species occurrences across sites while maintaining the observed
#' number of species per site (fixed row sums) and total occurrences per species
#' (fixed column sums). It computes mean trait distances for each permutation
#' to create null expectations for comparing against observed patterns.
#'
#' @param list_occurrence A list where each element represents a time slice and
#'     contains a data frame of species occurrences across sites or assemblages.
#'     Each data frame should have sites as rows and species as columns (plus
#'     one column for site identifiers). Presence is indicated by 1, absence by 0.
#' @param dist_matrix_trait A distance matrix (class \code{matrix} or \code{dist})
#'     containing pairwise trait distances between species. Row and column names
#'     must match species names in \code{list_occurrence}.
#' @param nearest.taxon Numeric or character. The number of nearest neighbors to
#'     consider when calculating mean distances. Use a numeric value (e.g., \code{1})
#'     for mean nearest neighbor distance (MNND), or \code{"all"} for mean pairwise
#'     distance (MPD).
#' @param nperm Integer. The number of permutations to generate for the null model.
#'     Default is 1000. Higher values provide more robust null distributions but
#'     increase computation time.
#'
#' @return A matrix with \code{nperm} rows and columns equal to the length of
#'     \code{list_occurrence}. Each row represents one permutation, and each column
#'     represents one time slice. Cell values contain mean trait distances calculated
#'     from the randomized occurrence data.
#'
#' @details
#' The function performs the following steps for each time slice:
#' \enumerate{
#'   \item Generates \code{nperm} random permutations of the occurrence matrix
#'         using \code{vegan::permatfull()} with fixed row and column margins
#'   \item For each permutation, creates a species co-occurrence matrix
#'   \item Filters the trait distance matrix to include only co-occurring species
#'   \item Calculates mean distances based on \code{nearest.taxon}:
#'         \itemize{
#'           \item If numeric: computes mean of the \code{nearest.taxon} nearest neighbors
#'           \item If "all": computes mean pairwise distance (MPD) among all co-occurring species
#'         }
#' }
#'
#' The null model maintains:
#' \itemize{
#'   \item Species richness at each site (row sums)
#'   \item Total number of occurrences for each species (column sums)
#'   \item Presence/absence structure (binary data)
#' }
#'
#' This approach tests whether observed trait distances deviate from expectations
#' under random community assembly constrained by species prevalence and site richness.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example occurrence data for two time slices
#' occ_slice1 <- data.frame(
#'   site = c("site1", "site2", "site3"),
#'   sp1 = c(1, 1, 0),
#'   sp2 = c(1, 0, 1),
#'   sp3 = c(0, 1, 1)
#' )
#'
#' occ_slice2 <- data.frame(
#'   site = c("site1", "site2", "site3"),
#'   sp1 = c(1, 0, 1),
#'   sp2 = c(1, 1, 0),
#'   sp4 = c(0, 1, 1)
#' )
#'
#' list_occ <- list(occ_slice1, occ_slice2)
#'
#' # Create trait distance matrix
#' traits <- c(1.2, 2.5, 3.1, 4.0)
#' names(traits) <- c("sp1", "sp2", "sp3", "sp4")
#' dist_trait <- as.matrix(dist(traits))
#'
#' # Calculate null model with 100 permutations for MNND
#' null_results_mnnd <- calc_null_model(
#'   list_occurrence = list_occ,
#'   dist_matrix_trait = dist_trait,
#'   nearest.taxon = 1,
#'   nperm = 100
#' )
#'
#' # Calculate null model for MPD
#' null_results_mpd <- calc_null_model(
#'   list_occurrence = list_occ,
#'   dist_matrix_trait = dist_trait,
#'   nearest.taxon = "all",
#'   nperm = 100
#' )
#'
#' # Calculate observed values (hypothetical)
#' observed_values <- c(2.3, 2.7)  # for two time slices
#'
#' # Calculate standardized effect sizes (SES)
#' ses <- (observed_values - colMeans(null_results_mpd)) /
#'        apply(null_results_mpd, 2, sd)
#'
#' # Test significance (two-tailed)
#' p_values <- colMeans(abs(null_results_mpd -
#'              rep(observed_values, each = nrow(null_results_mpd))) >=
#'              abs(observed_values - colMeans(null_results_mpd)))
#' }
calc_null_model <-
  function(list_occurrence, dist_matrix_trait, nearest.taxon, nperm = 1000){

    # data frame to receive the results
    df_null_res <- matrix(NA, nrow = nperm, ncol = length(list_occurrence),
                          dimnames = list(paste("perm", 1:nperm, sep = "_"),
                                          paste("slice", 1:length(list_occurrence), sep = "_")
                          )
    )

    for(i in 1:length(list_occurrence)){
      # i = 310
      # compute null matrix for each timeslice composition
      null_occ_site <-
        vegan::permatfull(m = list_occurrence[[i]][, -1],
                          fixedmar = "rows",
                          mtype = "prab",
                          times = nperm)

      # computing species co-occurrence matrix
      list_matrix_coocccur_null_site <-
        lapply(null_occ_site$perm, function(x){
          mat <- as.matrix(x)  # remove site column
          cooccur_matrix_sites <- t(mat) %*% mat  # species-by-species co-occurrence
          cooccur_matrix_sites
        })

      # matching species co-occurrence matrix and species distance matrix
      list_dist_matrix_null <-
        lapply(list_matrix_coocccur_null_site, function(x){
          species_row <- intersect(rownames(x), rownames(dist_matrix_trait))
          species_col <- intersect(colnames(x), colnames(dist_matrix_trait))
          dist_matrix <- dist_matrix_trait[species_row, species_col]
          dist_mat_filter <- ifelse(x == 1, dist_matrix, NA)
          dist_mat_filter
        })

      list_dist_matrix_null2 <-
        lapply(list_dist_matrix_null, function(x){
          diag(x) <- NA
          x
        })

      # in case of nearest taxon
      if(is.numeric(nearest.taxon)){
        # ordering the null distance matrix
        list_order_null <-
          lapply(list_dist_matrix_null2, function(x){
            apply(x, 1, simplify = TRUE, function(y) sort(y))
          })

        # this is for nearest taxon
        list_mean_nearest <-
          lapply(list_order_null, function(x){
            mean(unlist(lapply(x, function(x) x[1:nearest.taxon])), na.rm = T)
          })

        # unlisting and storing the null values
        df_null_res[, i] <- unlist(list_mean_nearest)
      } else{ # this is the case when the user set trait distance as all


        # removing lower triangle of the matrix to get rid of the duplicated values
        list_dist_lower_matrix_null <-
          lapply(list_dist_matrix_null, function(x){
            x[lower.tri(x)]
          })

        vec_null_mpd_slice <-
          unlist(lapply(list_dist_lower_matrix_null, function(x){
            mean(x, na.rm = TRUE)
          }))

        # mean_null_slice <- mean(vec_null_mpd_slice)
        df_null_res[, i] <- vec_null_mpd_slice
      }
    } # end loop for each timeslice
    return(df_null_res)
  }


