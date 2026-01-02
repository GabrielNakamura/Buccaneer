#' Calculate Mean Site-Level Species Coexistence Across Time Slices
#'
#' This function computes the mean number of co-occurring species at individual
#' sites across different time slices. For each time slice, it determines which
#' species co-occur at each site based on occurrence records and temporal ranges,
#' then calculates the mean and variance of coexistence across all species for
#' time slice.
#'
#' @param df.TS.TE A data frame containing species temporal data with at least
#'     three columns: species names, origination times (TS), and extinction
#'     times (TE). Additional columns may include group assignments.
#' @param df.occ A data frame containing fossil occurrence records with at least
#'     four columns: species names, minimum age, maximum age, and site location ID.
#'     Each row represents a single occurrence record at a specific site.
#' @param time.slice Numeric. The time interval (in the same units as TS and TE)
#'     between consecutive time slices for temporal binning.
#' @param round.digits Integer. The number of decimal places to round time slice
#'     values. Default is 1. This affects temporal binning precision.
#' @param species Character. The name of the column in \code{df.TS.TE} and
#'     \code{df.occ} containing species identifiers. Default is "species".
#' @param TS Character. The name of the column in \code{df.TS.TE} containing
#'     origination (first appearance) times for each species. Default is "TS".
#' @param TE Character. The name of the column in \code{df.TS.TE} containing
#'     extinction (last appearance) times for each species. Default is "TE".
#' @param Max.age Character. The name of the column in \code{df.occ} containing
#'     the maximum (oldest) age estimate for each occurrence record. Default is "Max.age".
#' @param Min.age Character. The name of the column in \code{df.occ} containing
#'     the minimum (youngest) age estimate for each occurrence record. Default is "Min.age".
#' @param site Character. The name of the column in \code{df.occ} containing
#'     site location identifiers. Default is "site".
#' @param remove.singletons Logical. Should singleton species (species occurring
#'     alone at a site with no co-occurring species) be excluded from mean and
#'     variance calculations? Default is TRUE. When TRUE, singletons are treated
#'     as NA; when FALSE, they contribute 0 to the mean.
#' @param group Character. The name of the column in \code{df.TS.TE} containing
#'     group assignments for species (e.g., clade, family). Required if using
#'     \code{group.focal.compare}. Default is NULL.
#' @param group.focal.compare Character vector of length 2. The first element
#'     specifies the focal group and the second specifies the comparison group.
#'     If NULL (default), coexistence is calculated across all species regardless
#'     of group membership.
#' @param type.comparison Character. Specifies the type of coexistence comparison:
#'     \itemize{
#'       \item \code{"between"}: Count only co-occurrences between species from
#'             the focal and comparison groups.
#'       \item \code{"within"}: Count only co-occurrences among species within
#'             the focal group.
#'       \item NULL (default): Count all co-occurrences regardless of group.
#'     }
#'
#' @return A data frame with three columns:
#'   \item{mean.coexistence}{Numeric. The mean number of co-occurring species in
#'       each time slice. This represents average local coexistence
#'       (excluding or including singletons based on \code{remove.singletons}).}
#'   \item{var.distance}{Numeric. The variance in the number of co-occurring
#'       species in each time slice.}
#'   \item{time.slice}{Numeric. The time point representing each slice, typically
#'       the upper (older) boundary of the time bin.}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Creates time slices from maximum TS to minimum TE
#'   \item Generates regional co-occurrence matrices using \code{aux_matrix_regional_coex()}
#'   \item Determines which species occur at which sites in each time slice
#'   \item Creates site-based co-occurrence matrices using \code{comp_site_cooccurr()}
#'   \item For each species at each site, counts the number of co-occurring species
#'   \item Calculates mean and variance of coexistence across species with
#'       cooccurrence in sites
#'   \item Optionally filters by group membership
#' }
#'
#' Coexistence calculation details:
#' \itemize{
#'   \item \strong{Self-coexistence}: Removed from counts (a species does not
#'         "co-occur" with itself)
#'   \item \strong{Singleton species}: Species with 0 co-occurring taxa
#'         \itemize{
#'           \item If \code{remove.singletons = TRUE}: Treated as NA (excluded from mean)
#'           \item If \code{remove.singletons = FALSE}: Contribute 0 to the mean
#'         }
#'   \item \strong{Group comparisons}: When using \code{group.focal.compare},
#'         only counts co-occurrences between/within specified groups
#' }
#'
#' This function calculates site-level (local) coexistence patterns, which differ
#' from regional-scale richness. For regional richness, see \code{clade_regional_coexistence()}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example fossil data
#' df_temporal <- data.frame(
#'   species = c("sp1", "sp2", "sp3", "sp4"),
#'   TS = c(100, 95, 90, 85),
#'   TE = c(50, 45, 40, 35),
#'   group = c("A", "A", "B", "B")
#' )
#'
#' df_occurrences <- data.frame(
#'   species = c("sp1", "sp1", "sp2", "sp3", "sp4", "sp4"),
#'   Max.age = c(100, 95, 95, 90, 85, 85),
#'   Min.age = c(90, 85, 85, 80, 75, 75),
#'   site = c("site1", "site2", "site1", "site1", "site2", "site3")
#' )
#'
#' # Calculate mean site coexistence through time
#' result <- clade_site_coexistence(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
#'   time.slice = 10
#' )
#'
#' # View results
#' head(result)
#'
#' # Plot mean coexistence through time
#' plot(result$time.slice,
#'      result$mean.coexistence,
#'      type = "l",
#'      xlab = "Time (Ma)",
#'      ylab = "Mean Site Richness",
#'      main = "Mean Species Coexistence at Sites Through Time")
#'
#' # Add variance as error bars
#' arrows(result$time.slice,
#'        result$mean.coexistence - sqrt(result$var.distance),
#'        result$time.slice,
#'        result$mean.coexistence + sqrt(result$var.distance),
#'        length = 0.05, angle = 90, code = 3)
#'
#' # Calculate coexistence including singletons
#' result_with_singletons <- clade_site_coexistence(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
#'   time.slice = 10,
#'   remove.singletons = FALSE
#' )
#'
#' # Calculate coexistence between groups
#' result_between <- clade_site_coexistence(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
#'   time.slice = 10,
#'   group = "group",
#'   group.focal.compare = c("A", "B"),
#'   type.comparison = "between"
#' )
#'
#' # Calculate coexistence within a single group
#' result_within <- clade_site_coexistence(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
#'   time.slice = 10,
#'   group = "group",
#'   group.focal.compare = c("A", "B"),
#'   type.comparison = "within"
#' )
#' }
clade_site_coexistence <-
  function(df.TS.TE,
           df.occ,
           time.slice,
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE",
           Max.age = "Max.age",
           Min.age = "Min.age",
           site = "site",
           remove.singletons = TRUE,
           group = NULL,
           group.focal.compare = NULL,
           type.comparison = NULL){

    # sub-setting the dataframe to keep or remove the group information

    if(!is.null(group) == TRUE){
      df.TS.TE <- df.TS.TE[, c(species, TS, TE, group)]
      colnames(df.TS.TE) <- c("species", "TS", "TE", "group")
    } else{
      df.TS.TE <- df.TS.TE[, c(species, TS, TE)]
      colnames(df.TS.TE) <- c("species", "TS", "TE")
    }

    # sub-setting and correcting the names of occurrence dataframe
    df_occ <-
      df.occ[, c(species, Max.age, Min.age, site)]
    vars <- list(species, Max.age, Min.age, site)
    name_vars <- c("species", "Max.age", "Min.age", "site")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df_occ) <- column.names
    df_occ$site <- as.factor(df_occ$site)

    # Generating time intervals used to compute temporal coexistence
    seq_interval <- seq(from = ceiling(max(df.TS.TE[, "TS"])),
                        to = ceiling(min(df.TS.TE[, "TE"])),
                        by = -time.slice)


    # Time coexistence matrix for all species
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE,
                               time.slice,
                               round.digits = round.digits,
                               species = "species",
                               TS = "TS",
                               TE = "TE")

    # modified coexistence matrix containing group comparison
    if(!is.null(group.focal.compare) == TRUE){
      focal <- group.focal.compare[1]
      compare <- group.focal.compare[2]
      spp_focal <- df.TS.TE[which(df.TS.TE$group == focal), "species"]$species
      spp_compare <- df.TS.TE[which(df.TS.TE$group == compare), "species"]$species

      if(type.comparison == "between"){# comparison between species of two groups
        matrix_coex <- lapply(matrix_coex, function(x) x[spp_focal, spp_compare]) # focal species in lines and comparison in columns
      }
      if(type.comparison == "within"){ # comparison only within the focal group
        matrix_coex <- lapply(matrix_coex, function(x) x[spp_focal, spp_focal])
      }
    } else{
      matrix_coex <- matrix_coex
    }

    # species composition at each timeslice
    spp_slice <-
      lapply(matrix_coex, function(x){
        names(which(rowSums(x) >= 1))
      })

    # naming list with time slices
    names(spp_slice) <- format(seq_interval, trim = TRUE, scientific = FALSE)

    # calculating matrix of species cooccurrence for site
    list_matrix_cooccur_site <-
      comp_site_cooccurr(spp_slice = spp_slice, df.occ = df_occ)
    # naming list with time slices
    names(list_matrix_cooccur_site) <- format(seq_interval, trim = TRUE, scientific = FALSE)

    # only presence-absence
    list_matrix_cooccur_site2 <- lapply(list_matrix_cooccur_site, function(x) ifelse(x >= 1, 1, 0))
    names(list_matrix_cooccur_site2) <- format(seq_interval, trim = TRUE, scientific = FALSE)

    # number of coexistences including self coexistence
    list_n_coex_all <- lapply(list_matrix_cooccur_site2, function(x) rowSums(x)) # including self coex
    list_n_coex_all2 <- lapply(list_n_coex_all, function(x) x - 1) # removing self coexistence

    if(remove.singletons != TRUE){
      mean_species_site_coex <- lapply(list_n_coex_all2, function(x) mean(x)) # including zeros (singletons)
      var_species_site_coex <- lapply(list_n_coex_all2, function(x) var(x)) # including zeros (singletons)
    } else{
      list_n_coex_all3 <- lapply(list_n_coex_all2, function(x) ifelse(x == 0, NA, x))
      mean_species_site_coex <- lapply(list_n_coex_all3, function(x) mean(x, na.rm = TRUE)) # removing singletons
      var_species_site_coex <- lapply(list_n_coex_all2, function(x) var(x, na.rm = TRUE))
    }

    # function output - dataframe
    df_res <-
      data.frame(mean.coexistence = unlist(mean_species_site_coex),
                 var.distance = unlist(var_species_site_coex),
                 time.slice = seq_interval)
    return(df_res)
  }
