#' Calculate Species-Level Co-occurrence at Sites Across Time Slices
#'
#' This function computes the number of co-occurring species for each individual
#' species at fossil sites across different time slices. For each species present
#' at sites during a time slice, it counts how many other species it co-occurs with,
#' providing species-specific co-occurrence patterns through time. The function
#' can perform comparisons between taxonomic groups or within a single group.
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
#'     alone at sites with no co-occurring species) be excluded from the output?
#'     Default is TRUE. When TRUE, species with \code{n.coexistence = 0} are removed.
#' @param group Character. The name of the column in \code{df.TS.TE} containing
#'     group assignments for species (e.g., clade, family). Required if using
#'     \code{group.focal.compare}. Default is NULL.
#' @param group.focal.compare Character vector of length 2. The first element
#'     specifies the focal group and the second specifies the comparison group.
#'     If NULL (default), co-occurrence is calculated across all species regardless
#'     of group membership.
#' @param type.comparison Character. Specifies the type of co-occurrence comparison:
#'     \itemize{
#'       \item \code{"between"}: Count only co-occurrences between species from
#'             the focal and comparison groups.
#'       \item \code{"within"}: Count only co-occurrences among species within
#'             the focal group.
#'       \item NULL (default): Count all co-occurrences regardless of group.
#'     }
#'
#' @return A data frame with three columns:
#'   \item{time.slice}{Character. The time slice identifier (e.g., "100" for
#'       the time slice at 100 Ma).}
#'   \item{species}{Character. The name of each species.}
#'   \item{n.coexistence}{Integer. The number of other species that co-occur
#'       with the focal species at sites during that time slice. Self-coexistence
#'       is excluded (a species does not count as co-occurring with itself).}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Creates time slices from maximum TS to minimum TE
#'   \item Generates regional co-occurrence matrices using \code{aux_matrix_regional_coex()}
#'   \item Determines which species occur at which sites in each time slice
#'   \item Creates site-based co-occurrence matrices using \code{comp_site_cooccurr()}
#'   \item For each species, counts the number of co-occurring species across all sites
#'   \item Optionally removes singleton species (those with zero co-occurrences)
#'   \item Optionally filters by group membership
#' }
#'
#' Co-occurrence calculation details:
#' \itemize{
#'   \item \strong{Self-coexistence}: Automatically excluded from counts
#'   \item \strong{Singleton species}: Species with \code{n.coexistence = 0}
#'         \itemize{
#'           \item If \code{remove.singletons = TRUE}: Excluded from output
#'           \item If \code{remove.singletons = FALSE}: Included with value of 0
#'         }
#'   \item \strong{Site aggregation}: Co-occurrence counts are summed across
#'         all sites where a species occurs in each time slice
#'   \item \strong{Group comparisons}: When using \code{group.focal.compare},
#'         only counts co-occurrences between/within specified groups
#' }
#'
#' This function differs from \code{clade_site_coexistence()} by returning
#' species-level results rather than time slice-level aggregated means. It is
#' useful for analyzing individual species' ecological associations through time
#' or identifying species with consistently high or low co-occurrence patterns.
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
#' # Calculate species-level co-occurrence through time
#' result <- IndivSpec_site_coexistence(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
#'   time.slice = 10
#' )
#'
#' # View results
#' head(result)
#'
#' # Plot co-occurrence patterns for a specific species
#' library(ggplot2)
#' sp1_data <- subset(result, species == "sp1")
#' ggplot(sp1_data, aes(x = time.slice, y = n.coexistence)) +
#'   geom_col() +
#'   labs(x = "Time Slice", y = "Number of Co-occurring Species",
#'        title = "Co-occurrence Pattern for sp1") +
#'   theme_minimal()
#'
#' # Compare co-occurrence between species
#' ggplot(result, aes(x = time.slice, y = n.coexistence, fill = species)) +
#'   geom_col(position = "dodge") +
#'   labs(x = "Time Slice", y = "Number of Co-occurring Species") +
#'   theme_minimal()
#'
#' # Include singleton species in results
#' result_with_singletons <- IndivSpec_site_coexistence(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
#'   time.slice = 10,
#'   remove.singletons = FALSE
#' )
#'
#' # Identify singleton species
#' singletons <- subset(result_with_singletons, n.coexistence == 0)
#' singletons
#'
#' # Calculate co-occurrence between groups
#' result_between <- IndivSpec_site_coexistence(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
#'   time.slice = 10,
#'   group = "group",
#'   group.focal.compare = c("A", "B"),
#'   type.comparison = "between"
#' )
#'
#' # Calculate co-occurrence within a single group
#' result_within <- IndivSpec_site_coexistence(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
#'   time.slice = 10,
#'   group = "group",
#'   group.focal.compare = c("A", "B"),
#'   type.comparison = "within"
#' )
#'
#' # Summary statistics: mean co-occurrence per species across time
#' library(dplyr)
#' species_summary <- result %>%
#'   group_by(species) %>%
#'   summarise(
#'     mean_coexistence = mean(n.coexistence),
#'     sd_coexistence = sd(n.coexistence),
#'     max_coexistence = max(n.coexistence)
#'   )
#' }
IndivSpec_site_coexistence <-
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
      aux_matrix_regional_coex(df.TS.TE, time.slice,
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

    # organizig result in a data frame

    list_n_coex_all3 <- list_n_coex_all2[ lengths(list_n_coex_all2) > 0 ] # removing slices with no species coexistence in sites

    # organizing in a data frame
    df_long_res <-
      do.call(rbind, lapply(names(list_n_coex_all3), function(nm) {
        data.frame(
          time.slice = as.numeric(nm),
          species = names(list_n_coex_all3[[nm]]),
          n.coexistence = unname(list_n_coex_all3[[nm]])
        )
      })
      )

    if(remove.singletons != TRUE){
      df_long_res2 <- df_long_res
    } else{
      df_long_res2  <-
        df_long_res |>
        filter(n.coexistence != 0)
    }

    return(df_long_res2)

  }
