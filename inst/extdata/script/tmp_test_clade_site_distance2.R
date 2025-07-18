
df.TS.TE = df_TS_TE_mass
df.occ = df_occ_faurby
time.slice = 0.1
round.digits = 1
species = "species"
TS = "TS"
TE = "TE"
dist.trait = dist(df_TS_TE_mass$mean.size)
nearest.taxon = 1
trait = NULL
group = "group"
# group.focal.compare = c("Caniformia", "Feliformia")
group.focal.compare = NULL
type.comparison = NULL
Max.age = "Max.age"
Min.age = "Min.age"
site = "site"

clade_site_distance <-
  function(df.TS.TE,
           df.occ,
           time.slice,
           dist.trait,
           nearest.taxon,
           group = NULL,
           group.focal.compare = NULL,
           type.comparison = NULL,
           trait = NULL,
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE",
           Max.age = "Max.age",
           Min.age = "Min.age",
           site = "site"){
    # subseting columns
    if(!is.null(group) == TRUE){
      df.TS.TE <- df.TS.TE[, c(species, trait, TS, TE, group)]
      if(is.null(trait) == TRUE){
        colnames(df.TS.TE) <- c("species", "TS", "TE", "group")
      } else{
        colnames(df.TS.TE) <- c("species", "trait", "TS", "TE", "group")
      }
    } else{
      df.TS.TE <- df.TS.TE[, c(species, trait, TS, TE)]
      if(is.null(trait) == TRUE){
        colnames(df.TS.TE) <- c("species", "TS", "TE")
      } else{
        colnames(df.TS.TE) <- c("species", "trait", "TS", "TE")
      }
    }

    df_occ <-
      df.occ[, c(species, Max.age, Min.age, site)]
    vars <- list(species, Max.age, Min.age, site)
    name_vars <- c("species", "Max.age", "Min.age", "site")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df_occ) <- column.names

    # Generating time intervals used to compute temporal coexistence
    seq_interval <- seq(from = max(df.TS.TE[, "TS"]), to = min(df.TS.TE[, "TE"]), by = -time.slice)
    seq_interval <- c(round(seq_interval, digits = round.digits))

    # coexistence matrix
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE, time.slice, round.digits = 1,
                               species = "species",
                               TS = "TS",
                               TE = "TE")

    # species composition at each timeslice
    spp_slice <-
      lapply(matrix_coex, function(x){
        names(which(rowSums(x) >= 1))
      })

    names(spp_slice) <- seq_interval

    # calculating trait distances for all clades

    # trait distance
    if(!is.null(dist.trait) == TRUE){
      matrix_dist_trait <- as.matrix(dist.trait)
    } else{
      matrix_dist_trait <-
        as.matrix(dist(x = df.TS.TE[, "trait"], method = "euclidean", upper = T, diag = T))
    }
    rownames(matrix_dist_trait) <- df.TS.TE$species
    colnames(matrix_dist_trait) <- df.TS.TE$species


    # modified trait matrix containing group comparison
    if(!is.null(group.focal.compare) == TRUE){
      focal <- group.focal.compare[1]
      compare <- group.focal.compare[2]
      spp_focal <- df.TS.TE[which(df.TS.TE$group == focal), "species"]$species
      spp_compare <- df.TS.TE[which(df.TS.TE$group == compare), "species"]$species

      if(type.comparison == "between"){# comparison between species of two groups
        matrix_dist_trait_comp <- matrix_dist_trait[spp_focal, spp_compare] # focal speices in lines and comparison in columns
      }
      if(type.comparison == "within"){ # comparison only within the focal group
        matrix_dist_trait_comp <- matrix_dist_trait[spp_focal, spp_focal]
      }
    } else{
      matrix_dist_trait_comp <- matrix_dist_trait
    }


    # calculating matrix of species cooccurrence for site

    list_matrix_cooccur_site <-
      comp_site_cooccurr(spp_slice = spp_slice, df.occ = df_occ)


    # Ensure same species and same order in both cooccurrence matrix and distance matrix
    list_dist_spp <-
      lapply(list_matrix_cooccur_site, function(x){
        species_row <- intersect(rownames(x), rownames(matrix_dist_trait_comp))
        species_col <- intersect(colnames(x), colnames(matrix_dist_trait_comp))
        cooccur_matrix <- x[species_row, species_col]
        dist_matrix <- matrix_dist_trait_comp[species_row, species_col]

        # Create matrix to store the results of mean pairwise distances
        mean_distances <- matrix(NA, nrow = length(species_row), ncol = 1,
                                 dimnames = list(species_row, paste("mean.dist.to.cooccur", nearest.taxon, sep = ".")))

        # calculating distances for all species
        for (sp in species_row) {

          #checking if there are no species in the slice or
          if(length(species_row) == 0 | length(species_col) == 0){
            mean_distances[sp, 1] <- NA
            if(length(species_row) != 0){
              mean_distances[sp, 1] <- "NA_singleton"
            }

          } else{
            if(length(species_row) == 1 & length(species_col) == 1){
              mean_distances[sp, 1] <- dist_matrix
            } else{
              cooccur_species <- names(which(cooccur_matrix[sp, ] > 0 & names(cooccur_matrix[sp, ]) != sp))

              # If there are co-occurring species, compute mean distance
              if (length(cooccur_species) > 0) {
                distance_sorted <- sort(dist_matrix[sp, cooccur_species], decreasing = FALSE)
                if(nearest.taxon == "all"){ # calculating for all taxon
                  mean_distances[sp, 1] <- mean(as.numeric(dist_matrix[sp, cooccur_species]), na.rm = TRUE)
                } else{ # using the threshold distance set by the user
                  mean_distances[sp, 1] <- mean(distance_sorted[1:nearest.taxon], na.rm = TRUE)
                }
              }
            }
          }
        }
        return(mean_distances)
      })

    mean_dist_timeslice <- lapply(list_dist_spp, function(x) mean(x, na.rm = TRUE))
    var_dist_timeslice <- lapply(list_dist_spp, function(x) var(x, na.rm = TRUE))

    df_res <-
      data.frame(mean.distance = unlist(mean_dist_timeslice),
                 var.distance = unlist(var_dist_timeslice),
                 time.slice = seq_interval)

    return(df_res)

  }
