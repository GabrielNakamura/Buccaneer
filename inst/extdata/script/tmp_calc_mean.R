
  for(i in 1:length(list_matrix_cooccur_site)){
    # i = 1
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
  }
