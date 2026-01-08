
matrix_coex <-
  aux_matrix_regional_coex(df.TS.TE = df_longevities_canidae,
                           time.slice = 0.1,
                           round.digits = 1,
                           species = "species",
                           TS = "TS",
                           TE = "TE")

# species composition at each timeslice
spp_slice <-
  lapply(matrix_coex, function(x){
    names(which(rowSums(x) >= 1))
  })

seq_interval <- seq(from = ceiling(max(df_longevities_canidae[, "TS"])),
                    to = ceiling(min(df_longevities_canidae[, "TE"])),
                    by = -0.1)

names(spp_slice) <- format(seq_interval, trim = TRUE, scientific = FALSE)


res_obs_assemblage_trait_site <-
  assemblage_site_trait_distance(df.TS.TE = df_longevities_canidae,
                               df.occ = df_occurrence_canidae,
                               time.slice = 0.1,
                               dist.trait = dist_body_mass,
                               nearest.taxon = 1,
                               Max.age = "max_T",
                               Min.age = "min_T",
                               site = "site.char")

list_occurrence <-
  comp_site_occurrence(spp_slice = spp_slice,
                       df.occ = df_occurrence_canidae,
                       species = "species",
                       Max.age = "max_T",
                       Min.age = "min_T",
                       site = "site.char")

dist_body_mass <- as.matrix(dist(traits_canidae$LD1))
rownames(dist_body_mass) <- df_longevities_canidae$species
colnames(dist_body_mass) <- df_longevities_canidae$species

list_occurrence = list_occurrence
names(list_occurrence) <- format(seq_interval, trim = TRUE, scientific = FALSE)
dist_matrix_trait = dist_body_mass
nearest.taxon = 1
nperm = 1000

calc_null_model <-
  function(list_occurrence, dist_matrix_trait, nearest.taxon, nperm = 1000){

    # data frame to receive the results
    list_res_timeslice <-
      vector(mode = "list", length = length(list_occurrence))

    for(i in 1:length(list_occurrence)){
      # i = 12
      # compute null matrix for each timeslice composition
      if(dim(list_occurrence[[i]])[1] <= 1){
        list_res_timeslice[[i]] <- NA
      } else{
        null_occ_site <-
          vegan::permatfull(m = list_occurrence[[i]][, -1],
                            fixedmar = "rows",
                            mtype = "prab",
                            times = nperm)

        # computing mpd
        list_matrix_coocccur_null_site <-
          lapply(null_occ_site$perm, function(x){
            picante::mpd(samp = x,
                         dis = dist_body_mass,
                         abundance.weighted = FALSE)
          })

        # joining all values of null mpd, columns are assemblages, rows are reps, only for one timeslice
        matrix_null_mpd <- do.call(rbind, list_matrix_coocccur_null_site)

        # getting only the observed value for that slice
        dist_obs_slice <-
          res_obs_assemblage_trait_site |>
          filter(time.slice == format(seq_interval[i], trim = TRUE, scientific = FALSE)) |>
          pull(mean_dist_to_cooccur) # the same timeslice used to calculate null matrix

        # mean value for all communities in the slice
        mean_null_slice <- apply(matrix_null_mpd, 2, mean)

        # variance for all communities in slice
        sd_null_slice <- apply(matrix_null_mpd, 2, sd)

        # standardized effect size for all slices - this is ses.mpd
        z_score_slice <- ((dist_obs_slice - mean_null_slice) / sd_null_slice)

        # absolute deviations of null
        dev_null <- abs(sweep(matrix_null_mpd, 2, mean_null_slice, "-"))

        # absolute deviations of observed
        dev_obs <- abs(dist_obs_slice - mean_null_slice)

        # count number of null >= observed for each column
        counts <- colSums(dev_null >= rep(dev_obs, each = nrow(matrix_null_mpd)))

        # two-tailed p-value
        p_two <- (counts + 1) / (nrow(matrix_null_mpd) + 1)

        # joining results for one timeslice in a dataframe
        df_res <-
          data.frame(dist.obs = dist_obs_slice,
                     dist.obs.z = z_score_slice,
                     p.value = p_two,
                     time.slice = names(list_occurrence[i]))

        # list to receive all results

        list_res_timeslice[[i]] <- df_res

      } # end of conditional
      print(i)
    }# end loop for each timeslice


    return(list_res_timeslice)
  }


res_all_mpd_timeslices <- do.call(rbind, list_res_timeslice)

res_all_mpd_timeslices$signif <- res_all_mpd_timeslices$p.value <= 0.05



library(ggplot2)



ggplot(res_all_mpd_timeslices, aes(x = as.numeric(time.slice))) +
  #geom_point(aes(y = dist.obs, color = signif), alpha = 0.3) +
  #geom_smooth(aes(y = dist.obs), method = "loess", se = FALSE) +
  geom_point(aes(y = dist.obs.z, color = signif), alpha = 0.3) +
  geom_smooth(aes(y = dist.obs.z), method = "loess", se = FALSE, color = "blue") +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "black"),
                     na.value = "grey70",
                     labels = c(`TRUE` = "p ≤ 0.05", `FALSE` = "NS")) +
  scale_x_reverse() +                                   # ⬅ IMPORTANT
  labs(y = "MPD / SES(MPD)",
       x = "Time slice (Ma)",
       color = "Significance") +
  theme_minimal()

