#' Title
#'
#' @param samp
#' @param dis
#' @param null.model
#' @param abundance.weighted
#' @param runs
#' @param iterations
#'
#' @return
#' @export
#'
#' @examples
ses.mpd.modif <-
  function (samp, dis, null.model = c("taxa.labels", "richness",
                                      "frequency", "sample.pool", "phylogeny.pool", "independentswap",
                                      "trialswap"), abundance.weighted = FALSE, runs = 999, iterations = 1000)
  {
    dis <- as.matrix(dis)
    mpd.obs <- picante::mpd(samp, dis, abundance.weighted = abundance.weighted)
    null.model <- match.arg(null.model)
    mpd.rand <- switch(null.model, taxa.labels = t(replicate(runs,
                                                             picante::mpd(samp, picante::taxaShuffle(dis), abundance.weighted = abundance.weighted))),
                       richness = t(replicate(runs, picante::mpd(randomizeMatrix(samp,
                                                                        null.model = "richness"), dis, abundance.weighted))),
                       frequency = t(replicate(runs, picante::mpd(randomizeMatrix(samp,
                                                                         null.model = "frequency"), dis, abundance.weighted))),
                       sample.pool = t(replicate(runs, picante::mpd(randomizeMatrix(samp,
                                                                           null.model = "richness"), dis, abundance.weighted))),
                       phylogeny.pool = t(replicate(runs, picante::mpd(randomizeMatrix(samp,
                                                                              null.model = "richness"), picante::taxaShuffle(dis), abundance.weighted))),
                       independentswap = t(replicate(runs, picante::mpd(randomizeMatrix(samp,
                                                                               null.model = "independentswap", iterations), dis,
                                                               abundance.weighted))), trialswap = t(replicate(runs,
                                                                                                              picante::mpd(randomizeMatrix(samp, null.model = "trialswap",
                                                                                                                                  iterations), dis, abundance.weighted))))
    mpd.rand.mean <- apply(X = t(mpd.rand), MARGIN = 2, FUN = mean,
                           na.rm = TRUE)
    mpd.rand.sd <- apply(X = t(mpd.rand), MARGIN = 2, FUN = sd,
                         na.rm = TRUE)
    mpd.obs.z <- (mpd.obs - mpd.rand.mean)/mpd.rand.sd
    mpd.obs.rank <- apply(X = rbind(mpd.obs, mpd.rand), MARGIN = 2,
                          FUN = rank)[1, ]
    mpd.obs.rank <- ifelse(is.na(mpd.rand.mean), NA, mpd.obs.rank)
    data.frame(ntaxa = vegan::specnumber(samp), mpd.obs, mpd.rand.mean,
               mpd.rand.sd, mpd.obs.rank, mpd.obs.z, mpd.obs.p = mpd.obs.rank/(runs +
                                                                                 1), runs = runs, row.names = row.names(samp))
  }
