#' Define MFET ----------------------------------------------------------------
#'
#' Define the non-randomised non-conservative MFET. This is a function whether
#' to reject or not based on some value gamma0. Requires n, m, alpha, df,
#' gamma1, gamma2 and 'gamma0'.
#'
#' @param z Vector containing (u, t).
#' @param .odds_ratio The null hypothesis odds ratio being tested. No default.
#' @param .m Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .n Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .df Testing frame (data frame) generated as part as construct_test_frame().
#' @param .gamma0 Some randomisation probability gamma0.
#' @param .alpha The nominal significance level α. Defaults to 0.05.
#' @param .precision Defines the precision by which confidence limits, p-values, and size is determined. Defaults to 1E-03.
#'
#' @keywords modified fisher exact test non conservative randomised randomized

mod_fe_test <- function(z, .df, .gamma0, .odds_ratio, .m, .n,
                        .alpha, .precision){

  #df <- construct_test_frame(.odds_ratio, .m, .n, .alpha, .precision)

  s <- z[[1]]
  i <- z[[2]]

  if(s > .df$c1[[i+1]] & s < .df$c2[[i+1]]){ opttest <- 0 }

  if(s == .df$c1[[i+1]]){ opttest <- .df$gamma1[[i+1]] > .gamma0}
  if(s == .df$c2[[i+1]]){ opttest <- .df$gamma2[[i+1]] > .gamma0}

  if(s < .df$c1[[i+1]] | s > .df$c2[[i+1]]){ opttest <- 1 }

  if(.df$c1[[i+1]] == .df$c2[[i+1]] & s == .df$c1[[i+1]]){
    opttest <- (.df$gamma1[[i+1]]+.df$gamma2[[i+1]]) > .gamma0
  }

  return(opttest)

}
