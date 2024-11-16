#' Define randomised conservative fisher exact test ---------------------------
#'
#' This is a function of z = (u, t), and is an indicator function whether to
#' reject or not. Requires n, m, alpha, and df. Mainly for comparing with MFET.
#'
#' @param z Vector containing (u, t).
#' @param .odds_ratio The null hypothesis odds ratio being tested. No default.
#' @param .m Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .n Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .df Testing frame (data frame) generated as part as construct_test_frame().
#' @param .alpha The nominal significance level α. Defaults to 0.05.
#' @param .precision Defines the precision by which confidence limits, p-values, and size is determined. Defaults to 1E-03.
#'
#' @keywords randomised randomized fisher exact test

random_fe_test <- function(z, .df, .odds_ratio, .m, .n, .alpha, .precision){

  # WIP!!

  randtest <- NULL

  s <- z[[1]]
  i <- z[[2]]

  if(s > .df$c1[[i+1]] & s < .df$c2[[i+1]]){ randtest <- 0 }

  #if(s == .df$c1[[i+1]]){ randtest <- .df$gamma1[[i+1]] > .gamma0}
  #if(s == .df$c2[[i+1]]){ randtest <- .df$gamma2[[i+1]] > .gamma0}

  if(s < .df$c1[[i+1]] | s > .df$c2[[i+1]]){ randtest <- 1 }

  #if(.df$c1[[i+1]] == .df$c2[[i+1]] & s == .df$c1[[i+1]]){
  #  randtest <- (.df$gamma1[[i+1]]+.df$gamma2[[i+1]]) > .gamma0
  #}

  return(randtest)

}
