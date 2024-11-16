#' Find the power of the MFET -------------------------------------------------
#'
#' Define power of the MFET as a function of z = (u, u+v) for the null OR=1.
#'
#' @param p Vector containing (pi1, pi2).
#' @param .gamma0 Some randomisation probability gamma0.
#' @param .odds_ratio The null hypothesis odds ratio being tested. No default.
#' @param .m Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .n Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .df Testing frame (data frame) generated as part as construct_test_frame().
#' @param .alpha The nominal significance level α. Defaults to 0.05.
#' @param .precision Defines the precision by which confidence limits, p-values, and size is determined. Defaults to 1E-03.
#' @param .superiority A logical. Defaults to FALSE. Setting this to TRUE will calculate the power for testing superiority.
#'
#' @keywords find power test

power_mfet <- function(p, .gamma0, .odds_ratio, .m, .n, .df, .alpha,
                       .precision, .superiority){

  p0 <- p[[1]]
  p1 <- p[[2]]

  z <- c(0, 0)

  power <- 0

  for(u in 0:.m){

    for(v in 0:.n){

      z <- c(u, u+v)

      if(.superiority == TRUE){

        power <- power + (v/.n > u/.m)*
          mod_fe_test(z, .df, .gamma0, .odds_ratio, .m, .n, .alpha, .precision)*
          stats::dbinom(x = u, prob = p0, size = .m)*
          stats::dbinom(x = v, prob = p1, size = .n)

      } else {

        power <- power +
          mod_fe_test(z, .df, .gamma0, .odds_ratio, .m, .n, .alpha, .precision)*
          stats::dbinom(x = u, prob = p0, size = .m)*
          stats::dbinom(x = v, prob = p1, size = .n)

      }

    }

  }

  return(power)

}
