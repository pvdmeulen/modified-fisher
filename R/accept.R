#' Accept ---- helper function ------------------------------------------------
#'
#' Determine if H : theta = theta0 is accepted based on MFET. Used within
#' modified fisher exaxt test function and not exported.
#'
#' @param z Vector containing (u, t).
#' @param .odds_ratio The null hypothesis odds ratio being tested. No default.
#' @param .m Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .n Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .df Testing frame (data frame) generated as part as construct_test_frame().
#' @param .alpha The nominal significance level α. Defaults to 0.05.
#' @param .precision Defines the precision by which confidence limits, p-values, and size is determined. Defaults to 1E-03.
#' @param .method Defines the numerical method used to find the optimum nuisance parameter (maximising actual size) of the test. The default is "zoom", with the second option being "trust" (this uses the trust() function from the trust region package).
#' @param .maze Number of points at each iteration to select the nuisance parameter with maximum size from.
#' @param .zoom_iter Number of iterations to zoom in with (in the "zoom" method).
#'
#' @keywords accept reject hypothesis modified fisher exact test

accept <- function(z, .odds_ratio, .m, .n, .df, .alpha, .precision,
                   .method, .maze, .zoom_iter){

  gamma0 <- optimise_gamma0(.odds_ratio, .m, .n, .alpha, .precision,
                            .method, .maze, .zoom_iter)

  accept <- 1 - mod_fe_test(z, .df, gamma0, .odds_ratio, .m, .n,
                            .alpha, .precision)

  return(accept)

}
