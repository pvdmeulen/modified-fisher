#' Calculate expected OR value ---- helper function ---------------------------
#'
#' Calculate expected value of theta (OR) given m, n, t and precision.
#'
#' @param .odds_ratio The null hypothesis odds ratio being tested. No default.
#' @param .m Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .n Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .t Integer input responses and sample sizes. t = u + v.
#' @param .precision Defines the precision by which confidence limits, p-values, and size is determined. Defaults to 1E-03.
#'
#' @keywords calculate expected value odds ratio theta
#' @importFrom BiasedUrn dFNCHypergeo

calc_expected_value <- function(.odds_ratio, .m, .n, .t, .precision){

  lower <- max(.t-.n, 0)
  upper <- min(.m, .t)
  support <- lower:upper

  # If odds ratios are at extreme ends, return min/max a:

  if(.odds_ratio == 0)
    return(lower)

  if(.odds_ratio == Inf)
    return(upper)

  # Else, return the expected value given our null odds ratio and the
  # non-centric hypergeometric distribution:

  exp_val <- support * BiasedUrn::dFNCHypergeo(
    x = support, m1 = .m, m2 = .n, n = .t,
    odds = .odds_ratio, precision = .precision)

  exp_val <- sum(exp_val)

  return(exp_val)

}
