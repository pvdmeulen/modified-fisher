## Calculate expected value ---------------------------------------------------

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
