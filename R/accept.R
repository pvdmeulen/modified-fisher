## Determine if H : theta = theta0 is accepted based on modified test ---------

accept <- function(z, .odds_ratio, .m, .n, .df, .alpha, .precision,
                   .method, .maze, .zoom_iter){

  gamma0 <- optimise_gamma0(.odds_ratio, .m, .n, .alpha, .precision,
                            .method, .maze, .zoom_iter)

  accept <- 1 - mod_fe_test(z, .df, gamma0, .odds_ratio, .m, .n,
                            .alpha, .precision)

  return(accept)

}
