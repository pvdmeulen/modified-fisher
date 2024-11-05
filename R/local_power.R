# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# POWER OF TEST AS FUNCTION OF NULL ===========================================
# /////////////////////////////////////////////////////////////////////////////

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
          dbinom(x = u, prob = p0, size = .m)*
          dbinom(x = v, prob = p1, size = .n)

      } else {

        power <- power +
          mod_fe_test(z, .df, .gamma0, .odds_ratio, .m, .n, .alpha, .precision)*
          dbinom(x = u, prob = p0, size = .m)*
          dbinom(x = v, prob = p1, size = .n)

      }

    }

  }

  return(power)

}
