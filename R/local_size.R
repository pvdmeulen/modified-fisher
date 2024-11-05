# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# LOCAL SIZE AND SIZE GRADIENT OF TESTS =======================================
# /////////////////////////////////////////////////////////////////////////////

# Expand in future to cover other tests too / any test fn that outputs (0, 1).

## Define the size of the modified FE test as a function of null OR -----------

# For some fixed p0 nuisance parameter

local_size <- function(nuisance, .gamma0, .odds_ratio, .m, .n,
                       .df, .alpha, .precision){

  p0 <- min(max(0, nuisance[[1]]), 1)
  p1 <- p0/(p0 + .odds_ratio*(1-p0))

  locsize <- 0

  for(i in 0:.m){

    for(j in 0:.n){

      x <- c(i, i+j)

      locsize <- locsize + mod_fe_test(x, .df, .gamma0, .odds_ratio, .m, .n,
                                       .alpha, .precision)*
        dbinom(x = i, prob = p0, size = .m)*
        dbinom(x = j, prob = p1, size = .n)

    }

  }

  return(locsize)

}

## Create local size gradient -------------------------------------------------

local_size_gradient <- function(nuisance, .gamma0, .odds_ratio, .m, .n,
                                .df, .alpha, .precision){

  p0 <- nuisance[[1]]
  p1 <- p0/(p0 + .odds_ratio*(1-p0))

  locsizegrad <- 0
  dp1dp0 <- .odds_ratio/(p0+.odds_ratio*(1-p0))^2

  for(u in 0:.m){

    for(v in 0:.n){

      x <- c(u, u+v)

      if(u == 0){ term1 <- -.m*(1-p0)^(.m-1) }
      if(u == .m){ term1 <- .m*p0^(.m-1) }

      if(u > 0 & u < .m){
        # x = r of successes, prob, size = number of trials
        term1 <- dbinom(x = u-1, prob = p0, size = .m-2)*.m*(.m-1)/
          (u*(.m-u))*(u-.m*p0)
      }

      if(v == 0){ term2 <- -.n*(1-p1)^(.n-1) }
      if(v == .n){ term2 <- .n*p1^(.n-1) }

      if(v > 0 & v < .n){
        term2 <- dbinom(x = v-1, prob = p1, size = .n-2)*.n*(.n-1)/
          (v*(.n-v))*(v-.n*p1)
      }

      total <- term1*dbinom(x = v, prob = p1, size = .n) +
        term2*dbinom(x = u, prob = p0, size = .m)*dp1dp0

      locsizegrad <- locsizegrad + mod_fe_test(x, .df, .gamma0, .odds_ratio,
                                               .m, .n, .alpha,
                                               .precision)*total

    }

  }

  return(locsizegrad)

}
