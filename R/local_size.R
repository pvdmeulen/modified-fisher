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

  for(i in 0:.m){

    for(j in 0:.n){

      x <- c(i, i+j)

      if(i == 0){ term1 <- -.m*(1-p0)^(.m-1) }
      if(i == .m){ term1 <- .m*p0^(.m-1) }

      if(i > 0 & i < .m){
        # x = r of successes, prob, size = number of trials
        term1 <- dbinom(x = i-1, prob = p0, size = .m-2)*.m*(.m-1)/
          (i*(.m-i))*(i-.m*p0)

      }

      if(j == 0){ term2 <- -.n*(1-p1)^(.n-1) }
      if(j == .n){ term2 <- .n*p1^(.n-1) }

      if(j > 0 & j < .n){
        term2 <- dbinom(x = j-1, prob = p1, size = .n-2)*.n*(.n-1)/
          (j*(.n-j))*(j-.n*p1)
      }

      total <- term1*dbinom(x = j, prob = p1, size = .n) +
        term2*dbinom(x = i, prob = p0, size = .m)*dp1dp0

      locsizegrad <- locsizegrad + mod_fe_test(x, .df, .gamma0, .odds_ratio,
                                               .m, .n, .alpha,
                                               .precision)*total

    }

  }

  return(locsizegrad)

}

## Create second derivative (Hessian) for use in trust method -----------------

# local_size_hessian <- function(nuisance, .gamma0, .odds_ratio, .m, .n,
#                                .df, .alpha, .precision){
#
#   p0 <- nuisance[[1]]
#   p1 <- p0/(p0 + .odds_ratio*(1-p0))
#
#   locsizehessian <- 0
#
#   d2p0_d2p1
#   d2p1_d2p0
#
#   d2p0_dp1dp0
#   d2p0_dp0dp1
#
#   d2p1_dp1dp0
#   d2p1_dp0dp1
#
#   for(i in 0:.m){
#
#     for(j in 0:.n){
#
#       x <- c(i, i+j)
#
#       if(i == 0){ term1 <- -.m*(1-p0)^(.m-1) }
#       if(i == .m){ term1 <- .m*p0^(.m-1) }
#
#       if(i > 0 & i < .m){
#         # x = r of successes, prob, size = number of trials
#         term1 <- dbinom(x = i-1, prob = p0, size = .m-2)*.m*(.m-1)/
#           (i*(.m-i))*(i-.m*p0)
#
#       }
#
#       if(j == 0){ term2 <- -.n*(1-p1)^(.n-1) }
#       if(j == .n){ term2 <- .n*p1^(.n-1) }
#
#       if(j > 0 & j < .n){
#         term2 <- dbinom(x = j-1, prob = p1, size = .n-2)*.n*(.n-1)/
#           (j*(.n-j))*(j-.n*p1)
#       }
#
#       total <- term1*dbinom(x = j, prob = p1, size = .n) +
#         term2*dbinom(x = i, prob = p0, size = .m)*dp1dp0
#
#       locsizehessian <- locsizehessian + mod_fe_test(x, .df, .gamma0, .odds_ratio,
#                                                .m, .n, .alpha,
#                                                .precision)*total
#
#     }
#
#   }
#
#   return(locsizehessian)
#
# }
