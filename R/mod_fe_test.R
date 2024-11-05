# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Define non-randomised non-conservative modified FET -------------------------
# /////////////////////////////////////////////////////////////////////////////

# This is a function of z = (u, t), and is basically an indicator
# function whether to reject or not based on some value gamma0.

# Requires n, m, alpha, df$c1,c2,gamma1,gamma2 and 'gamma0'

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
