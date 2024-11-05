# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# DEFINE RANDOMISED CONSERVATIVE FET ==========================================
# /////////////////////////////////////////////////////////////////////////////

# Just for comparison to MFET.

# This is a function of z = (u, t), and is basically an indicator
# function whether to reject or not. Requires n, m, alpha, and df.

random_fe_test <- function(z, .df, .gamma0, .odds_ratio, .m, .n,
                           .alpha, .precision){

  s <- z[[1]]
  i <- z[[2]]

  if(s > .df$c1[[i+1]] & s < .df$c2[[i+1]]){ randtest <- 0 }

  if(s == .df$c1[[i+1]]){ randtest <- .df$gamma1[[i+1]] > .gamma0}
  if(s == .df$c2[[i+1]]){ randtest <- .df$gamma2[[i+1]] > .gamma0}

  if(s < .df$c1[[i+1]] | s > .df$c2[[i+1]]){ randtest <- 1 }

  if(.df$c1[[i+1]] == .df$c2[[i+1]] & s == .df$c1[[i+1]]){
    randtest <- (.df$gamma1[[i+1]]+.df$gamma2[[i+1]]) > .gamma0
  }

  return(randtest)

}
