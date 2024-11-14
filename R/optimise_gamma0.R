# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# FIND OPTIMAL GAMMA0 BY MAXIMISING SIZE WRT ALPHA ===========================
# ////////////////////////////////////////////////////////////////////////////

optimise_gamma0 <- function(.odds_ratio, .m, .n, .alpha, .precision,
                            .method, .maze, .zoom_iter){

  df <- construct_test_frame(.odds_ratio, .m, .n, .alpha, .precision)

  gamind <- c(df$gamma1, df$gamma2)

  # Sort gammas:
  gamind <- sort(gamind)

  # Start bisecting on the index of gamind:
  g0 <- 0
  g1 <- 2*(.m+.n+1)
  g <- .m+.n+1

  sz0 <- mfet_size(.c = gamind[[g]], .odds_ratio, .m, .n, .df = df,
                   .alpha, .precision, .method, .maze, .zoom_iter)

  while(as.integer(g1-g0) > 1){

    if(sz0 > .alpha){ g0 <- g } else { g1 <- g }
    if(sz0 < .alpha){

      size_old <- sz0
      # Only keep this old size if it is correct, old and new
      # can both be wrong

    }

    g <- as.integer((g0+g1)/2)
    sz0 <- mfet_size(.c = gamind[[g]], .odds_ratio, .m, .n, .df = df,
                     .alpha, .precision, .method, .maze, .zoom_iter)

    if( g1-g0 == 1 ) {

      if( sz0 > .alpha ){
        g_opt <- g1
        size_opt <- size_old
      } else {
        g_opt <- g
        size_opt <- sz0
      }

      if( g == g0 & sz0 > .alpha ){

        g_opt <- g+1
        size_opt <- size_old

      }

    }

  }

  gamma0 <- gamind[[g_opt]]

  return(gamma0)

}
