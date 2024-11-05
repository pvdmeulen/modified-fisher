# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# FIND ACTUAL SIZE OF MOD FE TEST =============================================
# /////////////////////////////////////////////////////////////////////////////

# As a function of gamma0 by maximising over a nuisance parameter.

mfet_size <- function(.c, .odds_ratio, .m, .n, .df, .alpha, .precision,
                      .method, .maze, .zoom_iter){

  if(.method == "npltr"){

    size_old <- 0
    start <- 0

    gamma0 <- .c[[1]]

    for(i in 1:(.maze-1)){

      x0 <- i/.maze

      size_new <- local_size(x0, .gamma0 = gamma0, .odds_ratio, .m, .n,
                             .df, .alpha, .precision)

      if(size_new > size_old){

        size_old <- size_new

        start <- x0

      }

    }

    # Optimise:
    optn <- c(0, 1)
    con <- c(0, 1)
    gamma0 <- .c[[1]]

    xopt <- nloptr::nloptr(x0 = optn,
                           eval_f = ~local_size(x0, .gamma0 = gamma0, .m, .n,
                                                .df, .alpha, .precision),
                           eval_grad_f = ~local_size_gradient(x0,
                                                              .gamma0 = gamma0,
                                                              .m, .n, .df,
                                                              .alpha,
                                                              .precision),
                           lb = gamma0[[1]], ub = gamma0[[2]])

    size <- local_size(xopt, .gamma0 = gamma0, .m, .n, .df, .alpha, .precision)

    # Test trust package here

  } else {

    res <- rep(0, .maze-1)
    pt <- 0
    gamma0 <- .c[[1]]
    maxat <- 0.5

    for(j in 1:.zoom_iter){

      for(i in 1:(.maze-1)){

        pt <- maxat + (i-.maze/2)/(.maze^j)

        res[[i]] <- local_size(pt, .gamma0 = gamma0, .odds_ratio, .m, .n, .df,
                               .alpha, .precision)

      }

      #i <- max(which(res == max(res)))
      i <- which.max(res)

      old <- maxat
      maxat <- old + (i-.maze/2)/(.maze^j)

    }

    size <- max(res)

  }

  return(size)

}
