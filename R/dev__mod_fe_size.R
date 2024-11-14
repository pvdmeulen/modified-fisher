# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# FIND ACTUAL SIZE OF MOD FE TEST =============================================
# /////////////////////////////////////////////////////////////////////////////

# As a function of gamma0 by maximising over a nuisance parameter.

mfet_size <- function(.c, .odds_ratio, .m, .n, .df, .alpha, .precision,
                      .method, .maze, .zoom_iter){

  if(.method == "trust"){

    size_old <- 0
    start <- 0

    gamma0 <- .c[[1]]

    for(i in 1:(.maze-1)){

      x0 <- i/.maze

      size_new <- local_size(x0, gamma0, .odds_ratio, .m, .n,
                             .df, .alpha, .precision)

      if(size_new > size_old){

        size_old <- size_new

        start <- x0

      }

    }

    # Optimise:
    #optn <- c(1, 0)
    #con <- c(0, 1)
    #gamma0 <- .c[[1]]

    objfun <- function(x, ...){

      # Constrain x:
      if(x > 1 | x < 0) return(list(value = -Inf))

      # Size (value, and function):
      f <- local_size(nuisance = x, ...)

      fx <- function(x) local_size(nuisance = x, ...)

      # Derivative (value, and function):
      g <- local_size_gradient(nuisance = x, ...)

      gx <- function(x) local_size_gradient(nuisance = x, ...)

      # Hessian (= Jacobian of first deriv):
      #B <- numDeriv::jacobian(func = gx, x)
      # Hessian :
      B <- numDeriv::hessian(func = fx, x)

      return(list(value = f, gradient = g, hessian = B))

    }

    result <- trust::trust(function(x) objfun(x, .gamma0 = gamma0, .odds_ratio,
                                              .m, .n, .df, .alpha, .precision),
                         start, precision, 1/precision, minimize = FALSE)

    #xopt <- result$value
    size <- result$value

    #print(paste0("gradient: ", round(result$gradient, digits = 4)))
    #print(paste0("hessian: ", round(result$hessian, digits = 4)))

    #xopt <- nloptr::nloptr(x0 = optn,
    #                       eval_f = ~local_size(x0, .gamma0 = gamma0, .m, .n,
    #                                            .df, .alpha, .precision),
    #                       eval_grad_f = ~local_size_gradient(x0,
    #                                                          .gamma0 = gamma0,
    #                                                          .m, .n, .df,
    #                                                          .alpha,
    #                                                          .precision),
    #                       lb = gamma0[[1]], ub = gamma0[[2]])

    # Try without the following:

    #size <- local_size(xopt, .gamma0 = gamma0, .odds_ratio,
    #                   .m, .n, .df, .alpha, .precision)

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

  } # End of zoom method

  return(size)

} # End of function
