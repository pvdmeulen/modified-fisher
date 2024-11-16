#' Find actual size of MFET ---------------------------------------------------
#'
#' Find the actual size of the MFET as a function of gamma0 by maximising over
#' a nuisance parameter. Uses either a trust-region method, or a zoom method.
#'
#' @param .c Gamma0 to test.
#' @param .odds_ratio The null hypothesis odds ratio being tested. No default.
#' @param .m Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .n Integer input responses and sample sizes. Tests u/m versus v/n. No default.
#' @param .df Testing frame (data frame) generated as part as construct_test_frame().
#' @param .alpha The nominal significance level α. Defaults to 0.05.
#' @param .precision Defines the precision by which confidence limits, p-values, and size is determined. Defaults to 1E-03.
#' @param .method Defines the numerical method used to find the optimum nuisance parameter (maximising actual size) of the test. The default is "zoom", with the second option being "trust" (this uses the trust() function from the trust region package).
#' @param .maze Number of points at each iteration to select the nuisance parameter with maximum size from.
#' @param .zoom_iter Number of iterations to zoom in with (in the "zoom" method).
#'
#' @importFrom numDeriv hessian
#' @importFrom trust trust
#' @keywords find actual size modified fisher exact test trust zoom

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

    # Define objective function that gives local size value, derivative,
    # and second derivative (Hessian):

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

    # Use this + the trust-region method to find the maximum size as a function
    # of the nuisance parameter:

    result <- trust::trust(function(x) objfun(x, .gamma0 = gamma0, .odds_ratio,
                                              .m, .n, .df, .alpha, .precision),
                           start, .precision, 1/.precision, minimize = FALSE)


    size <- result$value

    # SAS code has following, but not needed here - we get the
    # resulting value (size):

    #xopt <- result$value
    #size <- local_size(xopt, .gamma0 = gamma0, .odds_ratio,
    #                   .m, .n, .df, .alpha, .precision)

  } else {

    # Else, use the zoom method as per SAS macro:

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

      i <- which.max(res)

      old <- maxat
      maxat <- old + (i-.maze/2)/(.maze^j)

    }

    size <- max(res)

  } # End of zoom method

  return(size)

} # End of function
