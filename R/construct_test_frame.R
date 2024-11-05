# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# CONSTRUCT TESTING DATASET ===================================================
# /////////////////////////////////////////////////////////////////////////////

# Also find gamma1 and gamma2 randomisation values.

construct_test_frame <- function(.odds_ratio, .m, .n, .alpha, .precision,
                                 .message = FALSE){

  # Set OR to arbitrarily small number if null is 0 (avoid errors)
  if(.odds_ratio < .precision){
    message(paste0("Odds ratio set to ", .precision, " for testing purposes"))
    or <- .precision
  } else { or <- .odds_ratio }

  # Test:
  #or <- 2

  df <- data.frame(
    "t" = rep(0, .m+.n+1),
    "c1" = rep(0, .m+.n+1),
    "d1" = rep(0, .m+.n+1),
    "gamma1" = rep(0, .m+.n+1),
    "c2" = rep(0, .m+.n+1),
    "d2" = rep(0, .m+.n+1),
    "gamma2" = rep(0, .m+.n+1)
  )

  # Loop version:

  # Test:
  #s <- 8

  for(s in 0:(.m+.n)){

    #message(paste0("Starting on s = ", s))
    df$t[[s+1]] <- s

    # If no successes, set gamma to chosen (two sided) .alpha:

    if(s == 0){

      #df$c1[[s+1]] <- s
      #df$c2[[s+1]] <- s

      df$gamma1[[s+1]] <- .alpha/2
      df$gamma2[[s+1]] <- .alpha/2

      if(.message){
        message(
          paste0("Solution for t = ", s, " found: ", .alpha/2, ", ", .alpha/2)
        )
      }

    }

    # If s > 0, set critical values to the relevant quantiles and find gamma
    # which mixes so that total probability mass is still 1-.alpha:

    if(s > 0){

      # Find critical values that are (conservatively) <= .alpha/2:

      d1 <- BiasedUrn::qFNCHypergeo(n = s, p = .alpha/2, m1 = .m, m2 = .n,
                                    odds = or, precision = .precision,
                                    lower.tail = TRUE)

      d2 <- BiasedUrn::qFNCHypergeo(n = s, p = .alpha/2, m1 = .m, m2 = .n,
                                    odds = or, precision = .precision,
                                    lower.tail = FALSE)

      df$d1[[s]] <- d1
      df$d2[[s]] <- d2

      # If total successes are 1, critical values are 0 and 1 and gamma is .alpha:

      # if(s == 1){
      #
      #   #df$c1[[s+1]] <- s-1
      #   #df$c2[[s+1]] <- s
      #
      #   #df$d1[[s]] <- 0
      #   #df$d2[[s]] <- 1
      #
      #   df$gamma1[[s+1]] <- .alpha
      #   df$gamma2[[s+1]] <- .alpha
      #
      #   if(message){
      #     message(
      #       paste0("Solution for t = ", s, " found: ", .alpha, ", ", .alpha)
      #     )
      #   }
      #
      # }

      # Find randomisation probabilities such that probability mass stays the
      # same:

      correct <- FALSE

      c1 <- d1
      c2 <- d2

      gamma <- find_gamma(.t = s, .odds_ratio = or, c1 = c1, c2 = c2,
                          .m = .m, .n = .n, .alpha = .alpha,
                          .precision = .precision)

      # Test if gamma is 0 < gamma < 1:

      correct <- (0 <= gamma[[1]] & gamma[[1]] <= 1 &
                    0 <= gamma[[2]] & gamma[[2]] <= 1)

      # If solution is admissible:

      if(correct){

        if(.message){
          message(
            paste0("Solution for t = ", s," found: ",
                   paste(round(gamma, 3), collapse = ", "))
          )
        }

        df$c1[[s+1]] <- c1
        df$c2[[s+1]] <- c2

        df$gamma1[[s+1]] <- round(gamma[[1]], digits = 7)
        df$gamma2[[s+1]] <- round(gamma[[2]], digits = 7)

        # If .alpha/2 quantiles are still both zero, set
        # gamma to .alpha:

        # } else if(correct == FALSE & d1 == d2){
        #
        #   df$c1[[s+1]] <- d1
        #   df$c2[[s+1]] <- d2+1
        #
        #   df$gamma1[[s+1]] <- .alpha
        #   df$gamma2[[s+1]] <- .alpha
        #
        #   if(message){
        #     message(
        #       paste0("Solution for t = ", s, " found: ", .alpha, ", ", .alpha)
        #     )
        #   }
        #
        #   correct <- TRUE

      } else {

        if(.message){message(paste0("Solution for t = ", s,
                                   " not found, starting iteration"))}

        # If not, search for one that is by looking around quantile square-wise
        # and expanding the square (kxk) until a solution is found:

        k <- 1

        while(correct == FALSE){

          for(i in max(0, d1-k):min(d1+k, min(s, .m))){

            for(j in max(0, d2-k):min(d2+k, min(s, .m))){

              if(max(d1-i, d2-j) == k | min(d1-i, d2-j) == -k){

                c1 <- i
                c2 <- j

                gamma <- find_gamma(.t = s, .odds_ratio = or,
                                    c1 = c1, c2 = c2, .m = .m, .n = .n,
                                    .alpha = .alpha, .precision = .precision)

                check <- (0 <= gamma[[1]] & gamma[[1]] <= 1 &
                            0 <= gamma[[2]] & gamma[[2]] <= 1)

                if(check == TRUE){

                  correct <- TRUE

                  df$c1[[s+1]] <- c1
                  df$c2[[s+1]] <- c2

                  df$gamma1[[s+1]] <- round(gamma[[1]], digits = 7)
                  df$gamma2[[s+1]] <- round(gamma[[2]], digits = 7)

                  if(.message){
                    message(paste0("Solution for t = ", s,
                                   " found at iteration ",
                                   k, ": ", paste(round(gamma, 3), collapse = ", "))
                    )
                  }

                } # End of check == TRUE

              } # End of min/max == k

            } # End of for j loop

          } # End of for i loop

          k <- k+1

          if(k == 1e4){
            message(paste0("Solution for t = ", s,
                           " not found, stopped iterating at k = ", k))
            correct <- TRUE
          }

        } # End of while correct == FALSE

        k <- 0

      } # End of iteration loop

    } # End of s > 1

  } # End of s loop

  df$d1[[.m+.n+1]] <- .m
  df$d2[[.m+.n+1]] <- .m

  return(df)

} # End of function
