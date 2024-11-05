## Testing dataset ------------------------------------------------------------

u <- 5
m <- 12

v <- 7
n <- 11

#a <- max(u, 0.5)
#b <- max(m-u, 0.5)
#c <- max(v, 0.5)
#d <- max(n-v, 0.5)

t <- u+v

# RECREATE SAS PROGRAMME ======================================================

# require:
# biasurn
# nlopt

require(BiasedUrn)
require(nloptr)
require(trust) # investigate this package for nloptr option

alpha <- 0.05
precision <- 1e-03
odds_ratio <- 1
message <- FALSE



maze <- 10 # Number of points at each iteration to select nuisance parameter with max size from
method <- "zoom" #c("zoom", "nlopt")
zoom_iter <- 6 # Number of iterations used for the zoom optimisation method.
superiority <- FALSE
power <- TRUE
power_at_pi1 <- 0.5
power_at_pi2 <- 0.75

conf_int <- TRUE
pvalue <- TRUE
plot_size <- TRUE

# Data inputs:

# u -- observed a (# of successes for column 1)
# m -- observed column sum 1
# v -- observed b (# of successes for column 2)
# n -- observed column sum 2

## Check size of input table (original test) --------------------------------

# Checks here

# If size of table is > 2x2, there's an option to simulate p-value
# (see original function)

# If table is 2x2, proceed:

RESULT_pval <- NULL

DATANAME <- deparse(substitute(data))
METHOD <- paste0("Non-Conservative Size-α Modified Fisher's Exact Test for Count Data")

## Input values -------------------------------------------------------------

## Create support - which tables do we need to calculate probability for ----

# Calculate all possible values for a, given the row and column totals:

#lower <- max(.t-.n, 0)
#upper <- min(.m, .t)
#support <- lower:upper

# HELPER FUNCTIONS ------------------------------------------------------------

## Calculate expected value ---------------------------------------------------

calc_expected_value <- function(.odds_ratio, .m, .n, .t, .precision){

  lower <- max(.t-.n, 0)
  upper <- min(.m, .t)
  support <- lower:upper

  # If odds ratios are at extreme ends, return min/max a:

  if(.odds_ratio == 0)
    return(lower)

  if(.odds_ratio == Inf)
    return(upper)

  # Else, return the expected value given our null odds ratio and the
  # non-centric hypergeometric distribution:

  exp_val <- support * BiasedUrn::dFNCHypergeo(
    x = support, m1 = .m, m2 = .n, n = .t,
    odds = .odds_ratio, precision = .precision)

  exp_val <- sum(exp_val)

  return(exp_val)

}

## Find randomisation values so that size is alpha ----------------------------

# Require that phi (fisher exact test) is of size alpha given some
# given critical values:

find_gamma <- function(c1, c2, .odds_ratio, .m, .n, .t, .alpha,
                       .precision){

   # TEST vars
   # .odds_ratio = 2
   # c1 = 1
   # c2 = 2
   # .t = 2

  # requires .m = m, .n = n, .t (t), .odds_ratio = odds_ratio,
  # and critical values c1 and c2

  A <- matrix(0, nrow = 2, ncol = 2)
  b <- c(0, 0)
  gamma <- c(-1, -1) # doesn't seem to be used below??

  lower <- max(.t-.n, 0)
  upper <- min(.m, .t)

  if(.odds_ratio == 1){
    exp_value <- .m/(.m+.n)*.t
  } else {
    exp_value <- calc_expected_value(.odds_ratio, .m, .n,
                                     .t, .precision)
  }

  for(s in lower:upper){

    # Test:
    # s = 2

    # If support (a) is leq c1, or geq c2, simply calculate
    # probability (except when OR is zero)

    if(s <= c1 | s >= c2){

      if(.odds_ratio <= .precision){

        if(s < c1){ p <- 1 }
        if(s >= c1){ p <- 0 }

      } else {

        p <- BiasedUrn::dFNCHypergeo(n = .t, x = s, m1 = .m, m2 = .n,
                                     odds = .odds_ratio, precision = .precision)
      }

      if(s == c1){
        A[1, 1] <- A[1, 1] + p
        A[2, 1] <- A[2, 1] + s*p
      }

      if(s<c1 | s>c2){
        b[1] <- b[1] + p
        b[2] <- b[2] + s*p
      }

      if(s == c2){
        A[1, 2] <- A[1, 2] + p
        A[2, 2] <- A[2, 2] + s*p
      }

      c = c(.alpha, .alpha*exp_value)

      # When solution is found, det(A) is not equal to (c.a.) 0:

      if(det(A) <= -1e-07 | det(A) >= 1e-07){

        #solution <- solve(A) %*% (c-b)
        solution <- solve(A, c - b)

      } else {

        # Fail state when solution is not found:
        solution <- c(-1, -1)

      }

      # If this is the final solution, go to the next value up in
      # random_fe_test function.

      # Unless we're at the extreme values of .t (and no randomisation
      # between values can occur, just 'downrating' the last probabilities
      # to meet our particular probability mass .alpha/2):

      if(.t == 0 | .t == (.m + .n)){

        solution <- c(.alpha/2, .alpha/2)

      }

    } # End of s <= c1, s >= c2
  } # End of s for loop

  return(solution)

}

#find_gamma(c1 = 0, c2 = 2, .odds_ratio = 2, .t = 2, .m = m, .n = n,
#              .alpha = alpha, .precision = precision)

## Randomised FE test ---------------------------------------------------------

randomised_fe_test <- function(.odds_ratio, .m, .n, .alpha, .precision,
                               .message = FALSE){

  # Set OR to arbitrarily small number if null is 0 (avoid errors)
  if(.odds_ratio < .precision){ or <- .precision } else { or <- .odds_ratio }

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

      if(message){
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

        if(message){
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

        if(message){message(paste0("Solution for t = ", s,
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

                  if(message){
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

#randomised_fe_test(.odds_ratio = 2, .m = m, .n = n, .alpha = alpha,
#                   .precision = precision, .message = message)

## Define non-randomised non-conservative modified FET ------------------------

# This is a function of z = (u, t), and is basically an indicator
# function whether to reject or not based on some value gamma0.

# Requires n, m, alpha, df$c1,c2,gamma1,gamma2 and 'gamma0'

modified_fe_test <- function(z, .df, .gamma0, .odds_ratio, .m, .n,
                             .alpha, .precision){

  #df <- randomised_fe_test(.odds_ratio, .m, .n, .alpha, .precision)

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

      locsize <- locsize + modified_fe_test(x, .df, .gamma0, .odds_ratio, .m, .n,
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

      locsizegrad <- locsizegrad + modified_fe_test(x, .df, .gamma0, .odds_ratio,
                                                    .m, .n, .alpha,
                                                    .precision)*total

    }

  }

  return(locsizegrad)

}

## Find actual size of of modified FE test ------------------------------------

# As a function of gamma0 by maximising over a nuisance parameter

mfet_size <- function(.c, .odds_ratio, .m, .n, .df, .alpha, .precision){

  if(method == "npltr"){

    size_old <- 0
    start <- 0

    gamma0 <- .c[[1]]

    for(i in 1:(maze-1)){

      x0 <- i/maze

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

  } else {

    res <- rep(0, maze-1)
    pt <- 0
    gamma0 <- .c[[1]]
    maxat <- 0.5

    for(j in 1:zoom_iter){

      for(i in 1:(maze-1)){

        pt <- maxat + (i-maze/2)/(maze^j)

        res[[i]] <- local_size(pt, .gamma0 = gamma0, .odds_ratio, .m, .n, .df,
                               .alpha, .precision)

      }

      #i <- max(which(res == max(res)))
      i <- which.max(res)

      old <- maxat
      maxat <- old + (i-maze/2)/(maze^j)

    }

    size <- max(res)

  }

  return(size)

}

## Find optimal solution and turn into function -------------------------------

optimise_gamma0 <- function(.odds_ratio, .m, .n, .alpha, .precision){

  df <- randomised_fe_test(.odds_ratio, .m, .n, .alpha, .precision)
  gamind <- c(df$gamma1, df$gamma2)

  # Sort gammas:
  gamind <- sort(gamind)

  # Start bisecting on the index of gamind:
  g0 <- 0
  g1 <- 2*(.m+.n+1)
  g <- .m+.n+1

  sz0 <- mfet_size(.c = gamind[[g]], .odds_ratio, .m, .n, .df = df,
                   .alpha, .precision)

  while(as.integer(g1-g0) > 1){

    if(sz0 > .alpha){ g0 <- g } else { g1 <- g }
    if(sz0 < .alpha){

      size_old <- sz0
      # Only keep this old size if it is correct, old and new
      # can both be wrong

    }

    g <- as.integer((g0+g1)/2)
    sz0 <- mfet_size(.c = gamind[[g]], .odds_ratio, .m, .n, .df = df,
                     .alpha, .precision)

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

## Determine if H : theta = theta0 is accepted based on modified test ---------

accept <- function(z, .odds_ratio, .m, .n, .df, .alpha, .precision){

  gamma0 <- optimise_gamma0(.odds_ratio, .m, .n, .alpha, .precision)
  accept <- 1 - modified_fe_test(z, .df, gamma0, .odds_ratio, .m, .n,
                                 .alpha, .precision)

  return(accept)

}

## Power of mod FE test as function of z for null -----------------------------

power_mfet <- function(p, .gamma0, .odds_ratio, .m, .n, .df, .alpha, .precision){

  p0 <- p[[1]]
  p1 <- p[[2]]

  z <- c(0, 0)

  power <- 0

  for(u in 0:.m){

    for(v in 0:.n){

      z <- c(u, u+v)

      if(superiority == TRUE){

        power <- power + (v/.n > u/.m)*
          modified_fe_test(z, .df, .gamma0, .odds_ratio, .m, .n, .alpha, .precision)*
          dbinom(x = u, prob = p0, size = .m)*
          dbinom(x = v, prob = p1, size = .n)

      } else {

        power <- power +
          modified_fe_test(z, .df, .gamma0, .odds_ratio, .m, .n, .alpha, .precision)*
          dbinom(x = u, prob = p0, size = .m)*
          dbinom(x = v, prob = p1, size = .n)

      }

    }

  }

  return(power)

}

# START MAIN FUNCTION =========================================================

# Require:
# m
# n
# alpha
# u
# t (u + v)

# and:
z <- c(u, t)

df <- randomised_fe_test(.odds_ratio = odds_ratio, .m = m, .n = n, .alpha = alpha,
                         .precision = precision)

message <- FALSE

if(power == TRUE){

  or <- odds_ratio

  gamma0 <- optimise_gamma0(.odds_ratio = or, .m = m, .n = n, .alpha = alpha,
                            .precision = precision)

  pi1 <- power_at_pi1
  pi2 <- power_at_pi2

  p <- c(pi1, pi2)

  power_mfet <- power_mfet(p, gamma0, or, .m = m, .n = n, .df = df,
                           .alpha = alpha, .precision = precision)*100

}

# Input values, with +0.5 adjustment if 0:

a <- max(u, 0.5)
b <- max(m-u, 0.5)
c <- max(v, 0.5)
d <- max(n-v, 0.5)

or0 <- (a*d)/(b*c)

z_alpha <- qnorm(alpha/2)

upper0 <- exp(log(or0) + z_alpha*sqrt(1/a+1/b+1/c+1/d))
lower0 <- exp(log(or0) - z_alpha*sqrt(1/a+1/b+1/c+1/d))

if(conf_int){

  # Starting bisecting upper limit

  f1 <- or0 + 0.75*(upper0-or0)
  f2 <- 1.25*upper0
  crit <- abs(f2-f1)

  if(b>0 & c>0){

    while(crit > precision){

      or <- (f1+f2)/2
      answer <- accept(z, .odds_ratio = or, .m = m, .n = n, .df = df,
                       .alpha = alpha, .precision = precision)

      if(answer == 1){ f1 <- or } else { f2 <- or }

      crit <- abs(f2-f1)

    }

    conf_int_upper <- f1

  } else {

    conf_int_upper <- Inf

  }

  # Starting bisecting lower limit:

  f1 <- 0.5*lower0
  f2 <- lower0+0.25*(or0-lower0)
  crit <- abs(f2-f1)

  precision_conf_int <- 0.1

  if(a>0 & d>0){

    while(crit > precision){

      or <- (f1+f2)/2
      answer <- accept(z, .odds_ratio = or, .m = m, .n = n, .df = df,
                       .alpha = alpha, .precision = precision)

      if(answer == 1){ f2 <- or } else { f1 <- or }

      crit <- abs(f2-f1)

    }

    conf_int_lower <- f2

  } else {

    conf_int_lower <- 0

  }

} # End of conf_int == TRUE



if(plot_size == TRUE){

  # PICK UP HERE

  plot_data <- data.frame(
    "pi1" = seq(0, 1, by = 1/100)
  )

  #df <- randomised_fe_test(.odds_ratio = or, .m = m, .n = n, .alpha = alpha,
  #                         .precision = precision)

  gam0 <- optimise_gamma0(.odds_ratio = odds_ratio, .m = m, .n = n,
                          .alpha = alpha, .precision = precision)

  for(row in 1:101){

    point <- plot_data$pi1[[row]]

    plot_data$size[[row]] <- local_size(point, .gamma0 = gam0,
                                        .odds_ratio = odds_ratio,
                                        .m = m, .n = n, .df = df,
                                        .alpha = alpha,
                                        .precision = precision)*100

  }

  plot_data$size <- as.numeric(plot_data$size)
  #plot_data$size <- as.numeric(plot_data$size)

}

plot_data$method <- "modified"

require(ggplot2)
require(latex2exp)

plot_subtitle <- paste0("$m = ", m, ", n = ", n, ", H_0",
                        " : ", odds_ratio, "$")

ggplot(plot_data, aes(x = pi1, y = size)) +
  geom_line(aes(colour = method), linewidth = 1) +
  #scale_y_continuous(limits = c(-4, 4*alpha*100)) +
  geom_hline(yintercept = alpha*100, linetype = 2, colour = "grey30",
             linewidth = 0.4) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Size of Modified Fisher Exact Tests",
    subtitle = latex2exp::TeX(paste0("\\textit{$m=", m, ",\\,n=", n, ",\\,H_{0}=",
                                     odds_ratio, "$}")),
    y = "Size (%)",
    x = bquote(pi[.("1")])
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(family = "STIX Two Text", face = "italic"),
    axis.title.x = element_text(family = "STIX Two Text", face = "italic"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

ggsave(last_plot(), width = 160, height = 120, units = "mm", dpi = "retina",
       filename = paste0("Example/size_plot_or", odds_ratio, "_m", m, "_n", n,
                         ".png"))


p0 <- 0
p1 <- 1
a0 <- 0.5
crit <- p1-p0

# p-value for H0 :

while(abs(crit) > precision){

  reject <- 1-accept(z, .odds_ratio = odds_ratio, .m = m, .n = n, .df = df,
                     .alpha = alpha, .precision = precision)

  if(reject == 1){
    p1 <- a0
    a0 <- (p0+p1)/2
  } else {
    p0 <- a0
    a0 <- (p0+p1)/2
  }

  crit <- p1-p0

  RESULT_pval <- a0
  #label <- paste0("p-value for H0: OR = ", odds_ratio,
  #                round(p_value, digits = nchar(1/precision)))

}

## Create output --------------------------------------------------------------

# Given u, .m, v, .n, or0:

RESULT_conf_int <- if(conf_int) c(conf_int_lower, conf_int_upper)

attr(odds_ratio, "names") <- "odds ratio"
#attr(RESULT_estimate, "names") <- "odds ratio"
attr(RESULT_conf_int, "conf.level") <- 1-alpha
#attr(RESULT_crit, "names")
#attr(RESULT_random, "names")

# Put results into list:

RESULTS <- list(
  p.value = RESULT_pval,
  #estimate = RESULT_estimate,
  conf.int = if(conf_int) RESULT_conf_int,
  null.value = odds_ratio,
  alternative = "two.sided",
  method = METHOD,
  data.name = DATANAME,
  #crit.values = RESULT_crit,
  #random.probs = if(randomise) RESULT_random,
  support.data = df,
  gamma0 = gam0,
  power = power_mfet
)

attr(RESULTS, "class") <- "htest"

RESULTS$gamma0

RESULTS$support.data

