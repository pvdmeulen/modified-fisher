
# TESTING FUNCTIONS ===========================================================

## Inputs ---------------------------------------------------------------------

# Cell values:

a <- 3
b <- 8
c <- 4
d <- 2

input_data <- matrix(data = c(a, c, b, d), nrow = 2, ncol = 2)

## Fisher's Exact Test --------------------------------------------------------

fisher_exact_test <- function(data, odds_ratio = 1, alpha = 0.05,
                              alternative = "two.sided", conf_int = TRUE,
                              randomise){

  #++ TEST ++#
  # data <- matrix(data = c(3, 4, 8, 2), nrow = 2, ncol = 2)
  data = matrix(c(troubleshoot$test_a,
                  troubleshoot$test_c,
                  troubleshoot$test_b,
                  troubleshoot$test_d), nrow = 2, ncol = 2)
  odds_ratio <- 1
  alpha <- 0.05
  alternative <- "two.sided"
  conf_int <- TRUE
  randomise <- TRUE

  ## Check size of input table (original test) --------------------------------

  # Checks here

  conf.level <- 1-alpha

  # If size of table is > 2x2, there's an option to simulate p-value
  # (see original function)

  # If table is 2x2, proceed:

  RESULT_pval <- NULL

  DATANAME <- deparse(substitute(data))
  METHOD <- "Fisher's Exact Test for Count Data"

  ## Input values -------------------------------------------------------------

  # Row and column sums:

  obs_a <- data[1, 1]
  obs_b <- data[1, 2]
  obs_c <- data[2, 1]
  obs_d <- data[2, 2]

  obs_r1 <- obs_a + obs_b
  obs_r2 <- obs_c + obs_d

  obs_c1 <- obs_a + obs_c
  obs_c2 <- obs_b + obs_d

  obs_n <- obs_r1 + obs_r2

  ## Create support - which tables do we need to calculate probability for ----

  # Calculate all possible values for a, given the row and column totals:

  min_a <- max(obs_r1-obs_c2, 0)
  max_a <- min(obs_c1, obs_r1)

  support <- min_a:max_a

  # In general, the conditional distribution of a given cell value with the
  # marginals r1, r2, c1, and c2 is a non-central hypergeometric distribution
  # with non-centrality parameter being the odds ratio.

  # Calculate the density of the central hypergeometric distribution on its
  # support (x = support, m = c1, n = c2, k (successes) = r1):

  logdc <- dhyper(support, obs_c1, obs_c2, obs_r1, log = TRUE)

  #RESULT_logdc <- c("logged density" = logdc)

  # Logged to make it easier to work with.

  # DEFINE FUNCTIONS ----------------------------------------------------------

  # Move functions to inside this function to avoid defining every argument
  # every time?

  ## Non-central hypergeometric distribution function -------------------------

  # Density of the non-central hypergeometric distribution on its support,
  # using logged OR as non-centrality parameter:

  dn_hypergeo <- function(odds_ratio){

    # Add non-centrality parameter:
    d <- logdc + log(odds_ratio) * support

    # Normalise and turn back to non-logged version:
    d <- exp(d - max(d))

    # Calculate probabilities:
    prob <- d / sum(d)

    return(prob)

    # Does not work for NCP boundary values 0 or Inf
    # (these are edge cases resolved later)

  }

  ## Determine critical values (no randomisation) -----------------------------

  # Generate table with cumulative probabilities:

  prob_df <- data.frame(
    "support" = support,
    "prob" = dn_hypergeo(odds_ratio = odds_ratio)
  )

  prob_df$cumsum1 = cumsum(prob_df$prob)
  prob_df$cumsum2 = rev(cumsum(rev(prob_df$prob)))

  # Using our density function, what are the critical values for a where
  # we reject the null hypothesis of no association (or some other null)?

  critical_values <- function(odds_ratio, alternative){

    # Return the critical values given our null odds ratio and the chosen
    # alternative hypothesis -- NO RANDOMISATION:

    critical_values <- switch(
      alternative,
      less = if(nrow(prob_df[prob_df$cumsum1 <= alpha, ]) == 0){
        c("c1" = min_a, "c2" = NA)
      } else {
        c("c1" = min(prob_df[prob_df$cumsum1 <= alpha, "support"]), "c2" = NA)
      },

      greater = if(nrow(prob_df[prob_df$cumsum2 <= alpha, ]) == 0){
        c("c1" = NA, "c2" = max_a)
      } else {
        c("c1" = NA, "c2" = max(prob_df[prob_df$cumsum2 <= alpha, "support"]))
      },

      two.sided = if(nrow(prob_df[prob_df$cumsum1 <= alpha/2, ]) == 0 &
                     nrow(prob_df[prob_df$cumsum2 <= alpha/2, ]) == 0){

        c("c1" = min_a, "c2" = max_a)

        } else if(nrow(prob_df[prob_df$cumsum1 <= alpha/2, ]) == 0){

        c("c1" = min_a,
          "c2" = max(prob_df[prob_df$cumsum2 <= alpha/2, "support"]))

        } else if(nrow(prob_df[prob_df$cumsum2 <= alpha/2, ]) == 0){

          c("c1" = min(prob_df[prob_df$cumsum1 <= alpha/2, "support"]),
            "c2" = max_a)

        } else {

          c(
            "c1" = max(prob_df[prob_df$cumsum1 <= alpha/2, "support"]),
            "c2" = min(prob_df[prob_df$cumsum2 <= alpha/2, "support"]))
        }
      )

    #   less = min(support[dn_hypergeo(odds_ratio = odds_ratio) <= alpha],
    #              min_a+1),
    #   greater = min(support[dn_hypergeo(odds_ratio = odds_ratio) <= alpha],
    #                 max_a-1),
    #   two.sided = unique(c(
    #     min(support[dn_hypergeo(odds_ratio = odds_ratio) <= alpha/2], min_a+1),
    #     min(support[dn_hypergeo(odds_ratio = odds_ratio) <= alpha/2], max_a-1)
    #     )
    #     )
    # )

    return(critical_values)

  }

  # Add crit values to prob_df:

  prob_df$c1 <- prob_df$support == critical_values(odds_ratio = odds_ratio, alternative)[["c1"]]
  prob_df$c2 <- prob_df$support == critical_values(odds_ratio = odds_ratio, alternative)[["c2"]]

  ## Determine critical values (with randomisation) ---------------------------

  gamma <- c("gamma1" = NULL, "gamma2" = NULL)
  alpha_n <- if(alternative == "two.sided") alpha/2 else alpha

  # Define critical values:
  c1 <- min(critical_values(odds_ratio, alternative), na.rm = TRUE)
  c2 <- max(critical_values(odds_ratio, alternative), na.rm = TRUE)

  # Cumulative probabilities of these critical values given OR null:
  c1_prob_cumsum <- prob_df[prob_df$support == c1, "cumsum1"]
  c2_prob_cumsum <- prob_df[prob_df$support == c2, "cumsum2"]

  # Adjustment for when we need to go back/forward 1 in critical values
  # to randomise between:
  c1_adjustment <- 0
  c2_adjustment <- 0

  #
  d1 <- qhyper(alpha_n, obs_c1, obs_c2, obs_r1, lower.tail = TRUE)
  d2 <- qhyper(alpha_n, obs_c1, obs_c2, obs_r1, lower.tail = FALSE)

  #gamma3 <- 1

  if(randomise){

    # WITH RANDOMISATION:

    # Find gamma1, gamma2, such that gamma1*(c1-1) + (1-gamma1)*c1 = alpha/2,
    # and gamma2*(c2+1) + (1-gamma2)*c2 = alpha/2 in two-sided case.

    c1_randomise_between <- if(c1_prob_cumsum <= alpha_n){
      # If our critical value c1 has prob less than alpha/2, randomise between
      # it and the value up from it:
      c(c1, c1+1)
    } else {
      # if not, use only the critical value c1 to
      # assign a factor that takes the prob to alpha_n:
      c1
    }

    # Add these to prob_df:
    prob_df$c1_random <- prob_df$support <= max(c1_randomise_between)

    c2_randomise_between <- if(c2_prob_cumsum <= alpha_n){
      # If our critical value c2 has prob less than alpha/2, randomise between
      # it and the value down from it:
      c(c2-1, c2)
    } else {
      # if not, use only the critical value c2 to
      # assign a factor that takes the prob to alpha_n:
      c2
    }

    # Add these to prob_df:
    prob_df$c2_random <- prob_df$support >= min(c2_randomise_between)

    # If both critical values c1 and c2 need randomised, do it jointly:

    # if(length(c1_randomise_between) == 2 & length(c2_randomise_between) == 2){
    #
    #   gamma3 <- uniroot(
    #     function(x) x*prob_df[prob_df$support %in% c1_randomise_between, "prob"][[2]] +
    #       (1-x)*prob_df[prob_df$support %in% c1_randomise_between, "cumsum1"][[1]] - alpha_n
    #     + x*prob_df[prob_df$support %in% c2_randomise_between, "prob"][[1]] +
    #       (1-x)*prob_df[prob_df$support %in% c2_randomise_between, "cumsum2"][[2]] - alpha_n,
    #     interval = c(0, 1)
    #   )$root
    #
    # }

    # If we need to randomise between two values:

    if(length(c1_randomise_between) == 2){

      gamma1 <- uniroot(
        function(x) x*prob_df[prob_df$support == c1+1, "prob"] +
          (1-x)*prob_df[prob_df$support == c1, "cumsum1"] - alpha_n,
        interval = c(0, 1)
      )$root

      c1_adjustment <- 1

      # Check:

      #if(alpha/2 != gamma1*prob_df[1, "prob"] + (1-gamma1)*prob_df[2, "prob"]){
      #  warning("Gamma1 error")
      #}
    }

    if(length(c2_randomise_between) == 2){

      gamma2 <- uniroot(
        function(x) x*prob_df[prob_df$support == c2-1, "prob"] +
          (1-x)*prob_df[prob_df$support == c2, "prob"] - alpha_n,
        interval = c(0, 1)
      )$root

      c2_adjustment <- -1

      # Check:

      #if(alpha/2 != gamma2*prob_df[1, "prob"] + (1-gamma1)*prob_df[2, "prob"]){
      #  warning("Gamma1 error")
      #}
    }

    # If not:

    if(length(c1_randomise_between) == 1){

      gamma1 <- alpha_n/c1_prob_cumsum

    }

    if(length(c2_randomise_between) == 1){

      gamma2 <- alpha_n/c2_prob_cumsum

    }

    gamma <- c("gamma1" = gamma1, "gamma2" = gamma2)

  }

  ## Calculate our expected value for a given our null odds ratio -------------

  calc_exp_value <- function(odds_ratio){

    # If odds ratios are at extreme ends, return min/max a:

    if(odds_ratio == 0)
      return(min_a)

    if(odds_ratio == Inf)
      return(max_a)

    # Else, return the expected value given our null odds ratio and the
    # non centric hypergeometric distribution:

    sum(support * dn_hypergeo(odds_ratio = odds_ratio))

  }


  ## Gamma calculation via paper's formula ------------------------------------

  gamma_test <- c(NA, NA)
  c1 <- min(critical_values(odds_ratio = odds_ratio, alternative = alternative),
            na.rm = TRUE)
  c2 <- max(critical_values(odds_ratio = odds_ratio, alternative = alternative),
            na.rm = TRUE)

  term1 <- matrix(c(1/prob_df[support == c1, "prob"], 0, 0,
                    1/prob_df[support == c2, "prob"]), nrow = 2, ncol = 2)

  term2 <- matrix(c(c2, -c1, -1, 1), nrow = 2, ncol = 2)

  term3 <- matrix(c(alpha - sum(prob_df[support < c1 | support > c2, "prob"]),
                    alpha*calc_exp_value(odds_ratio = odds_ratio)-sum(
                      prob_df[support < c1 | support > c2, "prob"]*
                        prob_df[support < c1 | support > c2, "support"])),
                  nrow = 2, ncol = 1)

  gamma_test <- c((1/(c2-c1))*term1%*%term2%*%term3)
  names(gamma_test) <- c("gamma1", "gamma2")


  ## Maximum Likelihood Estimator ---------------------------------------------

  # Determine the MLE for the odds ratio by solving E(X) = x, where the
  # expectation is with respect to our (nc) hypergeometric distribution.

  # NOTE the conditional MLE estimator is used rather than the unconditional
  # MLE estimator (sample odds ratio) as in R stats' fisher.exact.text()
  # function.

  # Sample odds ratio: (a*d)/(b*c)

  # When any of these values are 0, the sample odds ratio equals 0 or Inf (it is
  # undefined if both values in a row or in a column are zero). When these
  # outcomes have positive probability, the expected value and variance are not
  # defined. A better well-behaved estimator corrects for this by adding some
  # arbitrarily small number to each cell value, e.g. 0.5.

  # Both the sample odds ratio and the adjusted odds ratio estimators have the
  # same asymptotic normal distribution around the true odds-ratio. Unless the
  # sample size is quite large, however, the distributions of both of these
  # estimators are skewed. The log transformed-odds ratio converges more rapidly
  # to normality. An estimated standard error for log𝜃̂ is

  mle <- function(obs_a){

    # When observed a is already at the minimum or maximum possible value,
    # MLE estimate is a perfect negative/positive relationship:

    if(obs_a == min_a)
      return(0)

    if(obs_a == max_a)
      return(Inf)

    # When it is somewhere in-between, the MLE estimate is given by:

    mu <- calc_exp_value(odds_ratio = 1)

    if(mu > obs_a)
      estimated_odds_ratio <- uniroot(
        f = function(t) calc_exp_value(odds_ratio = t) - obs_a,
        interval = c(0, 1))$root

    else if(mu < obs_a)
      estimated_odds_ratio <- 1 / uniroot(
        f = function(t) calc_exp_value(odds_ratio = 1/t) - obs_a,
        interval = c(.Machine$double.eps, 1))$root

    else
      estimated_odds_ratio <- 1

    return(estimated_odds_ratio)

  }

  ## One-sided p-value --------------------------------------------------------

  # By calculating individual table probabilities using the NCHG distribution
  # we can sum these for our one-tailed p-values:

  one_sided_pvalue <- function(odds_ratio, upper_tail = FALSE){

    # If null hypothesis is OR=1:

    if(odds_ratio == 1){
      if(upper_tail){
        # obs_a1-1 since we want P[X>=x], not P[X>x] (see phyper docs):
        pval <- phyper(obs_a-1, obs_c1, obs_c2, obs_r1,
                                   lower.tail = FALSE)
      } else {
        pval <- phyper(obs_a, obs_c1, obs_c2, obs_r1,
                                   lower.tail = TRUE)
      }
    } else if(odds_ratio == 0){
      # If null hypothesis is OR=0 (perfect negative relationship) use
      # indicator function:
      if(upper_tail){
        pval <- as.numeric(obs_a <= min_a)
      } else {
        pval <- as.numeric(obs_a >= min_a)
      }
    } else if(odds_ratio == Inf){
      # If null hypothesis is OR->Inf (perfect positive relationship) use
      # indicator function:
      if(upper_tail){
        pval <- as.numeric(obs_a <= max_a)
      } else {
        pval <- as.numeric(obs_a >= max_a)
      }
    } else {
      # Or if null hypothesis is OR > 0 and < 1, or > 1:
      if(upper_tail){
        probs <- dn_hypergeo(odds_ratio)[support >= obs_a]
      } else {
        probs <- dn_hypergeo(odds_ratio)[support <= obs_a]
      }
      pval <- sum(probs)
    }

    return(pval)

  }

  ## Two-tailed and final p-value ---------------------------------------------

  # Once we have the one-sided p-value, we can calculate and return
  # our 'final' p-value (i.e. two-tailed if needed).

  calc_pvalue <- function(alternative){

    ## Calculate p-value:

    pvalue <- switch(

      alternative,

      less = one_sided_pvalue(odds_ratio = odds_ratio, upper_tail = FALSE),

      greater = one_sided_pvalue(odds_ratio = odds_ratio, upper_tail = TRUE),

      two.sided = {
        if(odds_ratio == 0){
          # Perfect negative relationship - p-value of 0 or 1, 0 when observed a
          # is exactly the minimum possible a value (not leq or geq):
          as.numeric(obs_a == min_a)
        } else if(odds_ratio == Inf){
          # Perfect positive relationship - p-value of 0 or 1, 0 when observed a
          # is exactly the maximum possible a value (not leq or geq):
          as.numeric(obs_a == max_a)
        } else {
          # For all other cases, it's the sum of all tables as extreme as our
          # observed table - in this case taken to mean tables with probabilities
          # <= our table's probability + some arbitrary small number:
          err <- 1 + 10e-7
          # Sum probabilities:
          probs <- dn_hypergeo(odds_ratio = odds_ratio)

          sum(probs[probs <= probs[obs_a - min_a + 1] * err])
        }
      }
    )

    return(pvalue)

  }

  ## Determine confidence intervals for the estimated odds ratio ----------------

  ### Upper OR ------------------------------------------------------------------

  or_upper <- function(odds_ratio, alpha){

    if(obs_a == max_a)
      return(Inf)

    obs_prob <- one_sided_pvalue(odds_ratio = 1, upper_tail = FALSE)

    if(obs_prob < alpha)
      uniroot(function(t) one_sided_pvalue(odds_ratio = t,
                                           upper_tail = FALSE) - alpha,
              c(0, 1))$root

    else if(obs_prob > alpha)
      1 / uniroot(function(t) one_sided_pvalue(odds_ratio = 1/t,
                                               upper_tail = FALSE) - alpha,
                  c(.Machine$double.eps, 1))$root

    else
      1
  }

  ### Lower OR ------------------------------------------------------------------

  or_lower <- function(odds_ratio, alpha){

    if(obs_a == min_a)
      return(0)

    obs_prob <- one_sided_pvalue(odds_ratio = 1, upper_tail = TRUE)

    if(obs_prob > alpha)
      uniroot(function(t) one_sided_pvalue(odds_ratio = t,
                                           upper_tail = TRUE) - alpha,
              c(0, 1))$root

    else if(obs_prob < alpha)
      1 / uniroot(function(t) one_sided_pvalue(odds_ratio = 1/t,
                                               upper_tail = TRUE) - alpha,
                  c(.Machine$double.eps, 1))$root

    else
      1
  }

  ## Confidence intervals -------------------------------------------------------

  calc_conf_int <- function(odds_ratio, alpha, alternative){

    conf_int <- switch(
      alternative,

      less = c(0, or_upper(odds_ratio, alpha)),

      greater = c(or_lower(odds_ratio, alpha), Inf),

      two.sided = {
        alpha2 <- alpha/2
        c(or_lower(odds_ratio, alpha2),
          or_upper(odds_ratio, alpha2))
      }
    )

    return(conf_int)

  }


  # RESULTS -------------------------------------------------------------------

  ## Calculate p-value --------------------------------------------------------

  RESULT_pval <- calc_pvalue(alternative = alternative)

  ## Critical values ----------------------------------------------------------

  RESULT_crit <- critical_values(odds_ratio = odds_ratio,
                                 alternative = alternative)

  # When randomised, adjust critical values to the value above/below (if
  # there is any):

  RESULT_crit[["c1"]] <- RESULT_crit[["c1"]] + c1_adjustment
  RESULT_crit[["c2"]] <- RESULT_crit[["c2"]] + c2_adjustment

  ## Randomised critical values -----------------------------------------------

  # And store randomisation values:

  RESULT_random <- gamma

  ## Estimate odds ratio ------------------------------------------------------

  RESULT_estimate <- mle(obs_a = obs_a)

  ## Estimate confidence intervals --------------------------------------------

  RESULT_conf_int <- calc_conf_int(odds_ratio = odds_ratio,
                                   alpha = alpha, alternative = alternative)

  ## Calculate test size ------------------------------------------------------

  ## Randomisation ------------------------------------------------------------

  # If t = 0 (R1), then critical values are c1=c2=0, and randomisation
  # parameters are simply alpha/2. If t = 1, critical values c1=0 and c2=1, and
  # randomisation parameters are alpha.

  # Similarly:

  # t = n + m, means gamma1 = gamma2 = alpha/2 and c1 = c2 = m
  # t = n + m - 1, means gamma1 = gamma2 = alpha and c1 = m - 1, c2 = m

  # General solution:

  # (1/c2-c1)()()()

  ## Format results -----------------------------------------------------------

  attr(odds_ratio, "names") <- "odds ratio"
  attr(RESULT_estimate, "names") <- "odds ratio"
  attr(RESULT_conf_int, "conf.level") <- conf.level
  #attr(RESULT_crit, "names")
  #attr(RESULT_random, "names")

  # Put results into list:

  RESULTS <- list(
    p.value = RESULT_pval,
    estimate = RESULT_estimate,
    conf.int = if(conf_int) RESULT_conf_int,
    null.value = odds_ratio,
    alternative = alternative,
    method = METHOD,
    data.name = DATANAME,
    crit.values = RESULT_crit,
    random.probs = if(randomise) RESULT_random,
    support.data = prob_df,
    test.gamma = gamma_test
  )

  attr(RESULTS, "class") <- "htest"

  return(RESULTS)

}


# Test:

test <- fisher_exact_test(data = input_data, alpha = 0.05, odds_ratio = 1,
                  alternative = "two.sided", conf_int = TRUE, randomise = TRUE)

test$crit.values
test$random.probs
test$test.gamma

fisher_exact_test(data = input_data, alpha = 0.05, odds_ratio = 1,
                  alternative = "less", conf_int = TRUE, randomise = FALSE)

fisher_exact_test(data = input_data, alpha = 0.05, odds_ratio = 1,
                  alternative = "greater", conf_int = TRUE, randomise = FALSE)


test$support.data


library(dplyr)
library(tidyr)

testing_frame <- tibble(
  "C1" = 6, # U(m = 6, pi1 = (0,1))
  "C2" = 4, # V(n = 4, pi2 = (0,1))
  "R1" = 0:10,
  "R2" = 10:0,
  #"R2" = C1+C2-R1,
  "N" = C1+C2) %>%
  rowwise() %>%
  mutate(
    min = max(R1-C2, 0),
    max = min(C1, R1)
    ) %>%
  ungroup()

# Expand possible values for a:

testing_frame <- testing_frame %>%
  pivot_longer(min:max, names_to = "flag", values_to = "test_a") %>%
  group_by_at(vars(C1:N)) %>%
  complete(test_a = min(test_a):max(test_a)) %>%
  ungroup()

# Create contingency table permutations:

testing_frame <- testing_frame %>%
  mutate(
    test_b = R1-test_a,
    test_c = C1-test_a,
    test_d = C2-test_b
  )

# FET:

testing_frame <- testing_frame %>%
  rowwise() %>%
  mutate(
    FET = fisher_exact_test(data = matrix(c(test_a, test_c, test_b, test_d), nrow = 2, ncol = 2), odds_ratio = 1,
                            alternative = "two.sided", alpha = 0.05, randomise = TRUE)$p.value,
    crit1 = fisher_exact_test(data = matrix(c(test_a, test_c, test_b, test_d), nrow = 2, ncol = 2), odds_ratio = 1,
                              alternative = "two.sided", alpha = 0.05, randomise = TRUE)$crit.values[["c1"]],
    crit2 = fisher_exact_test(data = matrix(c(test_a, test_c, test_b, test_d), nrow = 2, ncol = 2), odds_ratio = 1,
                              alternative = "two.sided", alpha = 0.05, randomise = TRUE)$crit.values[["c2"]],
    gamma1 = fisher_exact_test(data = matrix(c(test_a, test_c, test_b, test_d), nrow = 2, ncol = 2), odds_ratio = 1,
                               alternative = "two.sided", alpha = 0.05, randomise = TRUE)$random.probs[["gamma1"]],
    gamma2 = fisher_exact_test(data = matrix(c(test_a, test_c, test_b, test_d), nrow = 2, ncol = 2), odds_ratio = 1,
                               alternative = "two.sided", alpha = 0.05, randomise = TRUE)$random.probs[["gamma2"]],
    gamma1_test = fisher_exact_test(data = matrix(c(test_a, test_c, test_b, test_d), nrow = 2, ncol = 2), odds_ratio = 1,
                               alternative = "two.sided", alpha = 0.05, randomise = TRUE)$test.gamma[["gamma1"]],
    gamma2_test = fisher_exact_test(data = matrix(c(test_a, test_c, test_b, test_d), nrow = 2, ncol = 2), odds_ratio = 1,
                                    alternative = "two.sided", alpha = 0.05, randomise = TRUE)$test.gamma[["gamma2"]]
  ) %>%
  ungroup()

# Reduce number of rows again by removing contingency table iterations for
# every row/column sum:

testing_frame %>%
  select(C1:N, crit1, crit2, gamma1, gamma2, gamma1_test, gamma2_test) %>%
  distinct()

troubleshoot <- testing_frame %>%
  filter(R1 == 6) %>%
  filter(test_a == last(test_a))



library(ggplot2)

plotData <- testing_frame %>%
  select(C1, test_a:test_d, FET, -flag) %>%
  distinct() %>%
  mutate(
    size = test_a/C1
  )

ggplot(plotData, aes(x = test_a/C1, y = FET)) +
  geom_line()
