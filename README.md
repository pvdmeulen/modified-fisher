
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Non-Conservative Size-α Modified Fisher’s Exact Test

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R-CMD-check](https://github.com/pvdmeulen/modifiedfisher/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pvdmeulen/modifiedfisher/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Installation

You can install the latest version of this `modifiedfisher` package by
using `devtools`:

``` r
devtools::install_github("pvdmeulen/modifiedfisher")
```

## Examples

The modified Fisher Exact Test (see more information below) can be
called using the `modified_fisher_exact_test()` function:

``` r
library(modifiedfisher)
modified_fisher_exact_test(u = 5, m = 12, v = 7, n = 11, odds_ratio = 1, alpha = 0.05)
#> 
#>  Non-Conservative Size-α Modified Fisher's Exact Test
#> 
#> data:  u = 5, v = 7, m = 12, n = 11
#> p-value = 0.3208
#> alternative hypothesis: true odds ratio is not equal to 1
#> 95 percent confidence interval:
#>  0.072034 2.209549
#> sample estimates:
#> odds ratio 
#>  0.4081633
```

The above example uses $m = 12$ and $n = 11$, echoing Figure 1(a) and
Table 2 in the paper linked below. As you can see, the p-value,
estimated odds ratio, and confidence intervals match the paper’s. A
similar size plot can be constructed by setting the `local_size_data`
argument to `TRUE`, and using the resulting local size dataframe to
construct the plot:

``` r
test_result <- modified_fisher_exact_test(u = 5, m = 12, v = 7, n = 11, odds_ratio = 1,
    alpha = 0.05, local_size_data = TRUE)
```

The data is now stored in the `test_result$local.size.data` object:

|  pi1 |      size | method |
|-----:|----------:|:-------|
| 0.00 | 0.0000000 | zoom   |
| 0.01 | 0.0006922 | zoom   |
| 0.02 | 0.0092791 | zoom   |
| 0.03 | 0.0393075 | zoom   |
| 0.04 | 0.1038314 | zoom   |
| 0.05 | 0.2116501 | zoom   |

Plotting this leads to a similar plot as Figure 1(a):

<picture>
<source media="(prefers-color-scheme: dark)" srcset="man/figures/README-dark_plot_data-1.svg">
<source media="(prefers-color-scheme: light)" srcset="man/figures/README-light_plot_data-1.svg">
<img alt="Shows a graph with test size on the y-axis and nuisance parameter pi_1 between zero and one on the x-axis, with a line showing the relationship between the two. The line is shaped like two humps, with the maximum being rouhgly 4%.">
</picture>

## The Modified Fisher Exact Test

### [Read the article here](https://www.researchgate.net/publication/351111885_Consistent_Confidence_Limits_P_Values_and_Power_of_the_Non-Conservative_Size_-a_Modified_Fisher_Exact_Test)

- Add brief rationale and explanation of theory

## Next steps

- Implement checks and unit tests
- Document functions in a consistent way with article / SAS macro
- Add randomised Fisher Exact Test option
