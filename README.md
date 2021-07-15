
<!-- README.md is generated from README.Rmd. Please edit that file -->

# penfa

<!-- badges: start -->

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Last-changedate](https://img.shields.io/badge/last%20change-2021--07--15-brightgreen.svg)](https://github.com/egeminiani/penfa/commits/main)
[![Website](https://img.shields.io/badge/website-penfa-orange.svg?colorB=E91E63)](https://egeminiani.github.io/penfa/)
[![Licence](https://img.shields.io/badge/licence-GPL--3-orange.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![R-CMD-check](https://github.com/egeminiani/penfa/workflows/R-CMD-check/badge.svg)](https://github.com/egeminiani/penfa/actions)
<!-- badges: end -->

### Overview

An R package for estimating single- and multiple-group penalized factor
models via a trust-region algorithm with integrated automatic multiple
tuning parameter selection (Geminiani et al., 2021). Supported penalties
include lasso, adaptive lasso, scad, mcp, and ridge.

### Installation

You can install the released version of penfa from CRAN with:

``` r
install.packages("penfa")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("egeminiani/penfa")
```

### Example

This is a basic example showing how to fit a *PENalized Factor Analysis*
model with the alasso penalty and the automatic tuning procedure. A
shrinkage penalty is applied to the whole factor loading matrix.

Let’s load the data (see `?ccdata` for details).

``` r
library(penfa)
data(ccdata)
```

<font size="4">**Step 1**</font> : specify the model syntax

``` r
syntax = 'help  =~   h1 + h2 + h3 + h4 + h5 + h6 + h7 + 0*v1 + v2 + v3 + v4 + v5
          voice =~ 0*h1 + h2 + h3 + h4 + h5 + h6 + h7 +   v1 + v2 + v3 + v4 + v5'
```

<font size="4">**Step 2**</font>: fit the model

``` r
alasso_fit <- penfa(model  = syntax,
                    data   = ccdata,
                    std.lv = TRUE,
                    pen.shrink = "alasso")
#> Computing weights for alasso (ML estimates)... done.
#> 
#> Automatic procedure: 
#> Iteration  1 : 0.00298271 
#> Iteration  2 : 0.00452604 
#> 
#> Largest absolute gradient value: 12.76355181
#> Fisher information matrix is positive definite
#> Eigenvalue range: [180.2917, 9189645]
#> Trust region iterations: 15 
#> Factor solution: admissible 
#> Effective degrees of freedom: 27.12936
```

``` r
alasso_fit
#> penfa 0.1.1 reached convergence
#> 
#>   Number of observations                                    767
#>                                                                
#>   Estimator                                                PMLE
#>   Optimization method                              trust-region
#>   Information                                            fisher
#>   Strategy                                                 auto
#>   Number of iterations (total)                               58
#>   Number of two-steps (automatic)                             2
#>   Effective degrees of freedom                           27.129
#>                                                                
#>   Penalty function:                                            
#>     Sparsity                                             alasso
#>                                                                
#> 
```

<font size="4">**Step 3**</font>: inspect the results

``` r
summary(alasso_fit)
#> penfa 0.1.1 reached convergence
#> 
#>   Number of observations                                    767
#>   Number of groups                                            1
#>   Number of observed variables                               12
#>   Number of latent factors                                    2
#>                                                                
#>   Estimator                                                PMLE
#>   Optimization method                              trust-region
#>   Information                                            fisher
#>   Strategy                                                 auto
#>   Number of iterations (total)                               58
#>   Number of two-steps (automatic)                             2
#>   Influence factor                                            4
#>   Number of parameters:                                        
#>     Free                                                     13
#>     Penalized                                                22
#>   Effective degrees of freedom                           27.129
#>   GIC                                                 17222.980
#>   GBIC                                                17348.928
#>                                                                
#>   Penalty function:                                            
#>     Sparsity                                             alasso
#>                                                                
#>   Additional tuning parameter                                  
#>     alasso                                                    1
#>                                                                
#>   Optimal tuning parameter:                                    
#>     Sparsity                                                   
#>      - Factor loadings                                    0.005
#>                                                                
#> 
#> Parameter Estimates:
#> 
#> Latent Variables:
#>                     Type    Estimate  Std.Err     2.5%    97.5%
#>   help =~                                                      
#>     h1               pen       0.766    0.030    0.707    0.825
#>     h2               pen       0.858    0.028    0.803    0.913
#>     h3               pen       0.775    0.030    0.717    0.834
#>     h4               pen       0.921    0.038    0.847    0.995
#>     h5               pen       0.810    0.040    0.732    0.887
#>     h6               pen       0.782    0.044    0.696    0.868
#>     h7               pen       0.523    0.050    0.426    0.620
#>     v1             fixed       0.000             0.000    0.000
#>     v2               pen       0.000                           
#>     v3               pen       0.000                           
#>     v4               pen       0.000                           
#>     v5               pen      -0.000                           
#>   voice =~                                                     
#>     h1             fixed       0.000             0.000    0.000
#>     h2               pen      -0.000                           
#>     h3               pen       0.000                           
#>     h4               pen      -0.041                           
#>     h5               pen       0.053    0.031   -0.008    0.114
#>     h6               pen       0.104    0.038    0.029    0.180
#>     h7               pen       0.341    0.049    0.246    0.437
#>     v1               pen       0.851    0.028    0.795    0.906
#>     v2               pen       0.871    0.028    0.817    0.926
#>     v3               pen       0.842    0.029    0.786    0.898
#>     v4               pen       0.843    0.029    0.787    0.899
#>     v5               pen       0.805    0.029    0.747    0.862
#> 
#> Covariances:
#>                     Type    Estimate  Std.Err     2.5%    97.5%
#>   help ~~                                                      
#>     voice           free       0.877    0.011    0.855    0.900
#> 
#> Variances:
#>                     Type    Estimate  Std.Err     2.5%    97.5%
#>    .h1              free       0.388    0.021    0.346    0.429
#>    .h2              free       0.233    0.014    0.205    0.261
#>    .h3              free       0.372    0.021    0.332    0.413
#>    .h4              free       0.184    0.012    0.160    0.209
#>    .h5              free       0.235    0.014    0.207    0.263
#>    .h6              free       0.201    0.012    0.177    0.225
#>    .h7              free       0.264    0.015    0.235    0.293
#>    .v1              free       0.245    0.015    0.216    0.275
#>    .v2              free       0.208    0.014    0.182    0.235
#>    .v3              free       0.261    0.016    0.230    0.292
#>    .v4              free       0.259    0.016    0.228    0.290
#>    .v5              free       0.324    0.019    0.287    0.361
#>     help           fixed       1.000             1.000    1.000
#>     voice          fixed       1.000             1.000    1.000
```

### Vignettes and Tutorials

-   See `vignette("automatic-tuning-selection")` for the estimation of a
    penalized factor model with lasso and alasso penalties. The tuning
    parameter producing the optimal amount of sparsity in the factor
    loading matrix is found through the automatic tuning procedure.

-   See `vignette("grid-search-tuning-selection")` for the estimation of
    a penalized factor model with scad and mcp penalties. A grid search
    is conducted, and the optimal tuning parameter is the one generating
    the penalized model with the lowest GBIC (Generalized Bayesian
    Information Criterion).

-   See
    [“multiple-group-analysis”](https://egeminiani.github.io/penfa/articles/articles/multiple-group-analysis.html)
    for the estimation of a multiple-group penalized factor model with
    the alasso penalty. This model encourages sparsity in the loading
    matrices and cross-group invariance of loadings and intercepts. The
    automatic multiple tuning parameter procedure is employed for
    finding the optimal tuning parameter vector.

-   See
    [“plotting-penalty-matrix”](https://egeminiani.github.io/penfa/articles/articles/plotting-penalty-matrix.html)
    for details on how to produce interactive plots of the penalty
    matrices.

### Literature

-   Geminiani, E., Marra, G., & Moustaki, I. (2021). “Single- and
    Multiple-Group Penalized Factor Analysis: A Trust-Region Algorithm
    Approach with Integrated Automatic Multiple Tuning Parameter
    Selection.” Psychometrika, 86(1), 65-95.
    <https://doi.org/10.1007/s11336-021-09751-8>

-   Geminiani, E. (2020). “A Penalized Likelihood-Based Framework for
    Single and Multiple-Group Factor Analysis Models.” PhD thesis,
    University of Bologna. <http://amsdottorato.unibo.it/9355/>.

### How to cite

    #> 
    #> To cite penfa in publications use:
    #> 
    #>   Geminiani, E., Marra, G., & Moustaki, I. (2021). Single- and
    #>   Multiple-Group Penalized Factor Analysis: A Trust-Region Algorithm
    #>   Approach with Integrated Automatic Multiple Tuning Parameter
    #>   Selection.  Psychometrika, 86(1), 65-95.
    #>   https://doi.org/10.1007/s11336-021-09751-8
    #> 
    #> A BibTeX entry for LaTeX users is
    #> 
    #>   @Article{,
    #>     title = {Single- and Multiple-Group Penalized Factor Analysis: A Trust-Region Algorithm Approach with Integrated Automatic Multiple Tuning Parameter Selection},
    #>     author = {Geminiani Elena and Marra Giampiero and Moustaki Irini},
    #>     journal = {Psychometrika},
    #>     year = {2021},
    #>     volume = {86},
    #>     number = {1},
    #>     pages = {65-95},
    #>     url = {https://doi.org/10.1007/s11336-021-09751-8},
    #>   }
    #> 
    #>   Elena Geminiani, Giampiero Marra and Irini Moustaki (2021). penfa:
    #>   Single- And Multiple-Group Penalized Factor Analysis. R package
    #>   version 0.1.1.
    #> 
    #> A BibTeX entry for LaTeX users is
    #> 
    #>   @Manual{,
    #>     title = {penfa: Single- And Multiple-Group Penalized Factor Analysis},
    #>     author = {Elena Geminiani and Giampiero Marra and Irini Moustaki},
    #>     year = {2021},
    #>     note = {R package version 0.1.1},
    #>   }
