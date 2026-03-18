# Draw from a dirichlet distribution based on number of counts in transition. Adapted from brms::rdirichlet

Draw from a dirichlet distribution based on number of counts in
transition. Adapted from brms::rdirichlet

## Usage

``` r
rdirichlet(n = 1, alpha, seed = NULL)
```

## Arguments

- n:

  Number of draws (must be \>= 1). If n\>1, it will return a list of
  matrices.

- alpha:

  A matrix of alphas (transition counts)

- seed:

  An integer which will be used to set the seed for this draw.

## Value

A transition matrix. If n\>1, it will return a list of matrices.

## Examples

``` r
rdirichlet(n=1,alpha= matrix(c(1251, 0, 350, 731),2,2))
#>           [,1]      [,2]
#> [1,] 0.7901734 0.2098266
#> [2,] 0.0000000 1.0000000
rdirichlet(n=2,alpha= matrix(c(1251, 0, 350, 731),2,2))
#> [[1]]
#>           [,1]      [,2]
#> [1,] 0.7827959 0.2172041
#> [2,] 0.0000000 1.0000000
#> 
#> [[2]]
#>           [,1]      [,2]
#> [1,] 0.7739751 0.2260249
#> [2,] 0.0000000 1.0000000
#> 
```
