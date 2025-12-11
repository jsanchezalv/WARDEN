# Draw from a dirichlet distribution based on mean transition probabilities and standard errors

Draw from a dirichlet distribution based on mean transition
probabilities and standard errors

## Usage

``` r
rdirichlet_prob(n = 1, alpha, se, seed = NULL)
```

## Arguments

- n:

  Number of draws (must be \>= 1). If n\>1, it will return a list of
  matrices.

- alpha:

  A matrix of transition probabilities

- se:

  A matrix of standard errors

- seed:

  An integer which will be used to set the seed for this draw.

## Value

A transition matrix. If n\>1, it will return a list of matrices.

## Examples

``` r
rdirichlet_prob(n=1,alpha= matrix(c(0.7,0.3,0,0.1,0.7,0.2,0.1,0.2,0.7),3,3),
se=matrix(c(0.7,0.3,0,0.1,0.7,0.2,0.1,0.2,0.7)/10,3,3))
#>           [,1]       [,2]       [,3]
#> [1,] 0.8086972 0.09241709 0.09888572
#> [2,] 0.1995130 0.64954789 0.15093915
#> [3,] 0.0000000 0.21477410 0.78522590

rdirichlet_prob(n=2,alpha= matrix(c(0.7,0.3,0,0.1,0.7,0.2,0.1,0.2,0.7),3,3),
se=matrix(c(0.7,0.3,0,0.1,0.7,0.2,0.1,0.2,0.7)/10,3,3))
#> [[1]]
#>           [,1]      [,2]      [,3]
#> [1,] 0.7679586 0.1169805 0.1150609
#> [2,] 0.2439433 0.5863156 0.1697411
#> [3,] 0.0000000 0.2274894 0.7725106
#> 
#> [[2]]
#>           [,1]      [,2]       [,3]
#> [1,] 0.7925735 0.1160001 0.09142641
#> [2,] 0.2927343 0.5605694 0.14669634
#> [3,] 0.0000000 0.1854090 0.81459104
#> 
```
