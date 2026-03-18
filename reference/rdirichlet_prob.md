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
#>           [,1]       [,2]      [,3]
#> [1,] 0.7614474 0.09334006 0.1452125
#> [2,] 0.2574358 0.56109816 0.1814661
#> [3,] 0.0000000 0.22882178 0.7711782

rdirichlet_prob(n=2,alpha= matrix(c(0.7,0.3,0,0.1,0.7,0.2,0.1,0.2,0.7),3,3),
se=matrix(c(0.7,0.3,0,0.1,0.7,0.2,0.1,0.2,0.7)/10,3,3))
#> [[1]]
#>           [,1]       [,2]       [,3]
#> [1,] 0.8033745 0.09897899 0.09764656
#> [2,] 0.2676214 0.55231391 0.18006468
#> [3,] 0.0000000 0.23225018 0.76774982
#> 
#> [[2]]
#>           [,1]      [,2]      [,3]
#> [1,] 0.7831308 0.1044302 0.1124389
#> [2,] 0.2055862 0.6073687 0.1870451
#> [3,] 0.0000000 0.1772613 0.8227387
#> 
```
