# Calculate conditional multivariate normal values

Calculate conditional multivariate normal values

## Usage

``` r
cond_mvn(mu, Sigma, i, xi, full_output = FALSE)
```

## Arguments

- mu:

  mean vector

- Sigma:

  covariance matrix

- i:

  index of the known parameter (1-based index)

- xi:

  known value of the i-th parameter

- full_output:

  boolean indicating whether to return the full list of parameters

## Value

List of length 2, one with new mu and other with covariance parameters

## Details

Function to compute conditional multivariate normal values

## Examples

``` r
mu <- c(1, 2, 3)
Sigma <- matrix(c(0.2, 0.05, 0.1, 
                  0.05, 0.3, 0.05, 
                  0.1, 0.05, 0.4), nrow = 3)

i <- 1:2  # Index of the known parameter
xi <- c(1.2,2.3)  # Known value of the first parameter

cond_mvn(mu, Sigma, i, xi,full_output = TRUE)
#> $mean
#> [1] 1.200000 2.300000 3.121739
#> 
#> $covariance
#>      [,1] [,2]      [,3]
#> [1,]    0    0 0.0000000
#> [2,]    0    0 0.0000000
#> [3,]    0    0 0.3478261
#> 
```
