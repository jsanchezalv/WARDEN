# Calculate conditional dirichlet values

Calculate conditional dirichlet values

## Usage

``` r
cond_dirichlet(alpha, i, xi, full_output = FALSE)
```

## Arguments

- alpha:

  mean vector

- i:

  index of the known parameter (1-based index)

- xi:

  known value of the i-th parameter (should be \>0)

- full_output:

  boolean indicating whether to return the full list of parameters

## Value

List of length 2, one with new mu and other with covariance parameters

## Details

Function to compute conditional dirichlet values

## Examples

``` r
alpha <- c(2, 3, 4)
i <- 2  # Index of the known parameter
xi <- 0.5  # Known value of the second parameter

# Compute the conditional alpha parameters with full output
cond_dirichlet(alpha, i, xi, full_output = TRUE)
#> [1] 0.1666667 0.5000000 0.3333333
```
