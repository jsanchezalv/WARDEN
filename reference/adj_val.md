# Adjusted Value Calculation

This function calculates an adjusted value over a time interval with
optional discounting. This is useful for instances when adding cycles
may not be desirable, so one can perform "cycle-like" calculations
without needing cycles, offering performance speeds. See the vignette on
avoiding cycles for an example in a model.

## Usage

``` r
adj_val(
  curtime,
  nexttime,
  by,
  expression,
  discount = NULL,
  vectorized_f = FALSE
)
```

## Arguments

- curtime:

  Numeric. The current time point.

- nexttime:

  Numeric. The next time point. Must be greater than or equal to
  `curtime`.

- by:

  Numeric. The step size for evaluation within the interval.

- expression:

  An expression evaluated at each step. Use `.time` as the variable
  within the expression.

- discount:

  Numeric or NULL. The discount rate to apply, or NULL for no
  discounting.

- vectorized_f:

  boolean, FALSE by default. If TRUE, evaluates the expression once
  using `.time` as a vector. If FALSE, it repeatedly evaluates the
  expression with time as a single value (slower).

## Value

Numeric. The calculated adjusted value.

## Details

The user can use the `.time` variable to select the corresponding time
of the sequence being evaluated. For example, in
`curtime = 0, nexttime = 4, by = 1`, `.time` would correspond to
`0, 1, 2, 3`. If using `nexttime = 4.2`, `0, 1, 2, 3, 4`

## Examples

``` r
# Define a function or vector to evaluate
bs_age <- 1
vec <- 1:8/10

# Calculate adjusted value without discounting
adj_val(0, 4, by = 1, expression = vec[floor(.time + bs_age)])
#> [1] 0.25
adj_val(0, 4, by = 1, expression = .time * 1.1)
#> [1] 1.65
#same result since .time * 1.1 can be vectorized w.r.t time
adj_val(0, 4, by = 1, expression = .time * 1.1, vectorized_f = TRUE)
#> [1] 1.65

# Calculate adjusted value with discounting
adj_val(0, 4, by = 1, expression = vec[floor(.time + bs_age)], discount = 0.03)
#> [1] 0.2463061
```
