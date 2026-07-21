# Draw a categorical index from probabilities

Given a uniform random draw and a vector of probabilities, returns the
1-based integer index of the selected category.

## Usage

``` r
draw_categorical(luck, probs)
```

## Arguments

- luck:

  Numeric between 0 and 1 (inclusive). A uniform random draw.

- probs:

  Numeric vector of probabilities. Will be normalized to sum to 1
  internally. All values must be \>= 0.

## Value

Integer: 1-based index of the selected category.
