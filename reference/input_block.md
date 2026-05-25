# Build an input block for automatic parameter selection

Creates an unevaluated [`{}`](https://rdrr.io/r/base/Paren.html)
expression that calls
[`pick_val_v()`](https://jsanchezalv.github.io/WARDEN/reference/pick_val_v.md)
with the correct arguments for base case, PSA, DSA, and scenario
analyses. The expression is meant to be passed directly to
[`run_sim()`](https://jsanchezalv.github.io/WARDEN/reference/run_sim.md)
or
[`run_sim_parallel()`](https://jsanchezalv.github.io/WARDEN/reference/run_sim_parallel.md)
as a `*_inputs` argument.

## Usage

``` r
input_block(
  .data = NULL,
  base,
  psa,
  sens,
  names_out,
  psa_indicators = NULL,
  sens_indicators = NULL,
  indicator_sens_binary = FALSE,
  dsa_names = NULL,
  distributions = NULL,
  covariances = NULL
)
```

## Arguments

- .data:

  Optional existing [`{}`](https://rdrr.io/r/base/Paren.html) block to
  prepend (for pipe chaining).

- base:

  A list of base case values, one entry per parameter.

- psa:

  An unevaluated expression producing PSA draws (e.g. `pick_psa(...)`),
  evaluated at runtime per simulation.

- sens:

  The sensitivity data object (e.g. `l_inputs`). At runtime the engine
  indexes it as `sens[[sens_name_used]]` to pick the active column.

- names_out:

  Character vector or list of strings of output parameter names.

- psa_indicators:

  List of 0/1 indicators controlling which parameters draw from PSA.
  `NULL` means all parameters draw from PSA.

- sens_indicators:

  List of indicators, one entry per parameter (may be a scalar or a
  vector for vector-valued parameters). In binary mode: 0 = inactive, 1
  = active. In grouped mode: 0 = inactive, same integer = varied
  together. `NULL` means all parameters are active. Only in grouped mode
  can subparameters (e.g., elements of a vector parameter) be varied at
  different steps of a DSA.

- indicator_sens_binary:

  Logical. `TRUE` = binary (0/1) mode; each active parameter is varied
  independently via
  [`create_indicators()`](https://jsanchezalv.github.io/WARDEN/reference/create_indicators.md).
  `FALSE` (default) = grouped integer mode when `sens_indicators` and
  `distributions` are supplied; otherwise falls back to binary mode.

- dsa_names:

  Character vector of `sensitivity_names` values that are DSA directions
  (e.g. `c("DSA_min","DSA_max")`). Those names iterate through active
  parameters/groups one at a time. All other names passed to
  `sensitivty_names` in
  [`run_sim()`](https://jsanchezalv.github.io/WARDEN/reference/run_sim.md)
  are treated as scenarios. `NULL` (the default) means all sensitivity
  names are scenarios.

- distributions:

  List of distribution names (e.g. `"rnorm"`, `"mvrnorm"`,
  `"rdirichlet"`), required in grouped mode to determine which
  parameters are vectors and require conditional distributions.

- covariances:

  List of covariance matrices or scalars, required for `"mvrnorm"` /
  `"rdirichlet"` distributions.

## Value

A [`{}`](https://rdrr.io/r/base/Paren.html) language object (same type
as
[`add_item()`](https://jsanchezalv.github.io/WARDEN/reference/add_item.md))
with a `warden_block_meta` attribute used by the engine to build the
iteration schedule.

## Details

Supports two modes:

- **Binary** (`indicator_sens_binary = TRUE`): `sens_indicators` are 0/1
  per parameter. Each non-zero parameter is varied independently via
  [`create_indicators()`](https://jsanchezalv.github.io/WARDEN/reference/create_indicators.md).

- **Grouped** (`indicator_sens_binary = FALSE`, the default when
  `sens_indicators` and `distributions` are supplied): `sens_indicators`
  are integer group labels. Parameters sharing the same integer are
  varied together; 0 means the parameter is never varied.

### DSA vs scenario analyses

Supply `dsa_names` with the subset of `sensitivity_names` that
correspond to DSA directions (e.g. `c("DSA_min","DSA_max")`). Each DSA
name runs one iteration per active parameter/group; other names
(scenarios) apply all active parameters simultaneously in a single
iteration.

When `dsa_names = NULL` (the default) every sensitivity name is treated
as a scenario — all active parameters take their scenario value at once,
one iteration per `sensitivity_names` entry. This is the correct default
when there is no DSA.

### `sens_indicators` and zeros

Any parameter whose `sens_indicators` entry is `0` (or all-zero for
vector-valued parameters) is **never varied**: it always uses the base
or PSA value regardless of whether we are running a DSA or a scenario.
This applies to both modes.

## See also

[`pick_val_v()`](https://jsanchezalv.github.io/WARDEN/reference/pick_val_v.md),
[`add_item()`](https://jsanchezalv.github.io/WARDEN/reference/add_item.md),
[`run_sim()`](https://jsanchezalv.github.io/WARDEN/reference/run_sim.md)

## Examples

``` r
l_inputs <- list(
  parameter_name = list("util.sick", "util.sicker"),
  base_value     = list(0.8, 0.5),
  PSA_dist       = list("rnorm", "rbeta_mse"),
  a              = list(0.8, 0.5),
  b              = list(0.04, 0.025),
  n              = list(1, 1),
  DSA_min        = list(0.6, 0.3),
  DSA_max        = list(0.9, 0.7),
  psa_indicators = list(1, 1)
)

blk <- input_block(
  base                  = l_inputs[["base_value"]],
  psa                   = pick_psa(l_inputs[["PSA_dist"]], l_inputs[["n"]],
                                   l_inputs[["a"]], l_inputs[["b"]]),
  sens                  = l_inputs,
  names_out             = l_inputs[["parameter_name"]],
  psa_indicators        = l_inputs[["psa_indicators"]],
  indicator_sens_binary = TRUE,
  dsa_names             = c("DSA_min", "DSA_max")
)
stopifnot(is.call(blk))
stopifnot(!is.null(attr(blk, "warden_block_meta")))
```
