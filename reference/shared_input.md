# Shared input object

Constructor for a lightweight "shared or immutable" value holder.

## Usage

``` r
shared_input(expr, constrained = NULL)
```

## Arguments

- expr:

  A value or expression to initialize the shared input with. The
  expression is evaluated immediately.

- constrained:

  Logical. If `TRUE`, creates a shared environment-backed object. If
  `FALSE`, creates an immutable copy-on-modify object. If `NULL`
  (default), the function looks up `constrained` in the calling
  environment; only an explicit `TRUE` enables shared mode.

## Value

An object of class `shared_input_val` (immutable mode) or
`shared_input_env` (shared mode), both inheriting from class
`"shared_input"`. Each instance exposes the following user methods:

- \$value():

  Returns the current stored value.

- \$modify(new_v):

  In immutable mode: returns a new independent wrapper with updated
  value. In shared mode: updates the shared value by reference and
  returns a new wrapper pointing to the same shared state.

- \$clone():

  Returns a deep copy (independent wrapper and independent internal
  state). Subsequent modifications on clones do not affect the original
  object or its aliases.

- \$reset():

  Returns a new wrapper whose value is restored to the original
  initialization value. In both modes this creates an independent fresh
  state.

- \$fork(n):

  Creates `n` independent deep clones as a list. Useful for generating
  multiple isolated copies quickly.

## Details

`shared_input()` produces a simple object that wraps a value with
controlled mutability semantics. It can operate in two distinct modes:

- **Immutable (non-shared)**: every modification produces a fresh,
  independent copy of the object (safe for parallel or functional code).

- **Shared (constrained)**: the object’s value is stored in a common
  environment shared across all aliases (by-reference semantics). This
  allows coordinated updates across multiple handles.

The mode is determined either by the explicit argument `constrained`, or
by inheriting the value of a `constrained` variable in the parent frame.

- In **immutable mode**, each wrapper stores its value in closures
  (`make_val()`) and is fully copy-on-modify. No references are shared.

- In **shared mode**, all wrappers produced by `$modify()` or direct
  aliasing point to the same underlying environment (`state`). This
  means updating one updates all aliases until a `$clone()` or
  `$reset()` breaks the link.

The underlying `state` environments are internal. Users should rely only
on the public methods above.

Note: if the stored value itself is a reference type (e.g., environment,
external pointer, R6 object), those internal references remain shared
regardless of mode, following normal R semantics.

## Examples

``` r
# --- Immutable (default) mode ---
a <- shared_input(5)
a$value()                 # 5
#> [1] 5
a2 <- a$modify(a$value() + 7)
a$value()                 # 5
#> [1] 5
a2$value()                # 12
#> [1] 12

# Cloning and resetting
a3 <- a2$clone()
a4 <- a2$reset()
a3$value(); a4$value()    # 12, 5
#> [1] 12
#> [1] 5

# Forking
forks <- a$fork(3)
vapply(forks, function(x) x$value(), numeric(1))
#> [1] 5 5 5

# --- Shared (constrained) mode ---
constrained <- TRUE
b1 <- shared_input(10)
b2 <- b1        # alias (same state)
b1$modify(11)
#> <environment: 0x55ced9b51138>
#> attr(,"class")
#> [1] "shared_input_env" "shared_input"    
b1$value(); b2$value()  # both 11
#> [1] 11
#> [1] 11

b3 <- b1$clone()
b1$modify(99)
#> <environment: 0x55ced94fa030>
#> attr(,"class")
#> [1] "shared_input_env" "shared_input"    
b1$value(); b3$value()  # 99, 11
#> [1] 99
#> [1] 11

# Reset breaks sharing
b4 <- b1$reset()
b4$value()              # 10
#> [1] 10
```
