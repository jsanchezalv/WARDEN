## Progress

### Conversation summary (2026-05-20): `pick_val_v` simplification

**Goal**: Simplify the `pick_val_v()` workflow so users don't have to manually manage iterators, indicators, offsets, or `n_sensitivity`.

**Discussion outcome**: We agreed on an `input_block()` function that:
- Takes explicit arguments (user still specifies base, psa, sens columns, etc.)
- Returns a `{ }` expression (same type as `add_item()`) so it's composable via pipe
- Carries metadata as an attribute (inspected ONCE by `run_sim` at setup, not in hot loop)
- The engine auto-computes `n_sensitivity`, offsets, and generates the indicator/iterator logic

**Bug found**: In grouped mode (`indicator_sens_binary = FALSE`), `indicator_psa` is completely ignored. When `psa_ind = TRUE`, all parameters start from PSA regardless of `indicator_psa`. Fix: respect `indicator_psa` per-element before applying DSA overrides while keeping computation fast.

---

### Implementation plan

#### Phase 1: Fix `pick_val_v` grouped mode PSA bug

**File**: `R/input_f.R`, lines 282-315

**Problem**: When `indicator_sens_binary = FALSE` and `sens_ind = TRUE`, the code does:
```r
if(psa_ind) { temp_data <- psa } else { temp_data <- base }
output <- temp_data
```
This ignores `indicator_psa`, so ALL parameters start from PSA when `psa_ind = TRUE`.

**Fix**: Build `temp_data` element-by-element respecting `indicator_psa`:
- If `indicator_psa` is NULL, default to all-1s (current behavior, all PSA)
- Otherwise, for each parameter `i`: use `psa[[i]]` if `indicator_psa[[i]]` has any 1, else `base[[i]]`
- For vector parameters with partial PSA (e.g., `indicator_psa = c(1,0)`): element-wise selection between psa and base

**Performance**: Vectorize where possible. Avoid per-element R loop for the common case (all scalar indicators). Use `vapply` or direct indexing.

**Tests**: Add tests in `tests/testthat/test-input_f.R` covering:
- Grouped mode with `indicator_psa = list(0,0,0,1)` — only param 4 from PSA
- Grouped mode with vector `indicator_psa = list(1, c(1,0))` — partial PSA on vector param
- Grouped mode with `indicator_psa = NULL` — all PSA (backward compat)
- Confirm binary mode behavior unchanged

#### Phase 2: Add `input_block()` function

**File**: `R/input_f.R` (new function, placed after `pick_val_v`)

**Signature** (explicit arguments, user declares what's what):
```r
input_block(
  .data = NULL,          # for pipe chaining (prior add_item block)
  base,                  # list of base case values
  psa,                   # unevaluated expression for PSA draws (e.g., pick_psa(...))
  sens,                  # list OR symbol referencing sensitivity columns
  names_out,             # character vector of parameter names
  psa_indicators = NULL, # list of 0/1 (or per-element vectors) — which params in PSA
  dsa_indicators = NULL, # list of integers (grouped) or NULL (binary one-at-a-time)
  distributions = NULL,  # list of distribution names (needed for grouped mode mvrnorm/dirichlet)
  covariances = NULL     # list of covariance matrices (only for mvrnorm/dirichlet)
)
```

**What it does at creation time**:
1. Validates structure (lengths match, indicators make sense)
2. Determines mode: `dsa_indicators` present → grouped; absent → binary
3. Computes block metadata:
   - `n_params`: number of parameters in this block
   - `n_groups`: number of unique DSA groups (grouped mode) or same as n_params (binary)
   - `mode`: "grouped" or "binary"
4. Builds an unevaluated `{ pick_val_v(...) }` expression with symbolic references to engine variables (`psa_bool`, `sensitivity_bool`, `sens_name_used`, `sens`, `n_sensitivity`)
5. In binary mode, the generated expression includes the `create_indicators()` call with a placeholder for `n_elem_before` (engine patches this)
6. In grouped mode, the generated expression passes `dsa_indicators` directly + `sens_iterator(sens, n_sensitivity)`
7. Attaches metadata as `attr(result, "warden_block_meta")`
8. Returns combined `{ }` expression (piped with `.data` if provided)

**What it does NOT do**: Evaluate anything. No class dispatch in hot loop. The attribute is read once by `run_sim`.

#### Phase 3: Engine modifications in `run_sim` / `run_sim_parallel`

**Files**: `R/run_sim.R`, `R/run_sim_parallel.R`

**Changes at the start of `run_sim()` or `run_sim_parallel()` (before loop)**:
1. Scan all `*_inputs` arguments for `"warden_block_meta"` attributes
2. If found:
   - Compute total `n_sensitivity` across all blocks (sum of `n_groups`)
   - For binary mode blocks: assign `n_elem_before` offsets based on order, patch the expression (replace placeholder with actual offset)
   - Generate `sensitivity_inputs` (the `i_sensitivity` block) automatically — contains `sens_iterator()` call
   - User no longer needs to pass `n_sensitivity` or `sensitivity_inputs` (still accepted for backward compat / manual override)
3. If no metadata found: old behavior, user manages everything manually

**Backward compatibility**:
- If user passes `n_sensitivity` explicitly, it takes precedence (override)
- If user passes `sensitivity_inputs` explicitly, it takes precedence
- Old `add_item() + pick_val_v()` chains still work (no metadata attribute → old path)

#### Phase 4: Documentation and vignette update

- Roxygen docs for `input_block()`
- Update `vignettes/articles/inputs_selector.Rmd` with simplified examples
- Add `input_block` to `_pkgdown.yml`
- NEWS.md bullet

---

### Open questions before implementation

1. For `input_block(..., sens = ?)`: The `sens` column changes at runtime (e.g., `l_inputs[["DSA_min"]]` vs `l_inputs[["DSA_max"]]`). Should the user pass the full list and sensitivity column names, so the expression generates `l_inputs[[sens_name_used]]`? Or pass a symbol? Proposed: user passes the full data list + the engine uses `sensitivity_names` from `run_sim` to know which columns to reference.

2. For binary mode with multi-level offsets: the engine patches the expression after computing offsets. Is expression manipulation (substituting the `n_elem_before` placeholder) acceptable, or should we use a different mechanism (e.g., store offset in the eval environment)?

3. Should `input_block` also accept a `level` argument (e.g., "common_all", "common_pt") for self-documentation, or is the assignment to the `run_sim` argument sufficient?
