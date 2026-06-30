## R package development

### Key commands

```
# To run code
"/c/Program Files/R/R-4.2.3/bin/Rscript.exe" -e "devtools::load_all(); code"

# To run all tests
"/c/Program Files/R/R-4.2.3/bin/Rscript.exe" -e "devtools::test()"

# To run all tests for files starting with {name} (e.g., 'input_f')
"/c/Program Files/R/R-4.2.3/bin/Rscript.exe" -e "devtools::test(filter = '{name}')"

# To redocument the package
"/c/Program Files/R/R-4.2.3/bin/Rscript.exe" -e "devtools::document()"

# To check pkgdown documentation
"/c/Program Files/R/R-4.2.3/bin/Rscript.exe" -e "pkgdown::check_pkgdown()"

# To check the package with R CMD check
"/c/Program Files/R/R-4.2.3/bin/Rscript.exe" -e "devtools::check()"

```

### Coding

* Use the base pipe operator (`|>`) not the magrittr pipe (`%>%`) whenever possible.
* Use `\() ...` for single-line anonymous functions. For all other cases, use `function() {...}` 

### Testing

- New tests for `R/{name}.R` go in `tests/testthat/test-{name}.R`. 
- All new code should have an accompanying test.
- If there are existing tests, place new tests next to similar existing tests.
- Strive to keep your tests minimal with few comments.

### Documentation

- Every user-facing function should be exported and have roxygen2 documentation.
- Wrap roxygen comments at 80 characters.
- Whenever you add a new (non-internal) documentation topic, also add the topic to `_pkgdown.yml`. 
- Always re-document the package after changing a roxygen2 comment.
- Use `pkgdown::check_pkgdown()` to check that all topics are included in the reference index.

### `NEWS.md`

- Every user-facing change should be given a bullet in `NEWS.md`. Do not add bullets for small documentation changes or internal refactorings.
- Each bullet should briefly describe the change to the end user and mention the related issue in parentheses.
- A bullet can consist of multiple sentences but should not contain any new lines (i.e. DO NOT line wrap).
- If the change is related to a function, put the name of the function early in the bullet.
- Order bullets alphabetically by function name. Put all bullets that don't mention function names at the beginning.

### GitHub

- Do not make any interactions with Github or Git, leave those to me.

### Writing

- Use sentence case for headings.
- Use US English.

### Proofreading

If the user asks you to proofread a file, act as an expert proofreader and editor with a deep understanding of clear, engaging, and well-structured writing. 

Work paragraph by paragraph, always starting by making a TODO list that includes individual items for each top-level heading. 

Fix spelling, grammar, and other minor problems without asking the user. Label any unclear, confusing, or ambiguous sentences with a FIXME comment.

Only report what you have changed.

### Workflow Guidance

Populate and maintain .claude/CLAUDE.md within the with all relevant project-wide context so you can resume work efficiently without me repeating context each session. Include:
- Project summary & active features
- Tech stack
- Location of all functions in each file within R/ folder for easy retrieval with a one line summary
- Code style & naming conventions
- Known bugs and next TODOs
- Test scenarios we haven’t completed yet (if any)
Keep it under 5k tokens total.

The workflow will always follow the following steps: 1) design, with questions, brainstorming, 2) general specifications, 3) detailed implementation plan (saved within progress.md), 4) plan-review against specs, 4) test writing, 5) implementation, 6) testing. If test fails, modify implementation until tests work (do not modify tests at this stage without my explicit permission), 7) code-review against plan.

Do not add additional dependencies without explicitly asking me so. Speed is vital in this code, make sure your implementation tries to be as efficient as possible. 

Write high quality code as if the code failing could imply human life losses. Code should be human readable and clear, avoid ai-ism and wrapping too many functions, make the code as simple and intuitive as possible while keeping the other criteria true.

---

### Project summary & active features

**WARDEN** (Workflows for Health Technology Assessments in R using Discrete EveNts) v2.0.2 is an R package for discrete event simulation (DES) in health technology assessments (HTA). It supports cost-effectiveness modelling aligned with NICE TSD 15.

Key features:
- Standard DES engine with Rcpp-based event queue (priority queue per patient)
- Resource-constrained DES engine (shared resources/inputs across patients within an arm)
- Probabilistic sensitivity analysis (PSA), deterministic sensitivity analysis (DSA), and scenario analysis
- Parallel simulation via `run_sim_parallel()` using `future`/`doFuture`
- Rcpp implementations of discounting, conditional quantile, and luck adjustment functions
- Debug mode and continue-on-error functionality
- `shared_input()` and `resource_discrete()` for constrained simulations
- `random_stream()` for reproducible random number handling

### Tech stack

- **Language**: R (>= 2.10) with Rcpp (C++ extensions compiled via `src/`)
- **Key imports**: `data.table`, `purrr`, `foreach`, `future`, `doFuture`, `progressr`, `flexsurv`, `MASS`, `zoo`, `tidyr`, `lifecycle`, `magrittr`
- **Suggests**: `dplyr`, `ggplot2`, `knitr`, `rmarkdown`, `kableExtra`, `testthat (>= 3.0.0)`, `survival`
- **Tests**: testthat (edition 3), 4 test files covering inputs, queues, resource_discrete, and model runs
- **C++ source files**: `src/evt_queue.cpp`, `src/condq.cpp`, `src/adj_luck.cpp`, `src/disc_cycle_v.cpp`, `src/disc_instant_v.cpp`, `src/disc_ongoing_v.cpp`, `src/resource_constrained.cpp`, `src/rpoisgamma_rcpp.cpp`

### Function Index

The function index is stored in `.claude/function_index.md`. Make sure to update that file whenever you make changes to the code so the lines are always up to date.

### Test files

- `tests/testthat/test-input_f.R` — tests for `R/input_f.R` (add_item, queue functions, event wrappers, etc.)
- `tests/testthat/test-resource_discrete.R` — tests for `resource_discrete()` and constrained DES
- `tests/testthat/test-queues.R` — tests for event queue operations (new_event, modify_event, etc.)
- `tests/testthat/test-model_runs.R` — integration tests for full model runs (run_sim, run_sim_parallel)
- `tests/testthat/test-regression_ssd.R` — regression tests for SSD vignette (det/DSA/probDSA/PSA); 14 tests; `skip_on_cran()`
- `tests/testthat/test-regression_ssd_constrained.R` — regression tests for constrained SSD vignette (constrained/unconstrained/unbinding/DSA/probDSA/PSA); 18 tests; `skip_on_cran()`
- `tests/testthat/test-regression_eBC.R` — regression test for early breast cancer vignette (deterministic); 5 tests; `skip_on_cran()`
- `tests/testthat/test-regression_inputs_selector.R` — regression tests for inputs_selector vignette (det/PSA/DSA/scenario/split/cov/vec); 27 tests; `skip_on_cran()`

**Key constraint for regression test deferred expressions**: `add_item()`/`add_tte()`/`add_reactevt()` capture expressions unevaluated via `substitute()`. When `run_sim()` evaluates them, it looks up names in sim_env → `.GlobalEnv` only. Test-file-level and test-block-level variables are NOT accessible. All data values must be inlined directly in the deferred expressions; use `list(base=..., DSA_min=..., DSA_max=...)[[sens_name_used]]` for sensitivity dispatch.

### Code style & naming conventions

- Use base pipe `|>` not `%>%`
- Use `\() ...` for single-line anonymous functions; `function() {...}` for multi-line
- File naming: `*_f.R` for function files (e.g., `input_f.R`, `calculator_f.R`)
- Test files: `test-{name}.R` mirrors `R/{name}.R` (e.g., `test-input_f.R`)
- Rcpp wrapper functions end with `_cpp` (internal); R-level wrappers drop that suffix
- Vectorised Rcpp versions of discounting functions end with `_v` (e.g., `disc_cycle_v()`)
- `_old` suffix = legacy R implementations superseded by Rcpp
- Sensitivity/PSA iterations use `sens` (integer index) throughout

### Known bugs & next TODOs

- Update parallel engine with **Mori** package
- Update `"i"` iterator with `pat_i` to make it clearer and less bug-prone
- No known bugs remaining in resource engine (all identified issues from v2.0.5 review resolved in v2.0.5 patch)

### Test scenarios not yet completed

- No dedicated test file for `R/calculator_f.R` (distributional utilities untested via formal tests)
- No dedicated test file for `R/results_summary_f.R` (CEAC, EVPI, summary functions)
- No tests for `run_sim_parallel()` PSA mode edge cases
- No tests for `adj_val()` with `vectorized_f = TRUE`
- No tests for `qtimecov()` or `luck_adj()` in isolation