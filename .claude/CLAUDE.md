## R package development

### Key commands

```
# To run code
Files/R/R-4.2.3/bin/Rscript.exe -e "devtools::load_all(); code"

# To run all tests
Files/R/R-4.2.3/bin/Rscript.exe -e "devtools::test()"

# To run all tests for files starting with {name}
Files/R/R-4.2.3/bin/Rscript.exe -e "devtools::test(filter = '^{name}')"

# To redocument the package
Files/R/R-4.2.3/bin/Rscript.exe -e "devtools::document()"

# To check pkgdown documentation
Files/R/R-4.2.3/bin/Rscript.exe -e "pkgdown::check_pkgdown()"

# To check the package with R CMD check
Files/R/R-4.2.3/bin/Rscript.exe -e "devtools::check()"

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

- If you use `gh` to retrieve information about an issue, always use `--comments` to read all the comments.
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

If it becomes impossible to keep all relevant info under 5k tokens, split less critical sections into separate files under the .claude/ directory. For example, if there are details about a future system version (e.g., MK2), create a separate markdown file like .claude/CLAUDE-2.md instead of bloating CLAUDE.md.

Don't do large indexing files, as gargabe-collection occurs on the connection and exclude .html files.
options:
  files_search:
    extensions: ["R", "Rmd", "qmd"]
    exclusions: ["html", "docs/",".git/",".github/"]
    
Avoid using "btw_tool_files_search" as it takes a super long time to search and it does not work well, so just grab directly the files and read them.

    
### TODO

- Update parallel engine with Mori package

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
- **Documentation**: roxygen2 + pkgdown website at <https://jsanchezalv.github.io/WARDEN/>
- **Tests**: testthat (edition 3), 4 test files covering inputs, queues, resource_discrete, and model runs
- **C++ source files**: `src/evt_queue.cpp`, `src/condq.cpp`, `src/adj_luck.cpp`, `src/disc_cycle_v.cpp`, `src/disc_instant_v.cpp`, `src/disc_ongoing_v.cpp`, `src/resource_constrained.cpp`, `src/rpoisgamma_rcpp.cpp`

### Function index by file

**`R/run_sim.R`** — Main user-facing simulation runner
- `run_sim()` [L150–523]: Run a DES (standard or constrained) across arms/sensitivities/PSA

**`R/run_sim_parallel.R`** — Parallel simulation runner
- `run_sim_parallel()` [L164–543]: Like `run_sim()` but uses `future`/`doFuture` for parallel execution

**`R/run_engine.R`** — Internal standard DES engine
- `run_engine()` [L20–290]: Internal engine processing events per patient per arm

**`R/run_engine_constrained.R`** — Internal constrained DES engine
- `run_engine_constrained()` [L47–386]: Internal engine for resource-constrained DES

**`R/engine_helper_f.R`** — Internal engine utilities
- `load_inputs()` [L41–44]: Evaluate and load unevaluated input expressions into environment
- `debug_inputs()` [L60–91]: Compare old vs new inputs for debug logging
- `initiate_evt()` [L109–130]: Initialise events for a patient at arm start
- `react_evt()` [L153–194]: Process the next event and apply reactions
- `eval_reactevt()` [L216–277]: Evaluate reaction list for a given event
- `get_input()` [L300–332]: Retrieve input value with type coercion
- `interval_out()` [L356–365]: Format simulation output as confidence interval string
- `.set_last_ctx()` [L367–393]: Internal — update error context beacon
- `log_add()` [L396–402]: Append an entry to the debug log
- `on_error_check()` [L404–416]: Wrap expression for error handling/continue-on-error
- `with_write_flags_lang()` [L467–512] *(exported)*: Track assignments in expression for debug mode
- `transform_debug()` [L527–542]: Reshape raw debug data into display format
- `export_log()` [L557–613]: Write debug log to file
- `expand_evts_bwd()` [L627–683]: Expand event data to fill time points (backward accumulation)
- `expand_evts_fwd()` [L697–753]: Expand event data to fill time points (forward accumulation)
- `compute_outputs_timseq()` [L785–996]: Compute timed-frequency outputs
- `compute_outputs()` [L1035–1392]: Compute all discounted outputs for a patient

**`R/input_f.R`** — Input construction, event queue wrappers, and utilities
- `replicate_profiles()` [L31–47]: Replicate patient profiles across arms
- `create_indicators()` [L67–94]: Create sensitivity/scenario indicator variables
- `sens_iterator()` [L115–120]: Iterate through DSA/scenario combinations
- `pick_psa()` [L192–203]: Draw PSA samples from a distribution function
- `pick_val_v()` [L265–357]: Select values across base/sensitivity/PSA with a vector of options
- `add_item()` [L406–457]: Add named items to the input list (does not support native pipe `|>`)
- `add_item2()` [L475–489] *(deprecated)*: Older item-adding variant (now merged into `add_item()`)
- `modify_item()` [L521–536] *(deprecated)*: Modify items in input list (no longer needed)
- `modify_item_seq()` [L572–597] *(deprecated)*: Sequential item modification (no longer needed)
- `queue_create()` [L610–612]: Create Rcpp-backed priority event queue
- `new_event()` [L639–675]: Add event(s) for a patient to the queue
- `next_event()` [L686–689]: Peek at the next event(s) in the queue
- `next_event_pt()` [L701–706]: Peek at next event(s) for a specific patient
- `pop_event()` [L716–719]: Remove the top event from the queue
- `pop_and_return_event()` [L728–731]: Remove and return the top event
- `remove_event()` [L743–779]: Remove named event(s) for a patient
- `modify_event()` [L806–848]: Modify event time/name for a patient (creates if missing by default)
- `queue_empty()` [L857–860]: Check if queue is empty
- `queue_size()` [L869–872]: Get number of events in queue
- `has_event()` [L883–887]: Check if a patient has a named event
- `get_event()` [L898–902]: Get time of a named event for a patient
- `resource_discrete()` [L960–1159]: Create Rcpp-backed discrete resource object
- `print.resource_discrete()` [L1168–1175]: Print method for `resource_discrete`
- `shared_input()` [L1271–1337]: Create shared input object for constrained DES
- `add_reactevt()` [L1376–1399]: Add reactions (event-triggered logic) to the model
- `random_stream()` [L1439–1465]: Create a random number stream object for reproducibility
- `add_tte()` [L1494–1527]: Add time-to-event draws to the input list
- `adj_val()` [L1573–1613]: Integrate a value over a time interval (supports discounting)
- `disc_ongoing()` [L1635–1651]: Discount ongoing (flow) value between two time points
- `disc_instant()` [L1671–1681]: Discount an instantaneous value at a single time point
- `disc_cycle()` [L1709–1764]: Discount a cycle-based value
- `extract_from_reactions()` [L1799–1805]: Extract item/event references from reaction list
- `ast_as_list()` [L1861–1877]: Convert R AST expression to a nested list
- `extract_elements_from_list()` [L~2000+]: Extract assignments/references from AST list

**`R/calculator_f.R`** — Statistical/distributional utilities
- `draw_tte()` [L45–76]: Draw time-to-event from common parametric distributions (via flexsurv)
- `rdirichlet()` [L96–124]: Sample from a Dirichlet distribution
- `rdirichlet_prob()` [L148–176]: Sample Dirichlet with SE-based parameterisation
- `rbeta_mse()` [L197–208] / `qbeta_mse()` [L226–232]: Beta distribution parameterised by mean and SE
- `rgamma_mse()` [L253–267] / `qgamma_mse()` [L287–297]: Gamma distribution parameterised by mean and SE
- `rcond_gompertz_lu()` [L317–326] / `rcond_gompertz()` [L345–352] / `pcond_gompertz()` [L371–376]: Conditional Gompertz functions
- `rpoisgamma()` [L415–494]: Draw from Poisson-Gamma (negative binomial) distribution
- `cond_mvn()` [L522–564]: Conditional multivariate normal distribution
- `cond_dirichlet()` [L588–616]: Conditional Dirichlet distribution
- `discrete_resource_clone()` [L735–744]: Clone a `resource_discrete` object

**`R/RcppExports.R`** — Auto-generated Rcpp wrappers (do not edit manually)
- `luck_adj()` [L74–76]: Luck adjustment for correlated survival draws
- `qcond_exp()` [L91–93] / `qcond_weibull()` [L108–110] / `qcond_weibullPH()` [L125–127] / `qcond_llogis()` [L142–144] / `qcond_gompertz()` [L159–161] / `qcond_lnorm()` [L180–182] / `qcond_norm()` [L201–203] / `qcond_gamma()` [L222–224]: Vectorised conditional quantile functions
- `qtimecov()` [L378–380]: Time-to-event with time-varying covariates
- `disc_cycle_v()` [L422–424] / `disc_instant_v()` [L439–441] / `disc_ongoing_v()` [L457–459]: Vectorised Rcpp discounting functions
- Queue/resource C++ wrappers [L461–603]: `queue_create_cpp`, `new_event_cpp`, etc.

**`R/results_summary_f.R`** — Results summarisation
- `summary_results_det()` [L52–130]: Summarise deterministic simulation results
- `summary_results_sim()` [L160–242]: Summarise PSA simulation results
- `summary_results_sens()` [L273–371]: Summarise sensitivity/scenario analysis results
- `extract_psa_result()` [L401–406]: Extract a specific element from PSA results list
- `ceac_des()` [L445–482]: Compute cost-effectiveness acceptability curve (CEAC) data
- `evpi_des()` [L518–558]: Compute expected value of perfect information (EVPI)

**`R/old_R_nowincpp.R`** — Legacy R implementations (superseded by Rcpp, kept for reference)
- `luck_adj_old()` [L74–98] / `disc_ongoing_v_old()` [L117–130] / `disc_instant_v_old()` [L146–152] / `disc_cycle_v_old()` [L196–243]: Old vectorised R versions
- `qcond_gompertz_old()` [L260–268] / `qcond_exp_old()` [L285–291] / `qcond_weibull_old()` [L308–314] / `qcond_weibullPH_old()` [L330–340] / `qcond_llogis_old()` [L357–363] / `qcond_lnorm_old()` [L383–389] / `qcond_norm_old()` [L410–416] / `qcond_gamma_old()` [L437–443]: Old conditional quantile functions
- `qtimecov_old()` [L602–698]: Old time-varying covariate TTE function

**`R/data.R`** — Dataset documentation
- `tte.df`: Example time-to-event data frame for package examples

### Test files

- `tests/testthat/test-input_f.R` — tests for `R/input_f.R` (add_item, queue functions, event wrappers, etc.)
- `tests/testthat/test-resource_discrete.R` — tests for `resource_discrete()` and constrained DES
- `tests/testthat/test_queues.R` — tests for event queue operations (new_event, modify_event, etc.)
- `tests/testthat/test_model_runs.R` — integration tests for full model runs (run_sim, run_sim_parallel)

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

From `### TODO` section above:
- Update parallel engine with **Mori** package
- Simplify `pick_val_v()` approach so the user does not need to front-load iterators
- Update `"i"` iterator with `pat_i` to make it clearer and less bug-prone

### Test scenarios not yet completed

- No dedicated test file for `R/calculator_f.R` (distributional utilities untested via formal tests)
- No dedicated test file for `R/results_summary_f.R` (CEAC, EVPI, summary functions)
- No tests for `run_sim_parallel()` PSA mode edge cases
- No tests for `adj_val()` with `vectorized_f = TRUE`
- No tests for `qtimecov()` or `luck_adj()` in isolation