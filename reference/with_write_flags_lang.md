# Evaluate a language object with write-tracking via active bindings

Temporarily installs **active bindings** for a set of variable names in
`env` so that any **write** to those names during the evaluation of
`expr_lang` flips a corresponding `*_lastupdate` flag in `env`. This
flags assignments even when the value is overwritten with the same
value, and does not flag branches that are not executed.

## Usage

``` r
with_write_flags_lang(expr_lang, tracked, env, flag_value = 1L)
```

## Arguments

- expr_lang:

  A language object (call or expression) to evaluate. For example,
  `react_list[[i]][["react"]]` returned by
  [`add_reactevt()`](https://jsanchezalv.github.io/WARDEN/reference/add_reactevt.md).

- tracked:

  Character vector of **top-level** variable names to track (e.g.,
  `c("q_default","c_default","curtime")`). For each name `nm`, an active
  binding is installed at `env[[nm]]`, and writes set
  `env[[paste0(nm, "_lastupdate")]] <- flag_value`.

- env:

  Environment in which to install the bindings and evaluate `expr_lang`
  (your per-patient/per-arm working environment).

- flag_value:

  Scalar written into each `*_lastupdate` when a write occurs (commonly
  `1L`; you may use a timestamp like `env$curtime`).

## Value

Invisibly returns `NULL`. The function works by side effects (mutating
`env` and its `*_lastupdate` flags).

## Details

Typical use: when `accum_backwards = TRUE`, wrap the evaluation of a
reaction/input block so you can later accumulate only inputs that were
actually written in that event.

- The function creates a private backing store for `tracked` names and
  installs active bindings in `env`. Reads return the stored value.
  Writes (via `<-` or `=`) update the store **and** set the
  corresponding `*_lastupdate` flag in `env`.

- It **does not** zero the `*_lastupdate` flags; do that before calling
  this function (typically once per event).

- Only **top-level symbols** are tracked. To track `obj$el <- ...`,
  track `"obj"` (the container) rather than `"obj$el"`.

- Active bindings are removed and plain symbols restored on exit, even
  if an error occurs (via
  [`on.exit()`](https://rdrr.io/r/base/on.exit.html) teardown).

## Performance

Active bindings add a tiny overhead per read/write of the tracked names
(a function call). Keep `tracked` small and the binding body minimal.
Install only for the duration of the block and tear down immediately
(handled for you). If your block performs many reads and few writes,
consider a snapshot-and-diff approach instead
