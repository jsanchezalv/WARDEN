# Create a discrete resource

Creates a discrete resource management system for discrete event
simulations. This system manages a fixed number of identical resource
units that can be blocked (used) by patients and maintains a priority
queue for waiting patients.

## Usage

``` r
resource_discrete(
  n,
  discipline = "FIFO",
  max_queue = Inf,
  allow_multiple_queue = TRUE
)
```

## Arguments

- n:

  Integer \>= 0. Total capacity of the resource.

- discipline:

  `"FIFO"` (default) or `"LIFO"`. Queue discipline applied within the
  same priority level.

- max_queue:

  Non-negative integer or `Inf` (default). Maximum number of patients
  that can wait in the queue. Patients that arrive when the queue is
  full are rejected (`attempt_block()` returns `NA`).

- allow_multiple_queue:

  Logical (default `TRUE`). If `FALSE`, a patient that already has a
  queued entry will be rejected (`attempt_block()` returns `NA`) instead
  of being allowed to queue again.

## Value

An environment of class `"resource_discrete"` with the following
methods:

- `size()`:

  Total capacity.

- [`queue_size()`](https://jsanchezalv.github.io/WARDEN/reference/queue_size.md):

  Number of patients currently in the queue.

- `n_free()`:

  Number of free resource units.

- `n_using()`:

  Number of units currently in use.

- `utilization()`:

  Fraction of capacity in use (0–1).

- `patients_using()`:

  Vector of patient IDs currently using the resource.

- `patients_using_times()`:

  Vector of start times for patients using the resource.

- `queue_start_times()`:

  Vector of queue start times in queue order.

- `queue_priorities()`:

  Vector of priorities in queue order.

- `queue_info(n)`:

  Data frame with `patient_id`, `priority`, `start_time` for the queue.

- `is_patient_in_queue(patient_id)`:

  Check if patient is in queue. Defaults to `i` from the calling
  environment.

- `is_patient_using(patient_id)`:

  Check if patient is using the resource. Defaults to `i` from the
  calling environment.

- `attempt_block(patient_id, priority, start_time, amount)`:

  Attempt to block `amount` units. Returns `TRUE` (acquired), `FALSE`
  (queued), or `NA` (rejected). Defaults `patient_id` to `i` and
  `start_time` to `curtime` from the calling environment.

- `attempt_free(patient_id, remove_all, amount)`:

  Free the resource for a patient. Defaults `patient_id` to `i`.

- `attempt_free_if_using(patient_id, remove_all)`:

  Free only if patient is currently using. Defaults `patient_id` to `i`.

- `next_patient_in_line(n)`:

  Get next `n` patients in queue.

- `modify_priority(patient_id, new_priority)`:

  Modify patient priority in queue. When the patient has `k` active
  queue entries, complexity is O(k log n), not O(log n).

- `add_resource(n)`:

  Add `n` resource units.

- `remove_resource(n, current_time)`:

  Remove `n` resource units. Defaults `current_time` to `curtime`.

- `queue_wait_time(patient_id)`:

  Final queue wait time, set on successful dequeue. Returns `NA` if
  never queued or still waiting. Defaults `patient_id` to `i`.

- `queue_wait_time_current(patient_id, current_time)`:

  Current elapsed queue wait: `0` at entry, grows each event while
  waiting, returns final wait on acquisition, `NA` if never queued.
  Defaults to `i` / `curtime`.

- `had_to_queue(patient_id)`:

  Returns `1L` if patient ever queued, `0L` otherwise. Defaults
  `patient_id` to `i`.

- `time_in_use(patient_id, current_time)`:

  Time since patient acquired resource. Returns `NA` if not currently
  using. Defaults to `i` / `curtime`.

- `total_patients_blocked()`:

  Count of unique patients that ever successfully acquired.

- `total_patients_queued()`:

  Count of unique patients that ever entered the queue.

- `batch_seize(patient_ids, priority, start_time, amount_each)`:

  Seize for multiple patients in a single C++ call. Returns an integer
  vector (`1` = acquired, `0` = queued, `-1` = rejected).

## Examples

``` r
# Create a resource with 3 units
beds <- resource_discrete(3)

# Check initial state
beds$size()       # 3
#> [1] 3
beds$n_free()     # 3
#> [1] 3
beds$queue_size() # 0
#> [1] 0

# Block resources (reads i and curtime from calling environment)
i <- 101L; curtime <- 0.0
beds$attempt_block()
#> [1] TRUE

# Or explicitly
beds$attempt_block(patient_id = 102L, priority = 1L, start_time = 1.0)
#> [1] TRUE

# LIFO resource and limited queue
stack <- resource_discrete(2, discipline = "LIFO", max_queue = 5)
```
