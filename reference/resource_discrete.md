# Create a Discrete Resource

Creates a discrete resource management system for discrete event
simulations. This system manages a fixed number of identical resource
units that can be blocked (used) by patients and maintains a priority
queue for waiting patients.

## Usage

``` r
resource_discrete(n)
```

## Arguments

- n:

  Integer. The total capacity of the resource (must be \>= 1).

## Value

An environment with methods for resource management.

## Details

The returned environment has the following methods:

- `size()`: Returns the total capacity

- [`queue_size()`](https://jsanchezalv.github.io/WARDEN/reference/queue_size.md):
  Returns the number of patients in queue

- `n_free()`: Returns the number of free resource units

- `patients_using()`: Vector of patient IDs currently using the resource

- `patients_using_times()`: Vector of start times for patients using the
  resource

- `queue_start_times()`: Vector of queue start times parallel to queue
  order

- `queue_priorities()`: Vector of priorities parallel to queue order

- `queue_info(n)`: Data.frame with patient_id, priority, start_time for
  queue

- `is_patient_in_queue(patient_id)`: Check if patient is in queue

- `is_patient_using(patient_id)`: Check if patient is using resource

- `attempt_block(patient_id, priority, start_time)`: Attempt to block a
  resource unit

- `attempt_free(patient_id, remove_all)`: Free a resource unit

- `attempt_free_if_using(patient_id, remove_all)`: Free only if patient
  is using

- `next_patient_in_line(n)`: Get next n patients in queue

- `modify_priority(patient_id, new_priority)`: Modify patient priority
  in queue

- `add_resource(n)`: Add n resource units to total capacity

- `remove_resource(n, current_time)`: Remove n resource units from total
  capacity

## Examples

``` r
# Create a resource with 3 units
beds <- resource_discrete(3)

# Check initial state
beds$size()      # 3
#> [1] 3
beds$n_free()    # 3
#> [1] 3
beds$queue_size() # 0
#> [1] 0

# Block resources
i <- 101; curtime <- 0.0
beds$attempt_block()  # Uses i and curtime from environment
#> [1] TRUE

# Or explicitly
beds$attempt_block(patient_id = 102, priority = 1, start_time = 1.0)
#> [1] TRUE

# Check patient status
beds$is_patient_using(101)     # TRUE
#> [1] TRUE
beds$is_patient_in_queue(102)  # FALSE
#> [1] FALSE
```
