# WARDEN Explained

## Introduction

This document explains the logic behind WARDEN main approach to simulate
discrete event simulations, as well as explaining briefly the rationale
for certain design decisions.

### In a Nutshell

WARDEN main simulation engine at its core is nothing but a nested loop
at different levels. However, for this to work we need to delay the
execution of the inputs provided by the user, so the relevant inputs
provided through
[`add_tte()`](https://jsanchezalv.github.io/WARDEN/reference/add_tte.md),
[`add_item()`](https://jsanchezalv.github.io/WARDEN/reference/add_item.md)
and
[`add_reactevt()`](https://jsanchezalv.github.io/WARDEN/reference/add_reactevt.md)
are substituted for delayed execution and stored as lists.

The basic engine works as indicated below. This is the engine used for
unconstrained DES. For constrained DES, the only difference is that once
all inputs are loaded (analysis, simulation, arm, patient), they are all
saved, and we then iterate over the queue of events in order, loading
the corresponding inputs for the patient for which the event applies,
evaluating the event, and then moving on to the next event (and the
corresponding patient).

1.  **Per Analysis (DSA, scenarios) “sens”**
    1.  Load inputs sequentially and assign to the “sens” input
        environment.
    2.  **Per Simulation (PSA or deterministic) “simulation”**
        1.  Load inputs sequentially and store its components. The
            “sens” environment is integrated into a new environment with
            the “simulation” input environment.
        2.  **Per Patient “i”**
            1.  Load inputs sequentially and store its components. The
                “simulation” environment is integrated into a new
                environment with the “i” input environment
            2.  **Per Arm “arm”**
                1.  Load inputs sequentially and store its components.
                    The “i” environment is integrated into a new
                    environment with the “arm” input environment
                2.  Load initial time to events. First look into the
                    initial time to event expression declared by user;
                    if not found, look into the input list already
                    declared; if not found, set equal to `Inf`. Events
                    are added to an event queue object defined within
                    C++ for efficient management of the queue.
                3.  **While `curtime` (simulation time) is \< `Inf`**
                    1.  Select the next event by checking the event
                        queue; in case of ties, untie using the order
                        declared in `add_tte` for initial time to
                        events. If there are no events left, the event
                        time is `Inf` or if the user sets
                        `curtime = Inf` then the simulation ends.
                    2.  Evaluate the reaction of the event by looking at
                        the relevant expression from the list of event
                        reactions
        3.  Once the specific “simulation” is done, compute outputs
            vectorized (discount outcomes as relevant based on their
            type, aggregate data as relevant, obtain timed frequency
            outputs if needed, etc.)

The debug mode will store in a log the relevant data that is loaded or
changed by the event reactions, and will be exported when the simulation
stops (also on error). WARDEN allows to continue on error (though not
recommended)

WARDEN handles the random numbers automatically, setting the seeds
differently at the simulation, patient and arm level. WARDEN makes sure
that the starting seed is cloned for a patient across interventions.
However, it could be that conditional statements can alter the random
state of R if they conditional trigger random expressions (e.g.,
`if(arm==2){runif(1)}else{5}`) that change per intervention. To keep the
random number cloned as intended, it’s very strongly recommended to
pre-draw random numbers for each type of random object used and use
those (see the `vignette("example_ssd_stream")` vignette for more
information). WARDEN uses L’Ecuyer-CMRG random number generator.

``` r
#Code can be writing directly as an expression which will be evaluated at the right time
add_reactevt(name_evt = "event_1", input = {
  a  <- 1
  z  <- 2
  b  <- z + 1 
  c  <- b + 5 #b will be available 
}
```

### Storing Inputs, Making it Faster

Multiple ways of storing inputs and processing events can be thought of.
A few of these could be 1) data.frames, 2) lists, 3) environments, or 4)
utilize a C++ implementation (among others). WARDEN uses environments to
store inputs, and a queue with unordered maps in C++ to process events.

Data.frames can be slow and memory-intense to manipulate, so they were
avoided for this purpose.

\[Changed with WARDEN 1.0\] The limitation with the debugging mode has
been handled by extract the abstract syntax tree of the event reactions
and looking for any type of assignments. A limitation of this is that
“dynamic” assignments (e.g., `assign(paste0("x_",i), 5)` where `i` is
created by a loop) are NOT captured by the debugging engine, and
therefore will be excluded from the debugging log file. So the user
should try to assign variables explicitly whenever possible, e.g.,
`x_1 <- 5`.

\[Changed with WARDEN 2.0\] A C++ implementation was avoided for a long
time as the purpose of WARDEN is to be user-friendly and to give the
user as much as freedom as possible on how to define their inputs. While
it can make the constrained implementation doable and faster,
implementation in C++ requires careful consideration on how the user can
interact with the object in question. \[Changed with WARDEN 2.0\] In the
most recent version of WARDEN, several functions now use a C++ function
under the hood for speed improvement, but the user will not notice any
change relative to the R counterpart. The core engine has been revamped
so that events now use a C++ implementation of a queue, but it has been
designed so that users can interact as before with the queue using pure,
simple R (so
[`modify_event()`](https://jsanchezalv.github.io/WARDEN/reference/modify_event.md)
is still used, etc). Furthermore, a resource-constrained engine has been
created and a discrete constrained resource object has been created
using C++. Again, the user can interact with these objects using clear R
functions. Speed-wise, for the unconstrained set-up the user is unlikely
to see large speed gains, as there is an implicit cost of setting up
each event queue, modifying events, etc. The constrained approach can
achieve similar speeds to the unconstrained method.

Several objects now in WARDEN use an R6-like interaction, particularly
the random streams
([`random_stream()`](https://jsanchezalv.github.io/WARDEN/reference/random_stream.md)),
restricted resources
([`resource_discrete()`](https://jsanchezalv.github.io/WARDEN/reference/resource_discrete.md))
and shared inputs across patients
([`shared_input()`](https://jsanchezalv.github.io/WARDEN/reference/shared_input.md)).
The reason to select this type of object is due to the type of
interaction the user needs to perform with this objects. For example for
random streams we want to do two things at once: 1) draw a random number
from a pre-generated list of random numbers, and 2) also remove the last
used number from the list. For a discrete resource, we want a patient to
attempt to block the resource, and to also modify at the same time the
discrete resource object (add the patient status to be using the
resource, or to join the queue). This solution makes the code much
easier to implement from a user perspective.

### Parallel engine approach

Furthermore, a parallel core/thread implementation is also available at
the simulation level, i.e., it will perform the “simulation” loop in
parallel. The reason to select the simulation and not the patient is
that each patient normally takes a small amount of time to run, and the
simulation level offers the right balance in terms of time to run.

However, the user should expect it to be only slightly more efficient
(perhaps 20-40% speed increase for medium to large simulations), as
opposed to radically faster. Two factors will be important: the number
of simulations to be run (`n_sim`), and the size of each simulation
(given by the number of events and the number of patients and arms). If
`n_sim` is small, it may not be worth it to use a parallel approach as
there is a time loss to set up the different cores/threads (normally 2
to 5 seconds), so if each simulations runs fast because they are simple
(a couple or seconds or so) it may not be worth it. Even if `n_sim` is
large and each simulation is complex, the efficiency gain may be
~20-40%, even if using \>5 cores. The reason is that RAM use increases
fast as R creates new sessions with duplicated data (it’s not shared
among the cores/threads), and a medium to large simulation can easily
become \>2 GB of RAM use per simulation, so systems with large
processing power AND large RAM (e.g., 32 or 64GB) will benefit the most
from this approach.

The parallel implementation also has limitations in terms of exporting
logs if there is an error in a simulation (due to the parallel set-up),
so this approach is recommended when the user is quite confident that
the simulation will run without issues.
