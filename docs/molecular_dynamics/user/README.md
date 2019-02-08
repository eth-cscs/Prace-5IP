<nav aria-label="breadcrumb">
  <ol class="breadcrumb">
    <li class="breadcrumb-item"><a href="/molecular_dynamics/">Home</a></li>
    <li class="breadcrumb-item"><a href="/molecular_dynamics/dev/">Developer API</a></li>
    <li class="breadcrumb-item active" aria-current="page">User API</a></li>
    <li class="breadcrumb-item"><a href="/molecular_dynamics/gromacs/">GROMACS refactor</a></li>
  </ol>
</nav>

### User-Facing Code (app-agnostic)

This should be usable outside of GROMACS in custom MD codes as a non-bonded force calculator. In the same spirit as Scalapack routines can be called for basic matrix operations.

```c++
// Input Data Feed
.
.
Particles 		particles(x, m, q, v);
Topology 		topology(bondTables, localTopology);
Simulation 		simulaton(parameters, flags, box);
System 			system(particles, topology, simulation);
.
.
NBICalculator 	nbcalc(system, "gromacs");
.
for (i = 1:nSteps)
{
    ...
    forces = nbcalc.calculateForces();
    ...
    // Integration
    ...
    system.update(/* new particle coordinates */);
}
```

The interface only consists of data and objects relevant to the simulation from the user. User would basically read input and topology files and setup the simulation without being concerned about the underlying library.

The NBCalculator API would be external aiming to encapsulate the implementation details of the application from the end user. It caters to those who need a non-bonded force calculator much the same way as one might need a linear-system solver.

## Reusable Utilities for Schedules

The current implementation of `do_forcecutsVERLET(...)` performs many operation sequences that could serve as reusable utilities in other custom schedule as these would always be executed together.

The other benefit of doing so is to greatly improve code readability in the schedule and thereby encapsulate implementation details for developers assembling a customized schedule. It would also help shrink the lines of code and reduce code duplication.

Some examples can be:

```c++
void initShiftVectors(...);
void initPMEonDevice(...);
void
```

It can make use of objects that are pointed by the members of `ForceSchedule` class.

```c++
void doNSGridding() {
    if (neighbourSearch) {
        auto shiftVec = calcShiftVectors(system, bondGraph);
        system.applyTransformation(shiftVec);

        if (domainDecomposition) {
            // scattter particle data across ranks before gridding
			putParticlesOnGridNonLocal(...);
        }
        else {
            // perform gridding locally
            putParticlesOnGrid(...);
        }
        setAtomData(...);
    }
}
```
