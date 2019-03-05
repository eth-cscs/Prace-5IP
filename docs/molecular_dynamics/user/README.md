<nav aria-label="breadcrumb">
  <ol class="breadcrumb">
    <li class="breadcrumb-item"><a href="/molecular_dynamics/">Home</a></li>
    <li class="breadcrumb-item"><a href="/molecular_dynamics/dev/">Developer API</a></li>
    <li class="breadcrumb-item active" aria-current="page">User API</a></li>
    <li class="breadcrumb-item"><a href="/molecular_dynamics/gromacs/">GROMACS refactor</a></li>
  </ol>
</nav>

# User-Facing Code (app-agnostic)

The user facing API is meant to be used in custom user-written MD codes as a non-bonded force calculator. The user is meant to be a domain scientist with minimal knowledge of C++, which is willing to prototype a MD dynamics application with minimal effort. Below we show two use cases where NBLib can be used independently of any library implementation of the API, let it the GROMACS, Amber, LAMMPS or any other custom implementation, in the same spirit as Scalapack routines can be called for basic matrix operations.

Library implementations of this API are open to make use of any internal data structure that fits best the hardware and parallel programming paradigm supported. For now, the user-facing API provides simple input and output data structures that can be easily converted to more complex formats. Further extensions to support different input and output data structures are planned in the future, based on actual library implementations of this API, which would in turn decrease the input/output data conversion.

The user-facing API contains the

* NBICalculator
* Coordinates
* NBParams
* ForceOuput
* EnergyGroups

## The `NBICalculator`

The `NBICalculator` is the entry point for the user-facing API. It encompasses all the basic functionalities that a force calculator can perform, such as, compute forces, produce the forces and/or virial tensors and retrieve energy groups, given a set of coordinates and non-bonded force field parameters. Additionally, its constructor accepts a hardware related `enum` which controls in which hardware the forces should be computed, in case, the underlying hardware and library support different accelerators.

```c++
class NBICalculator
{
public:
    enum Accelerator { CPU, GPU };

    NBICalculator(Coordinates coords, NBParams params,
                  NBICalculator::Accelerator accelerator = NBICalculator::CPU);

    ForceOutput calculateForces();

    EnergyGroups getEnergyGroups();

    void setCoordinates();
}
```

> TODO: add text about the flexibility and portability of the `NBICalculator` focusing on how much it can help prototyping.

## The `Coordinates` object

> TODO: add text about the flexibility and portability of the `NBICalculator` focusing on how much it can help prototyping.


```c++
class Coordinates<T>
{
public:
    Coordinates();

    void setCoordinates(std::vector<float> &coords);
    void setCoordinates(std::vector<float> &x, std::vector<float> &y, std::vector<float> &z);

    void getCoordinates(std::vector<float> &coords);
    void getCoordinates(std::vector<float> &x, std::vector<float> &y, std::vector<float> &z);
}

class Coordinates<Simple>
{
public:
    Coordinates();

    void setCoordinates(std::vector<float> &coords);
    void setCoordinates(std::vector<float> &x, std::vector<float> &y, std::vector<float> &z);

    void getCoordinates(std::vector<float> &coords);
    void getCoordinates(std::vector<float> &x, std::vector<float> &y, std::vector<float> &z);

private:
    std:std::vector<float> coord;
}
```


```c++
class NBParams
{
public:
    NBParams();

    void setSRType(enum::nbSRtype type,
                   float rlist,
                   float rcut,
                   float epsilon_r,
                   float epsilon_rzero,
                   ...);
    void setLRType(enum::nbLRtype type,
                   float rlist,
                   float rcut,
                   float epsilon_r,
                   float epsilon_rzero,
                   ...);

    void setEnergyGroups(std::vector<float> atomIndexes);
    void setLJInteractionMatrix(std::vector<float> c6, std::vector<float> c12);
    void setParticleCharges(std::vector<float> q);
    void setParticleMasses(std::vector<float> m);
}
```

```c++
class ForceOutput
{
public:
    ForceOutput();
    std::vector<Matrix> getAtomicVirialArray();
    std::vector<float> getForceArray();
}

class EnergyGroups
{
public:
    ForceOutput();
    std::vector<float> getEnergyGroups();
}
```






```c++
// User provided functions
void integrate(Coordinates<Simple> &coord, std::vector<float>&velocities,
               ForceOutput const &forces, std::vector<float> const &masses);

ForceOutput restraint(Coordinates<Simple> &coord, EnergyGroups &energy_groups);

int main {

    // Input Data Feed
    InputData inputs(inputfilename);

    Coordinates<simple> coord = inputs.getCoordinates();
    std::vector<float> velocities = inputs.getVelocities();
    std::vector<float> masses = inputs.getMasses();

    NBParams NBParams = inputs.getNBParams();
    // End of Data Feed

    // Auxiliary data structures
    ForceOutput forces;

    // NB lib API
    NBICalculator<nblib> nbcalc(coord, NBParams);

    bool withVirial = false;
    for (int step = 0; i < nsteps; i++)
    {
        if(i % 10 == 0) {
            withVirial = true;
        } else {
            withVirial = false;
        }

        forces = nbcalc.calculateForces(withVirial);
        energy_groups =  nbcalc.getEnergyGroups();
        forces += restraint(coord, energy_groups);

        Coordinates<simple> current_pos = nbcalc.getCoordinates();

        integrate(coord, velocities, forces, masses);

        nbcalc.setCoordinates(current_pos);
    }
}
```

```c++

// User provided functions
void integrate(Coordinates<Simple> &coord, std::vector<float>&velocities,
               ForceOutput const &forces, std::vector<float> const &masses);

void integrate_constrained_atoms(Coordinates<Simple> &coord, std::vector<float>&velocities,
               ForceOutput const &forces, std::vector<float> const &masses);
ForceOutput const &forces computeBondedForces(Coordinates<Simple> &coord,
                                              BondedTopology &  bondedTopology,
                                              bool withVirial);


int main {

    // Input Data Feed
    InputData inputs(inputfilename);

    Coordinates<simple> coord = inputs.getCoordinates();
    std::vector<float> velocities = inputs.getVelocities();
    std::vector<float> masses = inputs.getMasses();
    BondedTopology bondedTopology = inputs.getBondedTopology();

    NBParams NBParams = inputs.getNBParams();
    // End of Data Feed

    // Auxiliary data structures
    ForceOutput forces;

    // NB lib API
    NBICalculator<nblib> nbcalc(coord, NBParams);

    bool withVirial = false;
    for (int step = 0; i < nsteps; i++)
    {
        if(i % 10 == 0) {
            withVirial = true;
        } else {
            withVirial = false;
        }

        forces = nbcalc.calculateForces(withVirial);
        energy_groups =  nbcalc.getEnergyGroups();
        forces += restraint(coord, energy_groups);

        Coordinates<simple> current_pos = nbcalc.getCoordinates();

        // multi-step integration to avoid constraint calculation
        for (int i = 0; i < 10; ++i)
        {
            bonded_forces = computeBondedForces(coord, bondedTopology, withVirial);
            integrate_constrained_atoms(coord, velocities,forces, masses);
        }
        integrate(coord, velocities, forces, masses);

        nbcalc.setCoordinates(current_pos);
    }
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


<!-- TODO: write an appropriate name -->
## Use case A

> Run vanila simualtion for nSteps and extract information MD info

```c++
int doMD(int nSteps) {
    // Get NB Calculator
    NBICalculator   nbcalc(system, "gromacs");

    // run md
    for (i = 1:nSteps)
    {
        ...
        forces = nbcalc.calculateForces(particles, topology, simulaton);
        ...
        // Integration
        ...
        system.update(/* new particle coordinates */);
    }
}

int main(int argc, char *argv[]) {
    // Input Data Feed
    ...
    Particles       particles(x, m, q, v);
    Topology        topology(bondTables, localTopology);
    Simulation      simulaton(parameters, flags, box);
    System          system(particles, topology, simulation);

    doMD(1000);

    my_new_analyses_tool(system, particles, topology, simulation);
}
```

<!-- TODO: write an appropriate name -->
## Use case B

> Investigate other aspects of the MD work flow, e. g., compute prototype new contraint method

```c++
int main(int argc, char *argv[]) {
    // Input Data Feed
    ...
    Particles       particles(x, m, q, v);
    Topology        topology(bondTables, localTopology);
    Simulation      simulaton(parameters, flags, box);
    System          system(particles, topology, simulation);

    // Get NB Calculator
    NBICalculator   nbcalc(system, "gromacs");

    // run md
    for (i = 1:nSteps)
    {
        ...
        forces = nbcalc.calculateForces(particles, topology, simulaton);
        ...
        // Integration
        ...
        computeMyNewConstraintForces(particles, topology, simulaton);
        system.update(/* new particle coordinates */);
    }
}
```

<!-- TODO: write an appropriate name -->
## Use case C

> Investigate other aspects of the MD work flow, e. g., add new external force


```c++
int main(int argc, char *argv[]) {
    // Input Data Feed
    ...
    Particles       particles(x, m, q, v);
    Topology        topology(bondTables, localTopology);
    Simulation      simulaton(parameters, flags, box);
    System          system(particles, topology, simulation);

    // Get NB Calculator
    NBICalculator   nbcalc(system, "gromacs");
    IForceProvider myNewExternalForce(/* TODO: fill in params */):

    nbcalc.registerNewForceProvider(myNewExternalForce);

    // run md
    for (i = 1:nSteps)
    {
        ...
        forces = nbcalc.calculateForces(particles, topology, simulaton);
        ...
        // Integration
        ...
        system.update(/* new particle coordinates */);
    }
}
```



<!-- TODO: write an appropriate name -->
## Use case D

> Investigate other force field equations


```c++
int main(int argc, char *argv[]) {
    // Input Data Feed
    ...
    Particles       particles(x, m, q, v);
    Topology        topology(bondTables, localTopology);
    Simulation      simulaton(parameters, flags, box);
    System          system(particles, topology, simulation);

    // Get NB Calculator
    NBICalculator   nbcalc(system, "gromacs");
    IForceProvider myNewForceField(/* TODO: fill in params */):

    nbcalc.registerForceFieldProvider(myNewForceField);

    // run md
    for (i = 1:nSteps)
    {
        ...
        forces = nbcalc.calculateForces(particles, topology, simulaton);
        ...
        // Integration
        ...
        system.update(/* new particle coordinates */);
    }
}
```


