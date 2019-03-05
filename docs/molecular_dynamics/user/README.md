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

class NBparams
{
public:
    NBparams();

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

class NBICalculator
{
public:
    NBICalculator(Coordinates coords, NBparams params,
                  enum::NBICalculator type = NBICalculator::AUTO);

    ForceOutput calculateForces();

    EnergyGroups getEnergyGroups();

    void setCoordinates();
}


// User provided functions
void integrate(Coordinates<Simple> &coord, std::vector<float>&velocities,
               ForceOutput const &forces, std::vector<float> const &masses);

ForceOutput restraint(Coordinates<Simple> &coord, EnergyGroups &energy_groups);

// addMyThermostat
// addMyBarostat

int main {

    // Input Data Feed
    InputData inputs(inputfilename);

    Coordinates<simple> coord = inputs.getCoordinates();
    std::vector<float> velocities = inputs.getVelocities();
    std::vector<float> masses = inputs.getMasses();

    NBparams nbparams = inputs.getNBparams();
    // End of Data Feed

    // Auxiliary data structures
    ForceOutput forces;

    // NB lib API
    NBICalculator<nblib> nbcalc(coord, nbparams);

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

        // mybarostat.compute(particles, topology, simulation)
        integrate(coord, velocities, forces, masses);

        nbcalc.setCoordinates(current_pos);
    }
}


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

    NBparams nbparams = inputs.getNBparams();
    // End of Data Feed

    // Auxiliary data structures
    ForceOutput forces;

    // NB lib API
    NBICalculator<nblib> nbcalc(coord, nbparams);

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


