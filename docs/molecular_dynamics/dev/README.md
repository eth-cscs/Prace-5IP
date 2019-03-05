<nav aria-label="breadcrumb">
  <ol class="breadcrumb">
    <li class="breadcrumb-item"><a href="/molecular_dynamics/">Home</a></li>
    <li class="breadcrumb-item active" aria-current="page">Developer API</a></li>
    <li class="breadcrumb-item"><a href="/molecular_dynamics/user/">User API</a></li>
    <li class="breadcrumb-item"><a href="/molecular_dynamics/gromacs/">GROMACS refactor</a></li>
  </ol>
</nav>


# API Abstractions Overview

## Developer-Facing Abstraction (application-specific)

GROMACS is the first contender to implement this API on. This needs significant functionality from GROMACS to be exposed. There needs to be an interface that would encapsulate the details of implementation and specific data structures used within GROMACS from a domain scientist dealing with more simplified representation of the system being simulated.

### NBICalculator Abstraction

The constructor would setup the GROMACS-specific data structures, domain decomposition, communicator, load balancing, etc concerning the force-calculations to initialize an `IForceSchedule` object which would then be called by the high level API. The aforementioned abstractions will be introduced in detail in the following section. However, this forms the basis of the user-facing API that aims to connect user provided data with the nuts and bolts of the application in focus.

```c++
class NBICalculator {
private:
    ScheduleBuilder 				scheduleBuilder_; 	// build & init schedules
    std::shared_ptr<IForceSchedule> schedule_; 			// compute forces/energy/potential
    MDBookKeeping 					mdBook_; 			// maintain flags on DD, NS, etc
    
public:
	NBICalculator(Coordinates coords, NBparams params,
                  enum::NBICalculator type = NBICalculator::AUTO) 
    {
        /* init ScheduleBuilder 
        --------------------------*/
        /*
        create GROMACS objects like DD, commrec,
        using the system data
        */
        ...
        scheduleBuilder_.selectSchedule(/* input/execution 			context*/);
        schedule_ = scheduleBuilder_.generateSchedule(/* gmx-specific args */);
        mdBook.init(schedule_, system);
	}
    
    ForceOutput calculateForces() {
		return schedule_->computeStep(mdBook_.flags);
    }
    
	EnergyGroups getEnergyGroups();
    void setCoordinates();
    
};
```

The schedule API underneath would have it's own high performance implementation that overlaps computation, communication and reduction tasks. There would be customized schedules for architecture. algorithms and physics being simulated. This would be a part of GROMACS itself and need not be exposed outside as its implementation details are not so important.

## Force Schedule API

The key idea is to have a collection of static force calculation schedules, that are hand-crafted for performance given the algorithm, execution context and required quantities. An abstract base class summarizes the minimum set of capabilities for a schedule.

```c++
class AbstractForceSchedule {
protected:
    // contains handles to app-specific system and state data structures
    SimulationState     simuStat_;    
	
    // app-specific data structures for internal tasks
    GMX_Communicator    cr_;        
    ForceVector         fr_;
    ForceVirial         vir_;
    PerfCounter         pc_;

    void init(/* app-specific args */);
public:
    virtual void computeStep(/* flags */);

    friend class ScheduleBuilder;
};
```

A collection of concrete schedules would target easy code maintainability, and provide hooks to add new features with minimal interference with the rest of the codebase. The implementation of the `computeStep(...)` function would sequence tasks of data communication, load balancing, computation, and ultimately reduction across nodes and with accelerators. It can optionally make use of private members specific to a schedule.

```c++
class FullGPU_SPME_Schedule : public AbstractForceSchedule {
private:
    GPUCommunicator gpuComm_;
public:
    void computeStep(/* flags */) override;
};
```

This abstraction may also be used to incorporate alternative force-calculation algorithms and their specific sequence of tasks. As a side benefit, it provides avenues to extend the MD library with new back-ends for force calculations.

```c++
class FMM_Schedule : public AbstractForceSchedule {
private:
    FMM_DomainTree tree_;
public:
    void computeStep(/* flags */) override;
};
```



## Schedule Builder

The task of selecting a particular schedule for a given node, or a given problem would be handled by this class. It encapsulates the selection and initialization of the implemented force schedules. It allows one to change schedules if required at runtime allowing more complex setups. The `selectSchedule` function essentially aggregates conditionality at the top level before a concrete schedule is initialized instead of doing in each step.

```c++
class ScheduleBuilder {
private:
    ScheduleType type_                     = None;
    bool          experimentalSchedules_ = false;
public:
    ScheduleBuilder(bool experimentalSchedules)
        : experimentalSchedules_(experimentalSchedules) {}

    void selectSchedule(/* input/execution context*/);
    void setScheduleType(ScheduleType type);

/*! Wraps around the protected init function in AbstractForceSchedule */
    AbstractForceSchedule* generateSchedule(/* app-specific args */);
};
```



## Decorated Schedules for Neighbour Search

In the current implementations, one conditionality that complicates maintainability of a force schedule is neighbor search (NS). As we see in the illustrated schedule, components of NS are interleaved with other tasks to maximize overlap. The order is different in the case where the MD step occurs without NS.

Separating these two flows would make the code more readable, but at the cost of code duplication. One way to mitigate this matter is by using decorated schedules. NS occurs once every 20-100 steps depending on the problem. At the initial cost of reduced overlap in an NS step, developers gain maintainability, reusability and modularity.

```c++
class NSSchedule : AbstractForceSchedule
{
public:
    NSSchedule(AbstractForceSchedule *noNS_schedule);

    void computeStep(/* flags */)
    {
        // Steps to calculate NS/pairlists

        inner_->computeStep(/* flags */);
    }

    void setInnerSchedule(AbstractForceSchedule *noNS_schedule);

protected:
    AbstractForceSchedule *inner_;
}
```

The added advantage is that the NS code is fully separated from the underlying force calculation routines. Irrespective of which concrete schedule is initialized (say forces, forces+virial, energy-only, etc), the NSSchedule object can potentially adapt at runtime accordingly.
