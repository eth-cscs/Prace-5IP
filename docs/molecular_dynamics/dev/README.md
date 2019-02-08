# Developer Facing API

# API Abstractions Overview

## Force Schedule API

The key idea is to have a collection of static force calculation schedules, that are hand-crafted for performance given the algorithm, execution context and required quantities. An abstract base class summarizes the minimum set of capabilities for a schedule.

```c++
class AbstractForceSchedule {
protected:
    SimulationState 	simuStat_;	// contains handles to system and state variables

    GMX_Communicator 	cr_;		// app-specific data structures for internal tasks
    ForceVector 		fr_;
    ForceVirial 		vir_;
    PerfCounter 		pc_;

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

The task of selecting a particular schedule for a given node, or a given problem would be handled by this class. It encapsulates the selection and initialization of the implemented force schedules. It allows one to change schedules if required at runtime allowing more complex setups.

```c++
class ScheduleBuilder {
private:
    AbstractForceSchedule root_;
public:
    void selectSchedule(/* input/execution context*/);
    void setScheduleType(scheduleType);
    void initSchedule(/* app-specific args */);  		// wraps around the protected init
    													// function in AbstractForceSchedule
    AbstractForceSchedule getSchedule;
};
```



## Multi-Level API

Make two aspects of the GMX API.

### Developer-Facing Abstraction (gmx-specific)

This part would setup the GROMACS-specific data structures, domain decomposition, communicator, load balancing, etc concerning the force-calculations to initialize an `IForceSchedule` object which would then be called by the high level API.

```c++
class NBICalculator {
private:
    ScheduleBuilder 				scheduleBuilder_; 	// build & init schedules
    std::shared_ptr<IForceSchedule> schedule_; 			// compute forces/energy/potential
    MDBookKeeping 					mdBook_; 			// maintain flags on DD, NS, etc

public:
	NBICalculator(/* args */) {
        /* init ScheduleBuilder
        --------------------------*/
        /*
        create GROMACS objects like DD, commrec,
        using the system data
        */
        ...
        scheduleBuilder_.selectSchedule(/* input/execution context*/);
        scheduleBuilder_.initSchedule(/* gmx-specific args */);
        schedule_ = scheduleBuilder_.getSchedule();
        mdBook.init(schedule_, system);
	}

    forceVector calculateForces() {
		return schedule_->computeStep(mdBook_.flags);
    }

    void update() {
        mdBook.sync(system); 						// updates if NS is needed, sets flags
        											// which quantities need to be computed
    }
};
```

The schedule API underneath would have it's own high performance implementation that overlaps computation, communication and reduction tasks. There would be customized schedules for architecture. algorithms and physics being simulated. This would be a part of GROMACS itself and need not be exposed outside as its implementation details are not so important.
