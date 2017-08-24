# PRACE 5IP Work package 7.4

Task 7.4, namely "Provision of Numerical Libraries for Heterogeneous/Hybrid Architectures"  is part of Work Package 7 (WP7)
of the PRACE-5IP EU funded project. WP7 will provide applications enabling services for HPC applications codes that are important
for European academic and/or industrial researchers to ensure that these applications can effectively exploit current
and future PRACE systems.

## Specific Goals

Task 7.4 (T7.4) is focused on the productisation and adoption of low-level numerical libraries targeting heterogeneous/hybrid architectures for the Exascale era. These numerical libraries will provide the foundation upon which scientific applications and tools can be built. The two main goals are

* to contribute to the development of scalable and accelerated libraries; and
* to facilitate the adoption of such libraries in community codes.

This work is targeted directly at existent supercomputing systems in Europe with architectures that are real and genuine precursors to Exascale computing, and as such, is complementary to work already being undertaken by the European CoE.

The work described in this task will be delivered to the CoEs, and close collaboration with the CoEs is planned for the future. Two specific problems will be targeted initially (i) Linear algebra libraries for distributed hybrid architectures; and (ii) Non-bonded interactions (Coulomb and Lennard-Jones potentials) in classical molecular dynamics.

## Linear Algebra Libraries

The use of the LAPACK and ScaLAPACK libraries for numerical linear algebra is ubiquitous in high-performance computing applications,
and many vendors and software providers use it as the basis of their own high-performance math libraries.
However, there currently exists no distributed GPU-enabled library.
Other packages for linear algebra with GPU support exists but they have different interfaces,
therefore major changes in the application code are required to adopt them.

We intend to develop a distributed linear algebra (DLA) interface which allows to choose the library used for each call at runtime.
This allows to choose the best performing implementation without no further major changes in the application code.

For more information follow the link [Linear Algebra](linear_algebra)

## Non-bonded interactions in classical molecular dynamics

In classical molecular dynamics simulations – which consume a significant fraction of compute cycles at HPC centres across Europe – the calculation of the non-bonded interactions is by far the most expensive component. We will develop a common API and tool for non-bonded interactions that can be plugged into widely-used community codes such as NAMD, GROMACS, and so on.

Similar to comparable tools in linear algebra, such a library would implement different algorithms in the backend, such as Particle Mesh Ewald and Fast Multiple Methods.

For more information follow the link [Molecular Dynamics](molecular_dynamics)

# Contributing partners

+ [ETHZ CSCS](http://www.cscs.ch)
+ [CaSToRC](https://www.cyi.ac.cy)
+ [CINECA](http:www.cineca.it)
+ [NCSA](http://www.scc.acad.bg)
