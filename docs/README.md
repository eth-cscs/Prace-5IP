# PRACE 5IP Work package 7.4

Task 7.4, namely "Provision of Numerical Libraries for Heterogeneous/Hybrid Architectures"  is part of Work Package 7 (WP7) 
of the PRACE-5IP EU funded project. WP7 will provide applications enabling services for HPC applications codes that are important 
for European academic and/or industrial researchers to ensure that these applications can effectively exploit current 
and future PRACE systems.

More specifically, Task 7.4 (T7.4), is focused on the productisation and adoption of low-level numerical libraries targeting heterogeneous/hybrid architectures for the Exascale era. These numerical libraries will provide the foundation upon which scientific applications and tools can be built. The two main goals are: to contribute to the development of scalable and accelerated libraries; and to facilitate the adoption of such libraries in community codes. This work is targeted directly at existent supercomputing systems in Europe with architectures that are real and genuine precursors to Exascale computing, and as such, is complementary to work already being undertaken by the European CoE. The work described in this task will be delivered to the CoEs, and close collaboration with the CoEs is planned for the future. Two specific problems will be targeted initially: (i) Linear algebra libraries for distributed hybrid architectures. The use of the LAPACK and ScaLAPACK libraries for numerical linear algebra is ubiquitous in high-performance computing applications, and many vendors and software providers use it as the basis of their own high-performance math libraries. However, there currently exists no distributed GPU-enabled library for dense linear algebra (cf. PBLAS/ScaLAPACK). We intend to systematically extend the MAGMA (Matrix Algebra on GPU and Multicore Architectures) library, in collaboration with its developers, to distributed architectures. Adoption of this library in community codes in materials science and electronic structure will be facilitated by providing access via a standard interface, hiding the underlying implementation from the user and thus avoiding the need for code modification. (ii) Non-bonded interactions (Coulomb and Lennard-Jones potentials) in classical molecular dynamics. In classical molecular dynamics simulations – which consume a significant fraction of compute cycles at HPC centres across Europe – the calculation of the non-bonded interactions is by far the most expensive component. We will develop a common API and tool for non-bonded interactions that can be plugged into widely-used community codes such as NAMD, GROMACS, and so on. Similar to comparable tools in linear algebra, such a library would implement different algorithms in the backend, such as Particle Mesh Ewald and Fast Multiple Methods.

The two main areas of interest are:

+ [Linear Algebra](linear_algebra)
+ [Molecular Dynamics](molecular_dynamics)

Contributing partners:
+ [http://www.cscs.ch](ETHZ CSCS)
+ CaSToRC
+ CINECA
+ NCSA
