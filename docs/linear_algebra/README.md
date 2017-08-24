# Motivation

The use of the LAPACK and ScaLAPACK libraries for numerical linear algebra is ubiquitous in high-performance computing applications,
and many vendors and software providers use it as the basis of their own high-performance math libraries.
However, there currently exists no distributed GPU-enabled library.
Other packages for linear algebra with GPU support (or without GPU support but with better performances) exist,
but they have different interfaces, therefore major change in the application code to adopt them.

# Goals

We intend to develop a distributed linear algebra (DLA) interface with the following features:
- The DLA package for each call can be decided at runtime.
- An interface that works with ScaLAPACK matrices.
- Performs matrix layout conversion if needed.

# Potential Benefit

Each Application which uses ScaLAPACK will be able to use the DLA interface with minor changes in the code,
and can choose the best pperforming DLA package to speed-up the application.

# DLA Interface

The DLA interface features:
- Communicator utilities
- Matrix class (Not yet supported)
- DLA routines C++ wrappers (Not yet supported)
- DLA routines Fortran wrappers (Not yet supported)

# List of DLA packages

- ScaLAPACK (MKL, Libsci, Netlib, ...)
- ELPA
- D-Plasma (ParSEC)
- Chameleon (StarPU)

# List of supported routines

# List of routines with planned support

The routines which will be available are (including the ScaLAPACK corresponding name):
- Matrix-matrix multiplication (p\*gemm)
- Matrix-vector multiplication (p\*gemm)
- Cholesky factorization (p\*potrf)
- LU factorization (p\*getrf)
- Upper/lower triangular matrix inversion (p\*trtri)
- Matrix inversion (LU) (p\*getri)
- Eigenvalue solver (p\*{sy,he}ev{d,x}, p\*geev)
- Solution of linear equations system (p\*gesv)
