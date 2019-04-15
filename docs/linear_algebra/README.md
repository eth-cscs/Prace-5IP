# Motivation

The use of the LAPACK and ScaLAPACK libraries for numerical linear algebra is ubiquitous in high-performance computing applications,
and many vendors and software providers use it as the basis of their own high-performance math libraries.
However, there currently exists no distributed GPU-enabled library.
Other packages for linear algebra with GPU support (or without GPU support but with better performances) exist,
but they have different interfaces, therefore major changes in the application code are needed to adopt them.

# Goals

The distributed linear algebra interface (DLA interface) has been developed with the following goals:
- Possibility to choose at runtime which computation library to use.
- C++ objects to simplify the development of new applications which provide distributed matrix functionalities.
- Interoperability with ScaLAPACK. Since the DLA interface will support only a small subset of the routines implemented in ScaLAPACK (the most relevant for scientific applications) the other routines have to be accessed using ScaLAPACK directly.
- Minimal changes to existing applications: Replacing ScaLAPACK with the DLA interface should require minimal changes to the source code of existing applications written in Fortran, C or C++.
- Possibility for adding support for new libraries without API changes.

# Potential Benefit

Each Application which uses ScaLAPACK will be able to use the DLA interface with minor changes in the code,
and can benefit from the use of the best performing DLA package to increase the performance of the application.

# DLA Interface

The source code of the library is available [here](https://github.com/eth-cscs/DLA-interface).
The DLA interface features:
- Communicator utilities
- Matrix class
- DLA routines C++ wrappers
- DLA routines Fortran wrappers

## List of DLA packages

The DLA library supported are:
- ScaLAPACK (MKL, Libsci, Netlib, ...)
- ELPA
- D-Plasma (ParSEC)

## List of routines with planned support

The routines which will be available are (including the ScaLAPACK corresponding name):
- Matrix-matrix multiplication (p\*gemm)
- Cholesky factorization (p\*potrf)
- LU factorization (p\*getrf)
- Upper/lower triangular matrix inversion (p\*trtri)
- Matrix inversion (Cholesky) (p\*potri)
- Eigenvalue solver (p\*{sy,he}ev{d,x})

# Interfaces and examples

## C++ interface
A full documentation is available in the [source code](https://github.com/eth-cscs/DLA-interface)

1. Initialization and finalization:
   ```cpp
     void CommunicatorManager::initialize(bool initialize_mpi = true);
     void CommunicatorManager::initialize(int nr_cores, bool initialize_mpi = true);
     void CommunicatorManager::initialize(int nr_cores, int* argc, char*** argv, bool initialize_mpi = true);

     void CommunicatorManager::finalize();
   ```

2. Communicator management:

   - Creations:
     ```cpp
       Communicator2DGrid& CommunicatorManager::createCommunicator2DGrid(MPI_Comm base_comm, int nr_rows, int nr_cols, Ordering comm_ordering);
     ```

   - Destruction:
     ```cpp
       void CommunicatorManager::free2DGrid(Communicator2DGrid& grid);
       void CommunicatorManager::free2DGridFromMPIComm(MPI_Comm comm);
       void CommunicatorManager::free2DGridFromBlacsContext(BlacsContextType id);
     ```

3. Linear algebra calls:
   - Matrix-matrix multiplication:
     ```cpp
        template <class ElType>
        void matrixMultiplication(OpTrans trans_a, OpTrans trans_b, ElType alpha, const DistributedMatrix<ElType>& mat_a, const DistributedMatrix<ElType>& mat_b, ElType beta, DistributedMatrix<ElType>& mat_c, SolverType solver, int print_timers = 0);
     ```
   - Cholesky factorization:
     ```cpp
        template <class ElType>
        void choleskyFactorization(UpLo uplo, DistributedMatrix<ElType>& mat, SolverType solver, int print_timers = 0);
     ```
   - LU factorization:
     ```cpp
        template <class ElType>
        void LUFactorization(DistributedMatrix<ElType>& mat, std::vector<int>& ipiv, SolverType solver, int print_timers = 0);

        template <class ElType>
        void LUFactorization(DistributedMatrix<ElType>& mat, int* ipiv, SolverType solver, int print_timers = 0);
     ```
   - Inverse (triangular)
     ```cpp
        template <class ElType>
        void triangularInverse(UpLo uplo, Diag diag, DistributedMatrix<ElType>& mat, SolverType solver, int print_timers = 0);
     ```
   - Inverse (Cholesky)
     ```cpp
        template <class ElType>
        void choleskyInverse(UpLo uplo, DistributedMatrix<ElType>& mat, SolverType solver, int print_timers = 0);
     ```
   - Eigenvalue solver:
     ```cpp
        template <class ElType>
        void hermitianEigenvectors(UpLo uplo, DistributedMatrix<ElType>& mat, std::vector<BaseType<ElType>>& evalues, DistributedMatrix<ElType>& evectors, SolverType solver, int print_timers = 0);

        template <class ElType>
        void hermitianEigenvectors(UpLo uplo, DistributedMatrix<ElType>& mat, BaseType<ElType>* evalues, DistributedMatrix<ElType>& evectors, SolverType solver, int print_timers = 0);
     ```

### Example
```cpp
#include "distributed_matrix.h"
#include "dla_interface.h"

using namespace dla_interface;

double element(int i, int j);
// Returns the (i, j)-th element of the matrix.

int main(int argc, char** argv) {

  int n = 10240; // matrix size
  int nb = 128;

  // Use a (3x2) rank grid 
  int p = 3;
  int q = 2;

  int nr_threads = 4; // number of threads per rank.

  const UpLo uplo = Lower;

  SolverType solver = ScaLAPACK; // Other solvers: DPlasma, ELPA, ...
  DistributionType dist = scalapack_dist; // Other distibution tile_dist

  // Initialize DLA interface
  comm::CommunicatorManager::initialize(nr_threads, &argc, &argv, true);

  // Set up the communicators and Blacs context
  auto& comm_grid = comm::CommunicatorManager::createCommunicator2DGrid(MPI_COMM_WORLD, p, q, RowMajor);

  // Set up the matrix
  DistributedMatrix<double> mat(n, n, nb, nb, comm_grid, dist);

  // Set the elements.
  for (int j = 0; j < mat.localSize().second; ++j) {
    for (int i = 0; i < mat.localSize().first; ++i) {
      Local2DIndex local_index(i, j);
      Global2DIndex global_index = mat.getGlobal2DIndex(local_index);
      mat(local_index) = element(global_index.row, global_index.col);

  choleskyFactorization(uplo, mat, solver);

  return 0;
}
```

## How to modify an existing Fortran application
> Note:
> In this sections the line that have to be replaced are prefixed with `-`, while the lines to be added with a `+`.
> E.g.
> ```
> -   <Line to be removed>
> +   <Line to be added>
> ```


1. Add fortran module
  ```fortran
  USE dla_interface
  ```

2. Initialize and finalize the DLA interface:

   - Initialization with internal `MPI_Init`
     ```fortran
     +   CALL dlai_initialize(nr_threads_per_rank, 1)
     ```
     or initialization where MPI has already been initialized by the application
     ```fortran
     +   CALL dlai_initialize(nr_threads_per_rank, 0)
     ```

   - Finalization
     ```fortran
     +   CALL dlai_finalize()
     ```

3. Replace the creation of the blacs context:

   ```fortran
   -     CALL blacs_gridinit(context, order, nprow, npcol)
   +     context = dlai_create_2d_grid(mpi_comm, nprow, npcol, order)
   ```

4. Replace linear algebra calls:
   - Matrix-matrix multiplication:
     ```fortran
     -   CALL pdgemm(trans_a, trans_b, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc, info)
     +   CALL dlai_d_matrix_multiplication(trans_a, trans_b, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc, solver, info)
     ```
   - Cholesky factorization:
     ```fortran
     -   CALL pdpotrf(uplo, n, a, ia, ja, desca, info)
     +   CALL dlai_d_cholesky_factorization(uplo, n, a, ia, ja, desca, solver, info)
     ```
   - LU factorization:
     ```fortran
     -   CALL pdgetrf(m, n, a, ia, ja, desca, ipiv, info)
     +   CALL dlai_d_lu_factorization(m, n, a, ia, ja, desca, ipiv, solver, info)
     ```
   - Inverse (triangular)
     ```fortran
     -   CALL pdtrtri(uplo, diag, n, a, ia, ja, desca, info)
     +   CALL dlai_d_triangular_inverse(uplo, diag, n, a, ia, ja, desca, solver, info)
     ```
   - Inverse (Cholesky)
     ```fortran
     -   CALL pdpotri(uplo, n, a, ia, ja, desca, info)
     +   CALL dlai_d_cholesky_inverse(uplo, n, a, ia, ja, desca, solver, info)
     ```
   - Eigenvalue solver:
     ```fortran
     -   CALL pdsyevd("V", uplo, n, a, ia, ja, desca, evals, v, iv, jv, descv, iwork, lwork, iwork, liwork, info)
     +   CALL dlai_d_hermitian_eigenvectors(uplo, n, a, ia, ja, desca, evals, v, iv, jv, descv, solver, info)
     ```
