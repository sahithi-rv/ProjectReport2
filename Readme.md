 ## Project Report 2 - Serial and Parallel CG for solving Discrete Poisson Equation
### File Description
#### Serial Codes directory - iter_serial.
*      cg.cpp - Pure CG serial code.
*     cg_pad.cpp - Modified CG code with padding described in section 3.1

#### Parallel Codes directory - iter_parallel.

*      cg_omp3.cpp - Main Parallel CG code.
*     cg_omp_next.cpp - modified CG with partial Matmul as discussed in section 3.2
*      cg_omp_tau.cpp - Added tau profiler clauses to obtain profiles of Parallel code.

#### Test Matrices directory - test_matrices
*   The file names have matrix sizes to indicate the file that needs to be passed as argument for a required matrix size.

* parallel_profile_output.txt - profile of parallel code for 4 threads
* serial_profile_output.txt - profile of serial code
* profile_serial_pad_output.txt - profile of modified padded serial code.



### Compile and Run Instructions
#### To compile and serial code to obtain profile
    
*         cd iter_serial
*         module load gcc/5.3.1
*          module load tau/2.27
*          tau_cxx.sh -std=c++11 -g -o tauperf_serial ./cg.cpp
*          tau_exec ./tauperf_serial  ../test_matrices/b1024 32
*          tau_exec ./tauperf_serial  ../test_matrices/b1600 40
*          tau_exec ./tauperf_serial  ../test_matrices/b10816 104
*          pprof
*          The arguments are Test Matrix file and the square root of matrix dimensions.
*     In order to simply compile the code run `make cg` and execute as `./cg <test_file_name>  <square root of matrix dimension>`. Note that the code should be compiled with -std=c++11 flag.
*      Similar instructions should be run for modified padded serial code - cg_pad.cpp
####      To compile and run parallel code and obtain the profiles
*          module purge
*          rm profile*
*          cd iter_parallel
*          module load gcc/5.3.1
*          module use /storage/work/a/awl5173/toShare/tauOMP/tau
*          module load adamsTauOMP_2.27
*          tau_cxx.sh -std=c++11 -fopenmp -g -o cg_omp_tau cg_omp_tau.cpp
*          ./cg_omp_tau ../test_matrices/b1024 32 4
*          ./cg_omp_tau ../test_matrices/b1600 40 4
*          ./cg_omp_tau ../test_matrices/b10816 104 4
*          pprof
*          The results in section 3 are shown for 4 thread-parallelism.  The arguments are the test matrix file, square root of matrix dimensions and the number of threads.
*      Similar instructions should be run for modified parallel code - cg_omp_next.cpp
*      In order to compile the code, run `make cg_omp3` and execute as `./cg_omp3 <test_file_name> <square root of matrix dimension> <number of threads>`.
    