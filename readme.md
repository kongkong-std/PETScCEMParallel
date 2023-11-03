# Userguide
## Objective
This software can solve <font color = blue>complex-type</font> time harmonic $Maxwell$ equation.
Linear system contains **coefficient matrix** and **right-hand side vector**, **preconditioner** as input.<br>
Before using it, there are some libraies needed to be installed.<br>
1. [PETSc](https://petsc.org/release/)
2. [AGMG](http://agmg.eu/)
3. [BLAS](https://netlib.org/blas/index.html)
4. [LAPACK](https://netlib.org/lapack/index.html)
5. [MUMPS](http://mumps-solver.org/index.php?page=home)
6. [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
7. [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)
8. [ScaLAPACK](https://netlib.org/scalapack/index.html)
9. [MPICH](https://www.mpich.org/)

## How to run it
First, compile source files and generate binary objective file.
1. cd exe
2. make

Second, attach path to matrix file, r.h.s. file, preconditioner file with command line parameter.
1. mpiexec -np
2. **modify** path_* in main.c
> ATTENTION: Matrix file must in <font color = red>csr</font> format!

Finally, setting solver parameter.
1. -ksp_type 
2. -ksp_max_it
3. -ksp_rtol
4. -ksp_monitor_true_residual
5. **-pc_type**
6. **-def_pc_agmg**
> choose 6 with AGMG solving preconditioning system