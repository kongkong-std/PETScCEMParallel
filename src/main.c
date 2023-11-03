#include "main.h"
#include "agmg_pc.h"

int main(int argc, char **argv)
{
    PetscMPIInt irank, nproc;
    Mat helmholtz_a, helmholtz_pc;
    Vec helmholtz_rhs, helmholtz_sol;

    PetscFunctionBeginUser;

    PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));
    PetscCallMPI(MPI_Comm_rank(MPI_COMM_WORLD, &irank));
    PetscCallMPI(MPI_Comm_size(MPI_COMM_WORLD, &nproc));

#if 0
    char path_mat[] = "../input_file/helmholtz/helmholtz_mat";
    char path_pc[] = "../input_file/helmholtz/helmholtz_pc";
    char path_rhs[] = "../input_file/helmholtz/helmholtz_rhs";
#endif
    char path_mat[] = "../input_file/helmholtz_test/helmholtz/helmholtz_mat";
    char path_pc[] = "../input_file/helmholtz_test/helmholtz/helmholtz_pc";
    char path_rhs[] = "../input_file/helmholtz_test/helmholtz/helmholtz_rhs";
    char path_file[128];

    CSR_Matrix *mat_data = NULL, *pc_data = NULL;
    CSR_Vector *rhs_data = NULL;

    // mat file process
    sprintf(path_file, "%s_%08d", path_mat, irank);
    MatrixFileProcess(path_file, &mat_data);

    // pc file process
    sprintf(path_file, "%s_%08d", path_pc, irank);
    MatrixFileProcess(path_file, &pc_data);

    // rhs file process
    sprintf(path_file, "%s_%08d", path_rhs, irank);
    RHSFileProcess(path_file, &rhs_data);

    // create matrix
#if 1
    PetscCall(MatCreate(PETSC_COMM_WORLD, &helmholtz_a));
    //PetscCall(MatSetSizes(helmholtz_a, mat_data->n_row, PETSC_DECIDE,
    PetscCall(MatSetSizes(helmholtz_a, mat_data->n_row, mat_data->n_column,
                          PETSC_DETERMINE, PETSC_DETERMINE));
    PetscCall(MatSetFromOptions(helmholtz_a));
    PetscCall(MatSetUp(helmholtz_a));
    PetscCall(MatCreate(PETSC_COMM_WORLD, &helmholtz_pc));
    //PetscCall(MatSetSizes(helmholtz_pc, pc_data->n_row, PETSC_DECIDE,
    PetscCall(MatSetSizes(helmholtz_pc, pc_data->n_row, pc_data->n_column,
                          PETSC_DETERMINE, PETSC_DETERMINE));
    PetscCall(MatSetFromOptions(helmholtz_pc));
    PetscCall(MatSetUp(helmholtz_pc));

    // assigning values for matrix helmholtz_a
    printf(">>>> rank %d: mat_data->n_row = %d\n", irank, mat_data->n_row);
    for (int index = 0; index < mat_data->n_row; ++index)
    {
        int index_start = mat_data->row_ptr[index];
        int index_end = mat_data->row_ptr[index + 1];
        for (int index_j = index_start; index_j < index_end; ++index_j)
        {
            PetscScalar val_tmp = mat_data->re_val[index_j] + mat_data->im_val[index_j] * PETSC_i;
            PetscCall(MatSetValues(helmholtz_a, 1, mat_data->row_idx + index,
                                   1, mat_data->col_idx + index_j, &val_tmp, INSERT_VALUES));
        }
    }
    PetscCall(MatAssemblyBegin(helmholtz_a, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(helmholtz_a, MAT_FINAL_ASSEMBLY));

    // assigning values for matrix helmholtz_pc
    printf(">>>> rank %d: pc_data->nrow = %d\n", irank, pc_data->n_row);
    for (int index = 0; index < pc_data->n_row; ++index)
    {
        int index_start = pc_data->row_ptr[index];
        int index_end = pc_data->row_ptr[index + 1];
        for (int index_j = index_start; index_j < index_end; ++index_j)
        {
            PetscScalar val_tmp = pc_data->re_val[index_j] + pc_data->im_val[index_j] * PETSC_i;
            PetscCall(MatSetValues(helmholtz_pc, 1, pc_data->row_idx + index,
                                   1, pc_data->col_idx + index_j, &val_tmp, INSERT_VALUES));
        }
    }
    PetscCall(MatAssemblyBegin(helmholtz_pc, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(helmholtz_pc, MAT_FINAL_ASSEMBLY));
#endif

    // create vector
#if 1
    int n_local_col;
    PetscCall(MatGetLocalSize(helmholtz_a, NULL, &n_local_col));
    printf(">>>> rank %d: n_local_col = %d\n", irank, n_local_col);
    PetscCall(VecCreate(PETSC_COMM_WORLD, &helmholtz_rhs));
    //PetscCall(VecSetSizes(helmholtz_rhs, n_local_col, mat_data->n_column));
    PetscCall(VecSetSizes(helmholtz_rhs, n_local_col, PETSC_DETERMINE));
    PetscCall(VecSetFromOptions(helmholtz_rhs));
    PetscCall(VecSetUp(helmholtz_rhs));

    // assigning values for rhs
    printf(">>>> rank %d: rhs->n = %d\n", irank, rhs_data->n);
    for (int index = 0; index < rhs_data->n; ++index)
    {
        PetscScalar val_tmp = rhs_data->re_val[index] + rhs_data->im_val[index] * PETSC_i;
        PetscCall(VecSetValues(helmholtz_rhs, 1, rhs_data->row_idx + index, &val_tmp, INSERT_VALUES));
    }
    PetscCall(VecAssemblyBegin(helmholtz_rhs));
    PetscCall(VecAssemblyEnd(helmholtz_rhs));

    PetscCall(VecCreate(PETSC_COMM_WORLD, &helmholtz_sol));
    //PetscCall(VecSetSizes(helmholtz_sol, n_local_col, mat_data->n_column));
    PetscCall(VecSetSizes(helmholtz_sol, n_local_col, PETSC_DETERMINE));
    PetscCall(VecSetFromOptions(helmholtz_sol));
    PetscCall(VecSetUp(helmholtz_sol));
    for (int index = 0; index < rhs_data->n; ++index)
    {
        PetscScalar val_tmp = 0 + 0 * PETSC_i;
        PetscCall(VecSetValues(helmholtz_sol, 1, rhs_data->row_idx + index, &val_tmp, INSERT_VALUES));
    }
    PetscCall(VecAssemblyBegin(helmholtz_sol));
    PetscCall(VecAssemblyEnd(helmholtz_sol));
#endif

#if 0 // view matrix and vector
    PetscCall(MatView(helmholtz_a, PETSC_VIEWER_STDOUT_WORLD));
    puts("\n\n");
    // PetscCall(VecView(helmholtz_rhs, PETSC_VIEWER_STDOUT_WORLD));
#endif

#if 1
    // using ksp
    KSP ksp;
    PC pc;

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetOperators(ksp, helmholtz_a, helmholtz_pc));
    PetscCall(KSPGetPC(ksp, &pc));

    // user-defined preconditioner
    PetscBool def_pc_agmg = PETSC_FALSE;
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-def_pc_agmg", &def_pc_agmg, NULL));
    if (def_pc_agmg)
    {
        PetscCall(PCSetType(pc, PCSHELL));
        PetscCall(PCShellSetSetUp(pc, AGMGShellPCSetup));
        PetscCall(PCShellSetApply(pc, AGMGShellPCApply));
        PetscCall(PCShellSetContext(pc, NULL));
    }

    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(KSPSolve(ksp, helmholtz_rhs, helmholtz_sol));
#endif

    // free memory
    free(mat_data->row_idx);
    free(mat_data->row_ptr);
    free(mat_data->col_idx);
    free(mat_data->re_val);
    free(mat_data->im_val);
    free(mat_data);
    free(pc_data->row_idx);
    free(pc_data->row_ptr);
    free(pc_data->col_idx);
    free(pc_data->re_val);
    free(pc_data->im_val);
    free(pc_data);
    free(rhs_data->row_idx);
    free(rhs_data->re_val);
    free(rhs_data->im_val);
    free(rhs_data);

    PetscCall(MatDestroy(&helmholtz_a));
    PetscCall(MatDestroy(&helmholtz_pc));
    PetscCall(VecDestroy(&helmholtz_rhs));
    PetscCall(VecDestroy(&helmholtz_sol));
    PetscCall(KSPDestroy(&ksp));

    PetscCall(PetscFinalize());
    return 0;
}

/*
 * command:
 * mpirun -np 4 ./app_agmg
 * -ksp_monitor_true_residual -ksp_type fgmres -ksp_gmres_restart 50
 * -ksp_max_it 2000 -ksp_rtol 1e-8
 * (optional) -pc_type none
 * (optional) -def_pc_agmg
 * (optional) | tee result.log*
 */
