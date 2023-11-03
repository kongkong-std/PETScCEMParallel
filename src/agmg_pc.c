#include "agmg_pc.h"

PetscErrorCode AGMGShellPCSetup(PC pc)
{
    PetscFunctionBeginUser;

    Mat mat_pc;
    PetscCall(PCGetOperators(pc, NULL, &mat_pc));

    int irank, nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);

    int ijob = 1; // setup
    int iprint = 0, nrest = 50, iter = 2;
    double tol = 0.01;
    int mpi_comm = MPI_Comm_c2f(MPI_COMM_WORLD);
    int ifirstlistrank = 0;
    int *listrank = NULL, *new_col_a = NULL;

    int n_loc = 0, nnz_loc = 0;
    double *a = NULL, *f = NULL, *x = NULL;
    int *ja = NULL, *ia = NULL;

    // get local size
    // PetscCall(MatGetLocalSize(mat_pc, &n_loc, NULL));
    PetscBool done = PETSC_TRUE;
    const PetscInt *csr_ia = NULL;
    const PetscInt *csr_ja = NULL;
    int loc_row_start = 0, loc_row_end = 0;
    PetscCall(MatGetRowIJ(mat_pc, 0, PETSC_FALSE, PETSC_TRUE, &n_loc,
                          &csr_ia, &csr_ja, &done));
    PetscCall(MatGetOwnershipRange(mat_pc, &loc_row_start, &loc_row_end));

    nnz_loc = csr_ia[n_loc];
    int *loc_row_idx = NULL;
    /*
     * memory allocation
     *     size a = ( 2 x nnz_loc ) x 1, complex value
     *     size f = ( 2 x n_loc ) x 1, complex value
     *     size x = ( 2 x n_loc ) x 1, complex value
     *     size ja = nnz_loc x 1, int type
     *     size ia = ( n_loc + 1 ) x 1
     *     size loc_row_idx = n_loc x 1
     */
    if ((a = (double *)malloc(2 * nnz_loc * sizeof(double))) == NULL ||
        (f = (double *)malloc(2 * n_loc * sizeof(double))) == NULL ||
        (x = (double *)malloc(2 * n_loc * sizeof(double))) == NULL ||
        (ja = (int *)malloc(nnz_loc * sizeof(int))) == NULL ||
        (ia = (int *)malloc((n_loc + 1) * sizeof(int))) == NULL ||
        (loc_row_idx = (int *)malloc(n_loc * sizeof(int))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed! \"agmgpar_ parameters\"\n");
        exit(EXIT_FAILURE);
    }

    // assigning value
    /*
     * loc_row_idx = loc_row_restart : loc_row_end - 1
     * ia = csr_ia[0] : csr_ia[n_loc]
     * ja = csr_ja[0] : csr_ja[nnz_loc - 1]
     * a = csr format assigning value
     * f = vec_rhs[loc_row_start] : vec_rhs[loc_row_end - 1]
     */
    for (int index = loc_row_start; index < loc_row_end; ++index)
    {
        loc_row_idx[index - loc_row_start] = index;
    }
    for (int index = 0; index < n_loc + 1; ++index)
    {
        ia[index] = csr_ia[index];
    }
    for (int index = 0; index < nnz_loc; ++index)
    {
        ja[index] = csr_ja[index];
    }
    for (int index = 0; index < n_loc; ++index)
    {
        // matrix element
        int index_start = csr_ia[index];
        int index_end = csr_ia[index + 1];
        for (int index_j = index_start; index_j < index_end; ++index_j)
        {
            PetscScalar val_tmp;
            PetscCall(MatGetValue(mat_pc, loc_row_idx[index], csr_ja[index_j], &val_tmp));
            a[2 * index_j] = PetscRealPart(val_tmp);
            a[2 * index_j + 1] = PetscImaginaryPart(val_tmp);
        }

        // vector element
        f[2 * index] = 0;
        f[2 * index + 1] = 0;

        // x initial guess
        x[2 * index] = 0;
        x[2 * index + 1] = 0;
    }

    ListRankGen(MPI_COMM_WORLD, ja, nnz_loc, n_loc, ia,
                irank, nproc,
                &new_col_a, &listrank, &ifirstlistrank);

    // updating 0-base to 1-base
    for (int index = 0; index < n_loc + 1; ++index)
    {
        ia[index] += 1;
    }
    for (int index = 0; index < nnz_loc; ++index)
    {
        new_col_a[index] += 1;
    }

    zagmgpar_(&n_loc, a, new_col_a, ia, f, x,
              &ijob, &iprint, &nrest, &iter, &tol,
              &mpi_comm, listrank, &ifirstlistrank);

    // free memory
    free(a);
    free(f);
    free(x);
    free(ja);
    free(ia);
    free(loc_row_idx);
    free(listrank);
    free(new_col_a);

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode AGMGShellPCApply(PC pc, Vec vec_rhs, Vec vec_sol)
{
    PetscFunctionBeginUser;

    // test pc shell
    // PetscCall(VecCopy(x, y));

    Mat mat_pc;
    PetscCall(PCGetOperators(pc, NULL, &mat_pc));

    int irank, nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);

    int ijob = 2; // solve
    int iprint = 0, nrest = 50, iter = 2;
    double tol = 0.01;
    int mpi_comm = MPI_Comm_c2f(MPI_COMM_WORLD);
    int ifirstlistrank = 0;
    int *listrank = NULL, *new_col_a = NULL;

    int n_loc = 0, nnz_loc = 0;
    double *a = NULL, *f = NULL, *x = NULL;
    int *ja = NULL, *ia = NULL;

    // get local size
    // PetscCall(MatGetLocalSize(mat_pc, &n_loc, NULL));
    PetscBool done = PETSC_TRUE;
    const PetscInt *csr_ia = NULL;
    const PetscInt *csr_ja = NULL;
    int loc_row_start = 0, loc_row_end = 0;
    PetscCall(MatGetRowIJ(mat_pc, 0, PETSC_FALSE, PETSC_TRUE, &n_loc,
                          &csr_ia, &csr_ja, &done));
    PetscCall(MatGetOwnershipRange(mat_pc, &loc_row_start, &loc_row_end));
#if 0 // check local row
    printf(">>>> pc shell rank %d: n_loc = %d\n", irank, n_loc);
    printf(">>>> ~~~~ pc shell rank %d: loc_row_start = %d, loc_row_end = %d\n",
           irank, loc_row_start, loc_row_end );
#endif

    nnz_loc = csr_ia[n_loc];
    int *loc_row_idx = NULL;
    /*
     * memory allocation
     *     size a = ( 2 x nnz_loc ) x 1, complex value
     *     size f = ( 2 x n_loc ) x 1, complex value
     *     size x = ( 2 x n_loc ) x 1, complex value
     *     size ja = nnz_loc x 1, int type
     *     size ia = ( n_loc + 1 ) x 1
     *     size loc_row_idx = n_loc x 1
     */
    if ((a = (double *)malloc(2 * nnz_loc * sizeof(double))) == NULL ||
        (f = (double *)malloc(2 * n_loc * sizeof(double))) == NULL ||
        (x = (double *)malloc(2 * n_loc * sizeof(double))) == NULL ||
        (ja = (int *)malloc(nnz_loc * sizeof(int))) == NULL ||
        (ia = (int *)malloc((n_loc + 1) * sizeof(int))) == NULL ||
        (loc_row_idx = (int *)malloc(n_loc * sizeof(int))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed! \"agmgpar_ parameters\"\n");
        exit(EXIT_FAILURE);
    }

    // assigning value
    /*
     * loc_row_idx = loc_row_restart : loc_row_end - 1
     * ia = csr_ia[0] : csr_ia[n_loc]
     * ja = csr_ja[0] : csr_ja[nnz_loc - 1]
     * a = csr format assigning value
     * f = vec_rhs[loc_row_start] : vec_rhs[loc_row_end - 1]
     */
    for (int index = loc_row_start; index < loc_row_end; ++index)
    {
        loc_row_idx[index - loc_row_start] = index;
    }
    for (int index = 0; index < n_loc + 1; ++index)
    {
        ia[index] = csr_ia[index];
    }
    for (int index = 0; index < nnz_loc; ++index)
    {
        ja[index] = csr_ja[index];
    }
    for (int index = 0; index < n_loc; ++index)
    {
        // matrix element
        int index_start = csr_ia[index];
        int index_end = csr_ia[index + 1];
        for (int index_j = index_start; index_j < index_end; ++index_j)
        {
            PetscScalar val_tmp;
            PetscCall(MatGetValue(mat_pc, loc_row_idx[index], csr_ja[index_j], &val_tmp));
            a[2 * index_j] = PetscRealPart(val_tmp);
            a[2 * index_j + 1] = PetscImaginaryPart(val_tmp);
        }

        // vector element
        PetscScalar val_tmp;
        PetscCall(VecGetValues(vec_rhs, 1, loc_row_idx + index, &val_tmp));
        f[2 * index] = PetscRealPart(val_tmp);
        f[2 * index + 1] = PetscImaginaryPart(val_tmp);

        // x initial guess
        x[2 * index] = 0;
        x[2 * index + 1] = 0;
    }

    ListRankGen(MPI_COMM_WORLD, ja, nnz_loc, n_loc, ia,
                irank, nproc,
                &new_col_a, &listrank, &ifirstlistrank);

    // updating 0-base to 1-base
    for (int index = 0; index < n_loc + 1; ++index)
    {
        ia[index] += 1;
    }
    for (int index = 0; index < nnz_loc; ++index)
    {
        new_col_a[index] += 1;
    }

    zagmgpar_(&n_loc, a, new_col_a, ia, f, x,
              &ijob, &iprint, &nrest, &iter, &tol,
              &mpi_comm, listrank, &ifirstlistrank);

    // assigning value to vec_sol
    for (int index = 0; index < n_loc; ++index)
    {
        PetscScalar val_tmp = x[2 * index] + x[2 * index + 1] * PETSC_i;
        PetscCall(VecSetValues(vec_sol, 1, loc_row_idx + index, &val_tmp, INSERT_VALUES));
    }
    PetscCall(VecAssemblyBegin(vec_sol));
    PetscCall(VecAssemblyEnd(vec_sol));

    // free memory
    free(a);
    free(f);
    free(x);
    free(ja);
    free(ia);
    free(loc_row_idx);
    free(listrank);
    free(new_col_a);

    PetscFunctionReturn(PETSC_SUCCESS);
}

// agmgpar_ interface
int CompareNum(const void *a, const void *b)
{
    int num_1 = *(int *)a;
    int num_2 = *(int *)b;

    if (num_1 > num_2)
    {
        return 1;
    }
    else if (num_1 < num_2)
    {
        return -1;
    }
    else
    {
        return 0;
    }
}

void RemoveDuplicateEle(int *array, int *size)
{
    // if *size <= 1, keep original array and size
    if (*size > 1)
    {
        int index_i = 0, index_j = 0;
        for (index_i = 0, index_j = 0; index_j < *size; ++index_j)
        {
            if (array[index_i] != array[index_j])
            {
                ++index_i;
                array[index_i] = array[index_j];
            }
        }
        *size = index_i + 1;
    }
}

int FindIndexEle(const int *array, int size, int value)
{
    for (int index = 0; index < size; ++index)
    {
        if (array[index] == value)
        {
            return index;
        }
    }

    return -1; // if value not in array
}

void ListRankGen(MPI_Comm mpi_comm, const int *ja, int nnz_loc, int n_loc, const int *ia,
                 int irank, int nproc,
                 int **new_col_a, int **listrank, int *ifirstlistrank)
{
    int *col_idx = NULL;
    int size_col_idx = nnz_loc;
    if ((col_idx = (int *)malloc(size_col_idx * sizeof(int))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed! \"listrank col_idx\"\n");
        exit(EXIT_FAILURE);
    }
    for (int index = 0; index < size_col_idx; ++index)
    {
        col_idx[index] = ja[index];
    }

    // col_idx sorted in increasing order
    qsort(col_idx, size_col_idx, sizeof(int), CompareNum);

    // removing duplicates in array
    RemoveDuplicateEle(col_idx, &size_col_idx);

    int host_start_index = 0;
    MPI_Scan(&n_loc, &host_start_index, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    host_start_index -= n_loc;

    int host_start_list[nproc], n_loc_list[nproc];
    MPI_Allgather(&host_start_index, 1, MPI_INT, host_start_list, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&n_loc, 1, MPI_INT, n_loc_list, 1, MPI_INT, MPI_COMM_WORLD);

    *ifirstlistrank = 1;
    int *map = NULL;
    if ((map = (int *)malloc(size_col_idx * sizeof(int))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\"listrank map\"\n");
        exit(EXIT_FAILURE);
    }
    if ((*listrank = (int *)malloc(size_col_idx * sizeof(int))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\"listank listrank\"\n");
        exit(EXIT_FAILURE);
    }

    int local_variable_index = 0, nonlocal_variable_index = n_loc;
    for (int i = 0; i < size_col_idx; ++i)
    {
        int min_index = host_start_list[irank];
        int max_index = host_start_list[irank] + n_loc - 1;
        if (col_idx[i] >= min_index && col_idx[i] <= max_index)
        {
            map[i] = local_variable_index;
            local_variable_index++;
            (*listrank)[map[i]] = -1;
        }
        else
        {
            map[i] = nonlocal_variable_index;
            nonlocal_variable_index++;

            int host_rank = 0;
            for (int rank = 0; rank < nproc; ++rank)
            {
                int min_id = host_start_list[rank];
                int max_id = host_start_list[rank] + n_loc_list[rank] - 1;
                if (col_idx[i] >= min_id && col_idx[i] <= max_id)
                {
                    host_rank = rank;
                    break;
                }
            }
            (*listrank)[map[i]] = host_rank;
        }
    }

    if ((*new_col_a = (int *)malloc(nnz_loc * sizeof(int))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\"listrank new_col_a\"\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < nnz_loc; ++i)
    {
        int index = FindIndexEle(col_idx, size_col_idx, ja[i]);
        (*new_col_a)[i] = map[index];
    }

    // free memory
    free(map);
    free(col_idx);
}
