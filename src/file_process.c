#include "main.h"

void MatrixFileProcess(const char *path_file, CSR_Matrix **mat)
{
    *mat = (CSR_Matrix *)malloc(sizeof(CSR_Matrix));

    FILE *fp = NULL;
    if ((fp = fopen(path_file, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file! \"%s\"\n", path_file);
        exit(EXIT_FAILURE);
    }

    int n_row = 0, n_column = 0, nnz = 0;
    fscanf(fp, "%d%d%d", &n_row, &n_column, &nnz);
    (*mat)->n_row = n_row;
    (*mat)->n_column = n_column;
    (*mat)->nnz = nnz;

    if (((*mat)->row_idx = (int *)malloc(n_row * sizeof(int))) == NULL ||
        ((*mat)->row_ptr = (int *)malloc((n_row + 1) * sizeof(int))) == NULL ||
        ((*mat)->col_idx = (int *)malloc(nnz * sizeof(int))) == NULL ||
        ((*mat)->re_val = (double *)malloc(nnz * sizeof(double))) == NULL ||
        ((*mat)->im_val = (double *)malloc(nnz * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed! \"*mat component\"\n");
        exit(EXIT_FAILURE);
    }

    for (int index = 0; index < n_row; ++index)
    {
        fscanf(fp, "%d", (*mat)->row_idx + index);
    }
    for (int index = 0; index < n_row + 1; ++index)
    {
        fscanf(fp, "%d", (*mat)->row_ptr + index);
    }
    for (int index = 0; index < nnz; ++index)
    {
        fscanf(fp, "%d", (*mat)->col_idx + index);
    }
    for (int index = 0; index < nnz; ++index)
    {
        fscanf(fp, "%lf%lf", (*mat)->re_val + index, (*mat)->im_val + index);
    }

    fclose(fp);
}

void RHSFileProcess(const char *path_file, CSR_Vector **vec)
{
    *vec = (CSR_Vector *)malloc(sizeof(CSR_Vector));

    FILE *fp = NULL;
    if ((fp = fopen(path_file, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file! \"%s\"\n", path_file);
        exit(EXIT_FAILURE);
    }

    int n = 0;
    fscanf(fp, "%d", &n);
    (*vec)->n = n;

    if (((*vec)->row_idx = (int *)malloc(n * sizeof(int))) == NULL ||
        ((*vec)->re_val = (double *)malloc(n * sizeof(double))) == NULL ||
        ((*vec)->im_val = (double *)malloc(n * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed! \"*vec component\"\n");
        exit(EXIT_FAILURE);
    }

    for (int index = 0; index < n; ++index)
    {
        fscanf(fp, "%d%lf%lf", (*vec)->row_idx + index,
                (*vec)->re_val + index, (*vec)->im_val + index);
    }

    fclose(fp);
}