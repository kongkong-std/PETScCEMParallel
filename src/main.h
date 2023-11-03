#ifndef MAIN_H_
#define MAIN_H_

// header file
#include <petscksp.h>

// struct definition
typedef struct csr_matrix
{
    /* data */
    int n_row;
    int n_column;
    int nnz;
    int *row_idx;   // size n x 1
    int *row_ptr;   // size (n + 1) x 1
    int *col_idx;   // size nnz x 1
    double *re_val; // size nnz x 1
    double *im_val; // size nnz x 1
} CSR_Matrix;

typedef struct csr_vector
{
    /* data */
    int n;
    int *row_idx;   // size n x 1
    double *re_val; // size n x 1
    double *im_val; // size n x 1
} CSR_Vector;

// function prototype
void MatrixFileProcess(const char *, CSR_Matrix **);
void RHSFileProcess(const char *, CSR_Vector **);

#endif