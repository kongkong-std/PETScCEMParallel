#ifndef AGMG_PC_H_
#define AGMG_PC_H_

#include <petscksp.h>

extern PetscErrorCode AGMGShellPCApply(PC, Vec, Vec);
extern PetscErrorCode AGMGShellPCSetup(PC);

// agmgpar interface function
#if 1
void zagmgpar_(int *n, double *a, int *ja,
               int *ia, double *f, double *x, int *ijob,
               int *iprint, int *nrest, int *iter, double *tol,
               int *MPI_COMM, int *listrank, int *ifirstlistrank);
#endif

/*
 * listrank and ifirstlistrank generation
 */
void ListRankGen(MPI_Comm, const int *, int, int, const int *,
                 int, int,
                 int **, int **, int *);

/*
 * comparing function
 *     if a > b, return 1
 *     if a < b, return -1
 *     if a == b, return 0
 */
int CompareNum(const void *, const void *);

/*
 * removing duplicate elements in array,
 * updating new size of array
 */
void RemoveDuplicateEle(int *, int *);

/*
 * find element index in array
 */
int FindIndexEle(const int *, int, int);

#endif