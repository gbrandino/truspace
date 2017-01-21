/*
 * set of c-programs that diagonalize matrices based on code found in
 *  Numerical Recipes 
 * authors: Robert Konik - BNL
 * last revision: November 2016
 */
#ifndef DIAG_H
#define DIAG_H
#include "tensors.h"
class diag
{
    public:
    /* declaration of routines for the diagonalization of matrices of long doubles.
     * These routines are used in diagonalization of the Kac matrix where high
     * precision is required
     */
        void tred2_ld(tensor2D<long double >& a, int n,long double* d,long double* e,int ev_flag);
        void tqli_ld(long double* d,long double* e,int n,tensor2D< long double >& z, int ev_flag);
        void shell_ld(unsigned long n,long double d[],long double *r[]);

    /* declaration of routines for the diagonalization of matrices of doubles.
     * These routines are used in the diagonalization of the Hamiltonian
     * matrices
     */
        void tred2_double(double *a[], int n,double d[],double e[],int ev_flag);
        void tqli_double(double d[],double e[],int n,double *z[], int ev_flag);
        void shell_double(unsigned long n,double d[],double *r[]);
};
#endif
