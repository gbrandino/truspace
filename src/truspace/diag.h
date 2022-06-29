/** @file diag.h
 * \brief Contains the declaration of the diag Class, which provides a set of c-programs that diagonalize matrices based on code found in Numerical Recipes, for both double and long double. 
 * \author Robert Konik - BNL, Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2016
 */

#ifndef DIAG_H
#define DIAG_H
#include "tensors.h"
/** Class containg a set of c-programs that diagonalize matrices based on code found in Numerical Recipes, for both double and long double. 
 */
class diag
{
    public:
    /* declaration of routines for the diagonalization of matrices of long doubles.
     * These routines are used in diagonalization of the Kac matrix where high
     * precision is required
     */
        /** Routine reduces a real symmetric matrix to a symmetric tridiagonal matrix using and accumulating orthogonal similarity transformations. Long double version.
          * @param[in]  a Reference of tensor2D<long double > of size n x n. The matrix to be reduced
          * @param[in]  n Size of the matrix 
          * @param[out] d Array of size n that will contain the diagonal element of the tridiagonal matrix
          * @param[out] e Array of size n that will contain the subdiagonal elements of the tridiagonal matrix in its last n-1 positions
          * @param[in]  ev_flag Set to 1 if eigenvectors are required, zero otherwise
          */
        void tred2_ld(tensor2D<long double >& a, int n,long double* d,long double* e,int ev_flag);

        /** Routine educes a real symmetric matrix to a symmetric tridiagonal matrix using and accumulating orthogonal similarity transformations. Double version.
          * @param[in/out] d Array of size n that contains the diagonal element of the tridiagonal matrix. On output, it contains the eigenvalues.
          * @param[in] e Array of size n that contains the subdiagonal elements of the tridiagonal matrix in its last n-1 positions
          * @param[in]  n Size of the matrix 
          * @param[out] e Reference of tensor2D<long double > of size n x n. It will contain the eigenvectors is ev_flag is set to 1.  
          * @param[in]  ev_flag Set to 1 if eigenvectors are required, zero otherwise
          */
        void tqli_ld(long double* d,long double* e,int n,tensor2D< long double >& z, int ev_flag);

        /** Routine to sort an array in place, also providing sorting indices. Long double version. 
          * @param[in]     n The size of the array to sort 
          * @param[in/out] d The array to sort. On output, it contains the sorted array
          * @param[out]    r Matrix of size n x n+1 contaning the sorting indeces  
          */
        void shell_ld(unsigned long n,long double d[],long double *r[]);

    /* declaration of routines for the diagonalization of matrices of doubles.
     * These routines are used in the diagonalization of the Hamiltonian
     * matrices
     */

        /** Routine reduces a real symmetric matrix to a symmetric tridiagonal matrix using and accumulating orthogonal similarity transformations. Long double version.
          * @param[in]  a Reference of tensor2D<long double > of size n x n. The matrix to be reduced
          * @param[in]  n Size of the matrix 
          * @param[out] d Array of size n that will contain the diagonal element of the tridiagonal matrix
          * @param[out] e Array of size n that will contain the subdiagonal elements of the tridiagonal matrix in its last n-1 positions
          * @param[in]  ev_flag Set to 1 if eigenvectors are required, zero otherwise
          */
        void tred2_double(double *a[], int n,double d[],double e[],int ev_flag);

        /** Routine educes a real symmetric matrix to a symmetric tridiagonal matrix using and accumulating orthogonal similarity transformations. Double version.
          * @param[in/out] d Array of size n that contains the diagonal element of the tridiagonal matrix. On output, it contains the eigenvalues.
          * @param[in] e Array of size n that contains the subdiagonal elements of the tridiagonal matrix in its last n-1 positions
          * @param[in]  n Size of the matrix 
          * @param[out] e Reference of tensor2D<long double > of size n x n. It will contain the eigenvectors is ev_flag is set to 1.  
          * @param[in]  ev_flag Set to 1 if eigenvectors are required, zero otherwise
          */
        void tqli_double(double d[],double e[],int n,double *z[], int ev_flag);

        /** Routine to sort an array in place, also providing sorting indices. Long double version. 
          * @param[in]     n The size of the array to sort 
          * @param[in/out] d The array to sort. On output, it contains the sorted array
          * @param[out]    r Matrix of size n x n+1 contaning the sorting indeces  
          */
        void shell_double(unsigned long n,double d[],double *r[]);
};
#endif
