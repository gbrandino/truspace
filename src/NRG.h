/* 
 * Class that handles exact diagonalization, NRG and Sweep in a model-agnostic way ( as much as possible )
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l., Robert Konik - BNL
 * last revision: December 2016
 */

#ifndef NRG_H
#define NRG_H

#include "tensors.h"
#include "parameters.h"
#include "Hamiltonian.h"
#include "filenames.h"

extern "C" {
void dsyevr_(char *jobz, char *range, char *uplo, int *n,double *a, int *lda, double *vl, double *vu, int *il, int * iu, double *abstol, int *m, double *w,double *z, int *ldz, int *isuppz, double *work_for,int *lwork,int *iwork,int *liwork, int *info);
}



class NRG
{
    public:
        NRG();
        ~NRG();
        void run(parameters&, Hamiltonian&, filenames&);
        tensor2D<double> EigenVectors;
        int NIters;
        int NSweeps;


};
#endif
