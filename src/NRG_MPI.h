/** @file NRG.h
 * \brief Class that handles exact diagonalization, NRG and Sweep in a model-agnostic way ( as much as possible ) 
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
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


/** Class that handles exact diagonalization, NRG and Sweep in a model-agnostic way ( as much as possible )
 */
class NRG
{
    public:

        /** Constructor 
         */        
        NRG();

        /** Destructor
         */
        ~NRG();
       
        /** Routine that perfoms exact diagonalization, NRG or Sweep according to the input file
         * @param[in] p Reference to a parameter object
         * @param[in] H Reference to an Hamiltonian object
         * @param[in] fnames Reference to a filenames object
         */
        void run(parameters& p, Hamiltonian& H, filenames& fnames);
        tensor2D<double> EigenVectors;
        int NIters; /**< Number of NRG iteration to perform. Calculated by the code. */
        int NSweeps; /**< Number of Sweep to perform. Read from input file. */


};
#endif
