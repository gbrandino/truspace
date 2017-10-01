/** @file worker.h 
 * \brief Contains driver Class for states, matrix element and NRG calculation. 
 *  It takes care of 
 *    - looping through the subspace found by parameters::analyze_subspace
 *    - finding chiral and non-chiral base states in each subspace
 *    - calling NRG routine on each subspace
 *    - calculating Matrix element on the eigenstates found, in each subspace. If the operator connects a set of subspaces, then those subspaces are diagonalized at the same time
 *
 * \author Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
 */

#ifndef WORKER_H
#define WORKER_H

#include "parameters.h"

/** Driver Class for states, matrix element and NRG calculation. 
 *  It takes care of 
 *    - looping through the subspace found by parameters::analyze_subspace
 *    - finding chiral and non-chiral base states in each subspace
 *    - calling NRG routine on each subspace
 *    - calculating Matrix element on the eigenstates found, in each subspace. If the operator connects a set of subspaces, then those subspaces are diagonalized at the same time
 *
 * \author Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
 */
class worker
{
    public:
        /** Default constructor 
        */
        worker();
        /** Default destructor 
        */
        ~worker();    
        /** 
         * @paramp[in] p paramerts object
         */    
        void run(parameters& p);

        // outfiles and dump files
        std::string EigenName;          /**< File name of eigenvalues file */
        std::string EigenVectName;      /**< File name of eigenvectors file */
        std::string MeOpName;           /**< File name of expectation value file */
//        std::string StoreName;          
        std::string OverlapName;        /**< File name of overlap file (Used when storing states or chiral matrix elements)*/
        std::string StateName;          /**< File name of states file (Used when storing states) */
        std::string MeChirName;         /**< File name of chiral matrix element file (Used when storing chiral matrix elements) */
};
#endif
