/* 
 * Class that takes care of 
 *  - loop through the subspace found by parameters::analyze_subspace
 *  - find chiral and non-chiral base states in each subspace
 *  - call NRG routine on each subspace
 *  - calculate Matrix element on the eigenstates found, in each subspace. If the operator connects a set of subspaces, then those subspaces are diagonalized at the same time
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l.
 * last revision: November 2016
 */

#ifndef WORKER_H
#define WORKER_H

#include "parameters.h"

class worker
{
    public:
        worker();
        ~worker();         
        void run(parameters&);

        // outfiles and dump files
        std::string EigenName;          // File name of eigenvalues file
        std::string EigenVectName;      // File name of eigenvectors file
        std::string MeOpName;           // File name
        std::string StoreName;          // File name
        std::string OverlapName;        // File name
        std::string StateName;          // File name of states file 
        std::string MeChirName;         // File name of mhiral matrix element file
};
#endif
