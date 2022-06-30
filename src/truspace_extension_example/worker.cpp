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
#include "worker.h"
#include "custom_parameters.h"
#include "harmonic_states.h"
#include "harmonic_filenames.h"
#include "tensors.h"
#include "harmonic_matrix_elements.h"
#ifdef MPI_PAR
#include "harmonic_Hamiltonian_MPI.h"
#include "NRG_MPI.h"
#else
#include "harmonic_Hamiltonian.h"
#include "NRG.h"
#endif
#include "custom_config.h"
#include <iostream>
#ifdef MPI_PAR
#include <mpi.h>
#endif

worker::worker()
{


}


worker::~worker()
{


}


void worker::run(custom_parameters& p)
{
#ifdef MPI_PAR
    int myrank_mpi,nprocs_mpi;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
#endif


#ifdef MPI_PAR 
    if (myrank_mpi == 0)
        {
#endif
        std::cout << std::endl;
        std::cout << "----------------"  << std::endl;
        std::cout << " Starting TCSA  "  << std::endl;
        std::cout << "----------------"  << std::endl;
        std::cout << std::endl;
#ifdef MPI_PAR 
        }
#endif

    harmonic_states st(p);
    harmonic_filenames fnames(p);

    // calculate states
    st.calculate(p);

    // calculate matrix elements
    harmonic_matrix_elements me(p,st);

    me.calculate(p,fnames,st);
        
    if ( p.calcExpVal )
        me.calculate_for_ExpVal(p, st);

#ifdef MPI_PAR 
    if (myrank_mpi == 0)
#endif
         std::cout << " Done with matrix elements calculation " << std::endl;
         harmonic_Hamiltonian H_iter(p, st, me);

         NRG nrg;
#ifdef MPI_PAR 
        if (myrank_mpi == 0)
            {
#endif
            std::cout << std::endl;
            std::cout << " --------------------------------------------- " << std::endl;
            std::cout << " Diagonalization of Hamiltonian " << std::endl;
#ifdef MPI_PAR 
            }
#endif

        for ( int r=0; r < p.rSteps; r++)
            {
            double R = p.rStart + r*p.rInc;
            H_iter.R = R;
#ifdef MPI_PAR 
            if (myrank_mpi == 0)
#endif                    
                std::cout << " Current system size : " << H_iter.R << std::endl;
                nrg.run(p, H_iter, fnames);
            }
            

}
