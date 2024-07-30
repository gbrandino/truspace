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
#include "tensors.h"
#include "minimal_model_parameters.h"
#include "characters.h"
#include "minimal_model_states.h"
#include "minimal_model_filenames.h"
#include "chirME.h"
#ifdef MPI_PAR
#include "minimal_model_Hamiltonian_MPI.h"
#include "NRG_MPI.h"
#else
#include "minimal_model_Hamiltonian.h"
#include "NRG.h"
#endif
#include "minimal_model_config.h"
#include <complex>
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


void worker::run(minimal_model_parameters& p)
{
#ifdef MPI_PAR
    int myrank_mpi,nprocs_mpi;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
#endif


    // iterate through the subspaces 
    for (int sub=0; sub < p.numspace; sub++)
        {
        p.currentSubspace = sub;
#ifdef MPI_PAR 
        if (myrank_mpi == 0)
            {
#endif
            std::cout << std::endl;
            std::cout << "--------------------------------"  << std::endl;
            std::cout << " Now working on subspace "<< sub << " of " << p.numspace << std::endl;
            std::cout << "--------------------------------"  << std::endl;
            std::cout << std::endl;
#ifdef MPI_PAR 
            }
#endif


        minimal_model_states st(p);
        minimal_model_filenames fnames(p);
        if (p.READ_STATES)
            st.read_chiral(p, sub, fnames);
        else
            st.calculate_chiral(p, sub, fnames);

        chirME me(p,st);
        me.calculate_or_read_chirME(p,fnames,st);
        
        if ( p.calcExpVal )
            me.calculate_chirME_for_ExpVal(p, st);
#ifdef MPI_PAR 
        if (myrank_mpi == 0)
#endif
        std::cout << " Done with chiral matrix elements calculation " << std::endl;
        st.calculate_nonchiral(p, sub, p.momentum);
        minimal_model_Hamiltonian H_iter(p, st, me);
        H_iter.Momentum = p.momentum;
        NRG nrg;
#ifdef MPI_PAR 
        if (myrank_mpi == 0)
            {
#endif
            std::cout << std::endl;
            std::cout << " --------------------------------------------- " << std::endl;
            std::cout << " Diagonalization of Hamiltonian for momentum " << p.momentum  << std::endl;
#ifdef MPI_PAR 
            }
#endif

        for ( int r=0; r < p.rSteps; r++)
             {
             double R = p.rStart + r*p.rInc;
             H_iter.R = R;//R*2.*pi;  //constant to match old results
#ifdef MPI_PAR 
             if (myrank_mpi == 0)
#endif                    
                std::cout << " R value: " << r << " (" << H_iter.R << ")" << std::endl;
            nrg.run(p, H_iter, fnames);
            }
            

        } //closing loop over subspaces
}
