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
#include "states.h"
#include "filenames.h"
#include "chirME.h"
#ifdef MPI_PAR
#include "Hamiltonian_MPI.h"
#include "NRG_MPI.h"
#else
#include "Hamiltonian.h"
#include "NRG.h"
#endif
#include "config.h"
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


        states st(p,sub);
        filenames fnames(p,sub);
        if (p.READ_STATES)
            st.read_chiral(p, sub, fnames);
        else
            st.calculate_chiral(p,sub, fnames);

        chirME me(p,st,sub);
        me.calculate_or_read_chirME(p,fnames,st, sub);
        
        if ( p.calcExpVal )
            me.calculate_chirME_for_ExpVal(p, st, sub);
#ifdef MPI_PAR 
        if (myrank_mpi == 0)
#endif
        std::cout << " Done with chiral matrix elements calculation " << std::endl;
        // if all momentum subspaces needs to be diagonalized together do so, otherwise loop through momentum spaces
        if (p.allMomSubspace)
            {    
            // all momentum subspaces diagonalized together
            st.calculate_nonchiral(p,sub,-1);
            Hamiltonian H_iter(p, st, me, sub);
            H_iter.Momentum = -1;
            NRG nrg;
#ifdef MPI_PAR 
            if (myrank_mpi == 0)
                {
#endif
            std::cout << std::endl;
            std::cout << " ------------------------------ " << std::endl;
            std::cout << " Diagonalization of Hamiltonian " << std::endl;
#ifdef MPI_PAR 
                }
#endif

            for ( int r=0; r < p.rSteps; r++)
                {
                double R = p.rStart + r*p.rInc;
                H_iter.R = R;//R*2.*pi; //constant to match old results
#ifdef MPI_PAR 
                if (myrank_mpi == 0)
#endif
                    std::cout << " R value: " << r << " (" << H_iter.R << ")" << std::endl;  
                nrg.run(p, H_iter,fnames);  
                }
            }
        else
            {
            // loop through momentum subspaces
            for ( int mom=0; mom <= p.maxMomentum; mom++)
                {
                st.calculate_nonchiral(p,sub, mom);
                Hamiltonian H_iter(p, st, me, sub);
                H_iter.Momentum = mom;
                NRG nrg;
#ifdef MPI_PAR 
            if (myrank_mpi == 0)
                {
#endif
                std::cout << std::endl;
                std::cout << " --------------------------------------------- " << std::endl;
                std::cout << " Diagonalization of Hamiltonian for momentum " << mom  << std::endl;
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
                }
            }
            

        } //closing loop over subspaces
}
