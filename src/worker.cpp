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

#include "worker.h"
#include "tensors.h"
#include "parameters.h"
#include "characters.h"
#include "states.h"
#include "filenames.h"
#include "chirME.h"
#include "Hamiltonian.h"
#include "NRG.h"
#include "config.h"
#include <complex>
#include <iostream>

worker::worker()
{


}


worker::~worker()
{


}


void worker::run(parameters& p)
{

    // checking if structure constants are real, not to have to implement additional hermitian matrix diagonalization
    int rflag=0;
    for (int i=0;i< p.maxp;i++)
        {
        for (int i1=0;i1<p.maxq;i1++)
            {
            for (int j=0;j< p.maxp;j++)
                {
                for (int j1=0;j1< p.maxq;j1++)
                    { 
                    for (int k=0;k< p.maxp;k++)
                        {
                        for (int k1=0;k1< p.maxq;k1++)
                            {
                            if (abs(p.cpt(i,i1,j,j1,k,k1) - std::conj(p.cpt(i,i1,j,j1,k,k1)))>0.0000000001)
                                { 
                                std::cout << std::endl;
                                std::cout<<" Warining - complex structure constant, the TCSA code can't handle that - quitting" << std::endl;
                                std::exit(1);
                                }
                            }
                        }
                    }
                }
            }
        }

    if (rflag==0)
        {
        std::cout << std::endl;
        std::cout << " Real structure constant, everything is fine :) " << std::endl;
        std::cout << std::endl;
        }


    // iterate through the subspaces 
    for (int sub=0; sub < p.numspace; sub++)
        {
        std::cout << std::endl;
        std::cout << "--------------------------------"  << std::endl;
        std::cout << " Now working on subspace "<< sub << " of " << p.numspace << std::endl;
        std::cout << "--------------------------------"  << std::endl;
        std::cout << std::endl;

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

        // if all momentum subspaces needs to be diagonalized together do so, otherwise loop through momentum spaces
        if (p.allMomSubspace)
            {    
            // all momentum subspaces diagonalized together
            st.calculate_nonchiral(p,sub,-1);
            Hamiltonian H_iter(p, st, me, sub);
            H_iter.Momentum = -1;
            NRG nrg;
            std::cout << std::endl;
            std::cout << " ------------------------------ " << std::endl;
            std::cout << " Diagonalization of Hamiltonian " << std::endl;
            for ( int r=0; r < p.rSteps; r++)
                {
                double R = p.rStart + r*p.rInc;
                H_iter.R = R;//R*2.*pi; //constant to match old results
                std::cout << " R value: " << r << " (" << H_iter.R << ")" << std::endl;  
                nrg.run(p, H_iter,fnames);  
                // exp val
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
                std::cout << std::endl;
                std::cout << " --------------------------------------------- " << std::endl;
                std::cout << " Diagonalization of Hamiltonian for momentum " << mom  << std::endl;
                for ( int r=0; r < p.rSteps; r++)
                    {
                    double R = p.rStart + r*p.rInc;
                    H_iter.R = R;//R*2.*pi;  //constant to match old results
                    
                    std::cout << " R value: " << r << " (" << H_iter.R << ")" << std::endl;
                    nrg.run(p, H_iter, fnames);
                    // exp val
                    }
                }
            }
            

        } //closing loop over subspaces
}
