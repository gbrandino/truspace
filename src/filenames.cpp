/** Output filenames helper class. 
 *   - Contains the filenames helper class 
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
 */
#include "filenames.h"
#include "parameters.h"
#include <string>
#include <sstream>
#include <iostream>
#ifdef MPI_PAR
#include <mpi.h>
#endif

filenames::filenames()
{


}


filenames::~filenames()
{

    delete [] MeChirName;

}


filenames::filenames(parameters& p, int subspace)
{

#ifdef MPI_PAR
    int myrank_mpi,nprocs_mpi;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
#endif


    std::stringstream temp;
    std::stringstream temp2;
    std::stringstream temp3;
    temp << "_p_" << p.p << "_q_" << p.q << "_minLev_" << p.minLev << "_maxLev_" << p.maxLev;
    for (int i=0; i< p.numPert; i++)
        {
        temp << "_Per_" << p.r[i] << "_" << p.s[i] << "_c_" << p.pConst[i]; 
        }
    temp << "_Subspace_" << subspace;

    temp2 << "Eigenvalues" << temp.str();
    EigenName = temp2.str();
    temp3 << "FINAL_" << temp2.str();
    EigenNameFinal = temp3.str();
    temp2.str(" ");    
   
    temp2 << "Eigenvectors" << temp.str();
    EigenVectName = temp2.str();
    temp2.str(" ");
    
    temp2 << "ExpVal_Op";
    for (int i=0; i<p.nModes_l_op; i++)
        temp2 << "_L_" << p.mode_l_op[i];
    temp2 <<  "_phi_" << p.rl_op << "_" << p.sl_op;
    for (int i=0; i<p.nModes_r_op; i++)
        temp2 << "_{bar}L_" << p.mode_r_op[i];
    temp2 <<  "_{bar}phi_" << p.rr_op << "_" << p.sr_op;
    temp2 << temp.str();
    MeOpName = temp2.str();
    temp2.str(" ");

    // resetting temp for states, overlaps and chiral matrix elements
    temp.str(" ");
    temp << "_p_" << p.p << "_q_" << p.q << "_minLev_" << p.minLev << "_maxLev_" << p.maxLev;

    temp2 << "States" << temp.str();
    StateName = temp2.str();
    temp2.str(" ");
    
    temp2 << "Overlaps" << temp.str();
    OverlapName  = temp2.str();
    temp2.str(" ");
    
    // resetting temp for states, overlaps and chiral matrix elements
    MeChirName = new std::string[p.numPert];
    temp.str(" ");
    temp << "_p_" << p.p << "_q_" << p.q << "_minLev_" << p.minLev << "_maxLev_" << p.maxLev;

    for (int i=0; i < p.numPert; i++)
        {
        temp2 << "ME_chir" << temp.str();
        temp2 << "_Op_phi_" << p.r[i] << "_" << p.s[i];
        MeChirName[i]= temp2.str();
        }

#ifdef MPI_PAR
    if ( myrank_mpi == 0)
    {
#endif
    std::cout << std::endl;    
    std::cout << " Output and dump file names prefixes " << std::endl;
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl << std::endl;
    std::cout << " Eigenvalue file name: " << EigenName << std::endl;
    std::cout << " Eigenvectors file name: " << EigenVectName << std::endl;
    if (p.NRG)
        {
        std::cout << " Eigenvalue file name for last NRG/Sweep iteration: " << EigenNameFinal << std::endl;
        }

    if (p.SAVE_STATES)
        {
        std::cout << " Saving chiral states to: " << StateName << std::endl;
        std::cout << " Saving overlaps ( to determine non-null vectors) to: " << OverlapName << std::endl;
        }
    if (p.READ_STATES)
        {
        std::cout << " Reading chiral states from: " << StateName << std::endl; 
        std::cout << " Reading overlaps ( to determine non-null vectors) from: " << OverlapName << std::endl;
        }

    if(p.SAVE_ME)
        {
        for (int i=0; i < p.numPert; i++)
            {
            std::cout << " Saving chiral matrix elements to: " << MeChirName[i] << std::endl;
            }
        }

    if(p.READ_ME)
        {
        for (int i=0; i < p.numPert; i++)
            {
            std::cout << " Reading chiral matrix elements from: " << MeChirName[i] << std::endl;
            }
        }


    if (p.calcExpVal)
        std::cout << " Expectation values file: " << MeOpName << std::endl; 
    std::cout << std::endl;
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl << std::endl;
    std::cout << std::endl;
#ifdef MPI_PAR
    }
#endif    

}
