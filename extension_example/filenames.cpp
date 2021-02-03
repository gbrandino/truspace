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

}


filenames::filenames(parameters& p)
{

#ifdef MPI_PAR
    int myrank_mpi,nprocs_mpi;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
#endif


    std::stringstream temp;
    std::stringstream temp2;
    std::stringstream temp3;
    temp << "_lowestState_" << p.minLev << "_HighestState_" << p.maxLev;
    for (int i=0; i< p.numPert; i++)
        {
        temp << "_Per_" << i <<  "_omega_" << p.pConst[i]; 
        }

    temp2 << "Eigenvalues" << temp.str();
    EigenName = temp2.str();
    temp3 << "FINAL_" << temp2.str();
    EigenNameFinal = temp3.str();
    temp2.str(" ");    
   
    temp2 << "Eigenvectors" << temp.str();
    EigenVectName = temp2.str();
    temp2.str(" ");
    
    temp2 << "ExpVal_Op";
    temp2 << temp.str();
    MeOpName = temp2.str();
    temp2.str(" ");


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


    if (p.calcExpVal)
        std::cout << " Expectation values file: " << MeOpName << std::endl; 
    std::cout << std::endl;
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl << std::endl;
    std::cout << std::endl;
#ifdef MPI_PAR
    }
#endif    

}
