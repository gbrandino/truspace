/*
 * Class that contains and calculates states
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l. , Robert Konik - BNL
 * last revision: October 2019
 */

#include "states.h"
#include "parameters.h"
#include "config.h"
#include <iostream>
#ifdef MPI_PAR
#include <mpi.h>
#endif


states::states()
{


}

states::states(parameters& p)
{
#ifdef MPI_PAR
    int myrank_mpi,nprocs_mpi;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
#endif

    nstate = (p.maxLev - p.minLev + 1)* 2 ;
    state_en=new double[nstate];
    state_quantum_number=new double[nstate];

    for (int j=0;j<nstate; j++)
        {
        state_en[j]=0.;
        state_quantum_number[j]=0;
        }
}

states::~states()
{
    delete [] state_en;
    delete [] state_quantum_number;
}

void states::calculate(parameters& p)
{
    // MaxPart and MaxNumGrade are (at least for now) global constant 
    // defined in config.h

#ifdef MPI_PAR
    int myrank_mpi,nprocs_mpi;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
#endif
       
#ifdef MPI_PAR
    if (myrank_mpi == 0)
#endif  
        std::cout << " Number of unpertubed states "<< nstate << std::endl;
    
    for (int i = p.minLev; i <= p.maxLev; i++)
    {
        state_quantum_number[ (i - p.minLev)* 2 ]     = ( 2 * i - 1) * PI;
        state_quantum_number[ (i - p.minLev)* 2 + 1 ] = ( 2 * i    ) * PI;
        double en =  p.hbar * p.hbar / 2 / p.m * (PI*(2*i-1))*(PI*(2*i-1)); // the size factor L^2 is put later, since we are looping on it
        double en2  = p.hbar * p.hbar / 2 / p.m * (2*PI*i)*(2*PI*i);
        state_en[ (i - p.minLev)* 2 ] =  en;
        state_en[ (i - p.minLev)* 2 + 1 ] = en2;
    }

}




