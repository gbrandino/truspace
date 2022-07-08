/** @file truspace_harmonic.cpp
 * \brief Main program of Truspace harmonic oscillator extension (driver routine).
 *         Acts as a driver routine, reading input files and calling worker class to perform TCSA if required
 *
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date October 2017
 * \copyright GNU Public License.
 */

#include <iostream>
#include "custom_parameters.h"
#include "worker.h"
#ifdef MPI_PAR
#include <mpi.h>
#endif

int main(int argc, char** argv)
{
#ifdef MPI_PAR
    int myrank_mpi,nprocs_mpi;
    MPI_Init( &argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
#endif

#ifdef MPI_PAR
    if (argc<3)
        {
        if (myrank_mpi==0)
            std::cout<<" Usage: " << argv[0] << " [model input file name]  [nrg input file name]" << std::endl;
        MPI_Finalize();
        std::exit(0);
        }
    if ( nprocs_mpi != 1 and nprocs_mpi != 4 and nprocs_mpi != 9 and nprocs_mpi != 25 and nprocs_mpi != 36 and nprocs_mpi != 49 and nprocs_mpi != 64)
        {
        if (myrank_mpi==0)
            std::cout<<" Number of MPI processes needs to be a perfect square - quitting " << std::endl;
        MPI_Finalize();
        std::exit(0);
        }
#else
    if (argc<3)
        {
        std::cout<<" Usage: " << argv[0] << " [input file name] [nrg input file name]" << std::endl;
        std::exit(0);
        }
#endif

    // Create parameter object from the given input file
    custom_parameters par(argv[1], argv[2]);

    // Print recap of input and subspace analysis
#ifdef MPI_PAR
    if (myrank_mpi==0)
#endif
        par.print();   

    worker sys;
    // If TCSA is selectd, run it, otherwise, exit
    if (par.TCSA)
        {
        sys.run(par);
        }

#ifdef MPI_PAR
    MPI_Finalize();
#endif

    return 0;

}
