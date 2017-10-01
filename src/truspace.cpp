/** @file truspace.cpp
 * \brief Main program of Truspace (driver routine).
 *         Acts a a driver routine, reading input files and calling worker class to perform TCSA if required
 *
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date October 2017
 * \copyright GNU Public License.
 */

#include <iostream>
#include "parameters.h"
#include "worker.h"

int main(int argc, char** argv)
{
    if (argc<2)
        {
        std::cout<<" Usage: " << argv[0] << " [input file name] " << std::endl;
        std::exit(0);
        }
    // Create parameter object from the given input file
    parameters par(argv[1]);
    // Analyse subspace (also calculates structure constans, necessary for subspace analysis 
    par.analyse_subspaces(); 
    // Print recap of input, structure constants, subspace analysis,  and Verma modules dimensions
    par.print();   

    worker sys;
    // If TCSA is selectd, run it, otherwise, exit
    if (par.TCSA)
        {
        sys.run(par);
        }

    return 0;

}
