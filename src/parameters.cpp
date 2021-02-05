/** Input parameter class. 
 *  Big fuzzy object that contains 
 *    - input paramenter
 *    - Kac table 
 *    - structure constant 
 *    - verma module in each subspace
 *  \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 *  \date November 2017
 */

#include "config.h"
#include "parameters.h"
#include "tensors.h"
#include "cpt.h"
#include "characters.h"
#include "block.h"
#include "sparse_tensor.h"
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>

parameters::parameters()
{
    

}


parameters::~parameters()
{
}


parameters::parameters(const char* nrg_input_file_name)
{

    std::ifstream nrg_input_file;
    int tmp;
    nrg_input_file.open(nrg_input_file_name);

    nrg_input_file>>maxEigen;  nrg_input_file.ignore(256, '\n');
    nrg_input_file>>Nsize; nrg_input_file.ignore(256, '\n');
    nrg_input_file>>Nstep;  nrg_input_file.ignore(256, '\n');
    nrg_input_file>>tmp; nrg_input_file.ignore(256, '\n');
    if (tmp == 0)
        NRG=false;
    else
        NRG=true;

    nrg_input_file>>nSweep;  nrg_input_file.ignore(256, '\n');
    nrg_input_file>>saveEigenvec;  nrg_input_file.ignore(256, '\n');

    nrg_input_file>>calcExpVal;  nrg_input_file.ignore(256, '\n');
    nrg_input_file>>onlyDiagExpV;  nrg_input_file.ignore(256, '\n');

    nrg_input_file.close();    

    // done reading nrg input file. 
} 

