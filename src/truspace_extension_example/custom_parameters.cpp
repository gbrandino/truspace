/** Custom input parameter class. 
 *  Big fuzzy object that contains custom input paramenter
 *  \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 *  \date October 2019
 */

#include "custom_config.h"
#include "custom_parameters.h"
#include <iostream>
#include <fstream>
#include <vector>

custom_parameters::custom_parameters()
{
    

}


custom_parameters::~custom_parameters()
{
    delete [] pConst;
}


custom_parameters::custom_parameters(const char* model_input_file_name, const char* nrg_input_file_name)
:
parameters(nrg_input_file_name)
{
    std::ifstream input_file;
    int tmp;

    input_file.open(model_input_file_name);

    // reading input paramters from file. Only first elements of each line is read, the rest of the line is ignored and can be used for comments        
    numPert = 1; // in this example we hard-code the number of perturbation to 1
    input_file>>maxLev;  input_file.ignore(256, '\n');
    input_file>>minLev;  input_file.ignore(256, '\n');

    pConst=new double[numPert];
    for (int i=0;i<numPert;i++)
        {
        input_file>>pConst[i];  input_file.ignore(256, '\n');
        }

    input_file>>rSteps;  input_file.ignore(256, '\n');
    input_file>>rInc; input_file.ignore(256, '\n');   
    input_file>>rStart; input_file.ignore(256, '\n'); 


    input_file>>tmp; input_file.ignore(256, '\n');
    if (tmp == 0) 
        TCSA=false;
    else
        TCSA=true;

    input_file.close();


}
 

void custom_parameters::print() const
{
    // Printing input configuration, for easy review
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "  TruSpace - Perturbed free particle - " << VERSION        << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << " Version " << VERSION << std::endl;
    std::cout << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::cout << " Reporting input parameter" << std::endl;
    std::cout << std::endl;
    std::cout << " Number of perturbations ---> " << numPert << std::endl;
    std::cout << " Highest unperturbed state to include ---> " << maxLev << std::endl;
    std::cout << " Lowest unperturbed state to include ---> " << minLev << std::endl;


    if (maxEigen==0)
        { 
        std::cout << " Printing all available eigenvalues" << std::endl;
        }
    else
        {
        std::cout << " Maximum number of eigenstates in eigen file(s) ---> " << maxEigen << std::endl;
        }

    for (int i=0;i<numPert;i++)
        {
        std::cout << " Constant for perturbation "<< i+1 << " (omega) ---> " << pConst[i] << std::endl;;
        }

    std::cout << " Number of values of system size ---> " << rSteps << std::endl;
    std::cout << " Incremet at each step of the system size ---> " << rInc << std::endl;
    std::cout << " Starting value of the system size ---> " << rStart << std::endl;
    std::cout << " Number of states in the basic matrix for NRG ---> " << Nsize << std::endl;
    std::cout << " Number of states added in every NRG step ---> " << Nstep << std::endl;

    if (TCSA)
       {
       std::cout << " Full calculation including TCSA " << std::endl;
       }
    else 
        {
        std::cout <<" No TCSA "<< std::endl;
        }

    if (NRG)
        {
        std::cout << " NRG calcutation enabled" << std::endl;
        }
    else 
        {
        std::cout<<" No NRG calculation " << std::endl;
        }

    if (NRG and nSweep>0)
        {
        std::cout << " Performing " << nSweep << " NRG sweeps" << std::endl; 
        }


    if (calcExpVal)
        {
        std::cout << " Calculation of matrix elements of x^2.";
        }

    if(calcExpVal)
        {
        if (onlyDiagExpV)
            {
            std::cout << " Only diagonal matrix elements (<i|x^2|i>) will be calculated" << std::endl;
            }
        else
            {
            std::cout<<" Full matrix elements calculation" << std::endl;
            }
        }

    if (saveEigenvec)
        {
        std::cout << " Saving Eigenvector to file is enable" << std::endl;
        } 
    std::cout << std::endl;
    std::cout << "--  Perturbation of the CFT hamiltonian  = m omega^2 x^2/2" << std::endl;
    std::cout << std::endl;
}
    
