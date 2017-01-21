/* 
 * Class that calculates and holds the Hamiltonian matrix and the function to calculates the elements needed for NRG and Sweep
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l., Robert Konik - BNL
 * last revision: December 2016
 */ 
#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "parameters.h"
#include "tensors.h"
#include "chirME.h"
#include "states.h"
#include <string>

class Hamiltonian
{
    public: 
        Hamiltonian(parameters&, states&, chirME&, int);
        ~Hamiltonian();
        void init();
        void print_eigenValues(int, int, double*, std::string, int );
        void print_eigenVectors(double*, std::string, int);
        void NRG_block_A(double*);
        void NRG_block_BC(double*, int, double*);
        void NRG_block_D(int);
        void sweep_block_BC(double*, int, double*, double*, int, int);

        
        tensor2D< double > H;
        int Msize;
        double R;
        int Momentum;
        parameters& p;
        states& st;
        chirME& me;
        int subspace;


};

#endif 
