/*
 * Classes that contains and calculates chiral states by removing null vectors from each verma module
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l., Robert Konik - BNL
 * last revision: November 2016
 */

#ifndef STATES_H
#define STATES_H
#include "parameters.h"
#include "filenames.h"

class states
{
    public:
        states(parameters&, int );
        ~states();
    
        void calculate_chiral(parameters&, int, filenames& );        
        void calculate_nonchiral(parameters&, int, int );
        void read_chiral(parameters&, int, filenames& );
        

        int MaxChiralStates;
        int MaxChiralVerma;
        int MaxStates;
        int nstate;
        double* state_en;
        double* state_eff_en;
        int* state_left;
        int* state_right;
        int* state_vmod;
        long double* h; // array that contains the conformal dimension of the primary field labelling the verma module for each verma module in the subspace
        tensor3D< long double > cstate_ortho_modecoeff;
        tensor3D< long double > me_chi_0;
        tensor2D< long double > norm; 
        tensor3D< long double > cstate_mode; 
        tensor2D< long double > cstate_norm;
        int* num_chir; 
        tensor2D< int > cstate_lev; 
        tensor2D< int > cstate_nmode;
        tensor2D< int > cstate_ortho_lev; 
        tensor2D< int > cstate_ortho_nmode; 
        tensor3D< int > cstate_ortho_mode; 
        int MaxKac;
        tensor2D< int > nstates_lev;
        tensor2D< int > nstates_lev_each;
        tensor3D< int > state_right_lev; // Used only with sweep. Tensor that given the verma module, the level and the index inside a level return the global index of the states of the chiral level in the verma module, for the anti-chiral part
        tensor3D< int > state_left_lev; // Used only with sweep. Tensor that given the verma module, the level and the index inside a level return the global index of the states of the chiral level in the verma module, for the chiral part
        tensor4D< int > state_nonchiral_lev; // Used only with sweep. Tensor that given the verma module, the level of a non chiral states, and the index of its chiral and antichiral components gives the index on the non chiral state in the global indexing
        char overlapname[100];
        char statename[100];

    private:
        states();
};
#endif 
