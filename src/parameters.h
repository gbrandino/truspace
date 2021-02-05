/** @file parameters.h
 * \brief Base input parameter class. 
 * Big fuzzy object that contains input parameters relavant for the NRG algorithm
 *
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date January 2021
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "tensors.h"
#include "sparse_tensor.h"
#include <complex>
#include <vector>
#include <string>
#include <map>

/** Input parameter class. 
 * Big fuzzy object that contains input parameters relavant for the NRG algorithm 
 *
*/
class parameters
{
    public:
        // input parameters
        int maxEigen;        /**< Maximum number of eigenvalues printed in the output file. */
        int Nsize;           /**< basic sizes the NRG matrix. */
        int Nstep;           /**< number of states added at each NRG iteration. The matrix diagonalized in NRG is Nsize+Nstep */
        bool NRG;            /**< True for NRG, false for full diagonalization */
        int nSweep;          /**< number of NRG sweep. Set 0 for basic NRG */
        bool saveEigenvec;   /**< True to save eigenvectors to file, false not to save them */
        bool calcExpVal;     /**< True to calculate expectation values on the eigenstates of the perturbed system of a primary operator, false not to */
        bool onlyDiagExpV;   /**< True to calculate only diagonal expectation values, false for the full matrix */



        // class methods
        /** Constructor
         * @param[in] The name of the nrg input file
         */
        parameters(const char*);
        
        /** Prints the input paramters, Verma modules size, Hilbert space size, subspace analysis and structure constant. <em> Needs to be called after analyse_subspaces() </em>
         */
        void print() const;

        /** Destructor
         */
        ~parameters();
    protected:

        /** Default constructor is private, to prevent instantiation of invalid object.
         */
        parameters();
};

#endif    
