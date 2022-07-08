/** @file custom_parameters.h
 * \brief Input parameter class with model-specific parameters
 * Big fuzzy object that contains input paramenter
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date October 2019
 */

#ifndef CUSTOM_PARAMETERS_H
#define CUSTOM_PARAMETERS_H
#include "parameters.h"
/** Input parameter class. 
 * Big fuzzy object that contains input paramenter
*/
class custom_parameters : public parameters
{
    public:
        // input parameters
        int numPert;         /**< number of perturbations */
        int maxLev;          /**< Highest unperturbed state to include */
        int minLev;          /**< Lowest  unperturbed state to include */
        double* pConst;      /**< array containing the constant for the perturbation */
        int rSteps;          /**< number of values of the system size*/
        double rInc;         /**< increment of system size at each step */
        double rStart;       /**< initial value of the system size */
        bool TCSA;           /**< True for diagonalization, false to have just the size of the Hilbert space */
        bool saveEigenvec;   /**< True to save eigenvectors to file, false not to save them */
    
        // we hard code here the value of m and hbar
        const double m = 1;
        const double hbar = 1;

        // class methods
        /** Constructor
         * @param[in] modelfilename The name of the model input file
         * @param[in] nrgfilename The name of the nrg input file
         */
        custom_parameters(const char* modelfilename, const char* nrgfilename);
        /** Prints the input paramters, Verma modules size, Hilbert space size, subspace analysis and structure constant. <em> Needs to be called after analyse_subspaces() </em>
         */
        void print() const;
        /** Calculates Verma modules size, Hilbert space size, subspace analysis and structure constant.
         */
        /** Destructor
         */
        ~custom_parameters();
    private:
        /** Default constructor is private, to prevent instantiation of invalid object.
         */
        custom_parameters();
};

#endif    
