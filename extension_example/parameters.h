/** @file parameters.h
 * \brief Input parameter class. 
 * Big fuzzy object that contains input paramenter
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date October 2019
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

/** Input parameter class. 
 * Big fuzzy object that contains input paramenter
*/
class parameters
{
    public:
        // input parameters
        int numPert;         /**< number of perturbations */
        int maxLev;          /**< Highest unperturbed state to include */
        int minLev;          /**< Lowest  unperturbed state to include */
        int maxEigen;        /**< Maximum number of eigenvalues printed in the output file. */
        double* pConst;      /**< array containing the constant for the perturbation */
        int rSteps;          /**< number of values of the system size*/
        double rInc;         /**< increment of system size at each step */
        double rStart;       /**< initial value of the system size */
        int Nsize;           /**< basic sizes the NRG matrix. */
        int Nstep;           /**< number of states added at each NRG iteration. The matrix diagonalized in NRG is Nsize+Nstep */
        bool TCSA;           /**< True for diagonalization, false to have just Kac table structure constants and subspace analysis */
        bool NRG;            /**< True for NRG, false for full diagonalization */
        int nSweep;          /**< number of NRG sweep. Set 0 for basic NRG */
        bool calcExpVal;     /**< True to calculate expectation values on the eigenstates of the perturbed system of a primary operator, false not to */
        bool onlyDiagExpV;   /**< True to calculate only diagonal expectation values, false for the full matrix */
        bool saveEigenvec;   /**< True to save eigenvectors to file, false not to save them */
    
        // we hard code here the value of m and hbar
        const double m = 1;
        const double hbar = 1;

        // class methods
        /** Constructor
         * @param[in] The name of the input file
         */
        parameters(const char*);
        /** Prints the input paramters, Verma modules size, Hilbert space size, subspace analysis and structure constant. <em> Needs to be called after analyse_subspaces() </em>
         */
        void print() const;
        /** Calculates Verma modules size, Hilbert space size, subspace analysis and structure constant.
         */
        /** Destructor
         */
        ~parameters();
    private:
        /** Default constructor is private, to prevent instantiation of invalid object.
         */
        parameters();
};

#endif    
