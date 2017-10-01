/** @file parameters.h
 * \brief Input parameter class. 
 * Big fuzzy object that contains 
 *   - input paramenter
 *   - Kac table 
 *   - structure constant 
 *   - verma module in each subspace
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "tensors.h" 
#include <complex>

/** Input parameter class. 
 * Big fuzzy object that contains 
 *   - input paramenter
 *   - Kac table 
 *   - structure constant 
 *   - verma module in each subspace
*/
class parameters
{
    public:
        // input parameters
        int p;               /**< p index in Kac matrix */
        int q;               /**< q index in Kac matrix */
        int numPert;         /**< number of pertubations */
        int maxLev;          /**< Maximum conformal level to include */
        int minLev;          /**< Minimum conformal level to include (usually 0) */
        int maxMomentum;     /**< Maximum momentum considered */
        bool allMomSubspace; /**< True to diagonalize all momentum subspaces together. Useful for calculation of momentum-carrying operator matrix elements.  False diagonalize momentum subspaces separately. */
        int maxEigen;        /**< Maximum number of eigenvalues printed in the output file. */
        int* r;              /**< array containing Kac r index of perturbations */
        int* s;              /**< array containing Kac s index of perturbations */
        double* deltas;      /**< array containing chiral conformal dimensions of each perturbation */
        double* pConst;      /**< array containing the constant for each perturbation */
        int rSteps;          /**< number of values of the radius to be considered */
        double rInc;         /**< increment of radius at each step */
        double rStart;       /**< initial value of the radius */
        int Nsize;           /**< basic sizes the NRG matrix. */
        int Nstep;           /**< number of states added at each NRG iteration. The matrix diagonalized in NRG is Nsize+Nstep */
        bool TCSA;           /**< True for diagonalization, false to have just Kac table structure constants and subspace analysis */
        bool NRG;            /**< True for NRG, false for full diagonalization */
        int nSweep;          /**< number of NRG sweep. Set 0 for basic NRG */
        bool SAVE_STATES;    /**< True to save computationale basis to file, false not to save them */
        bool READ_STATES;    /**< True to read computationale basis from file, false to calculate them */ 
        bool SAVE_ME;        /**< True to save chiral matrix element of perturbing fields to file, false not to save them */
        bool READ_ME;        /**< True to read chiral matrix element of perturbing fields to file, false to calculate them */
        bool calcExpVal;     /**< True to calculate expectation values on the eigenstates of the perturbed system of a primary operator, false not to */
        int rl_op;           /**< r index of holomorphic component of the operator to calculate expectation values */
        int sl_op;           /**< s index of holomorphic component of the operator to calculate expectation values */
        int nModes_l_op;     /**< number of virasoro operator applied to the holomorphic component of the primary operator (set 0 to evaluate operator itself) */
        int* mode_l_op;      /**< array containing the virasoro modes applied to the holomorphic component of the field */
        int rr_op;           /**< r index of anti holomorphic component of the operator to calculate expectation values */
        int sr_op;           /**< s index of anti holomorphic component of the operator to calculate expectation values */
        int nModes_r_op;     /**< number of virasoro operator applied to the anti holomorphic component of the primary operator (set 0 to evaluate operator itself) */
        int* mode_r_op;      /**< array containing the virasoro modes applied to the anti holomorphic component of the field */
        bool onlyDiagExpV;   /**< True to calculate only diagonal expectation values, false for the full matrix */
        bool saveEigenvec;   /**< True to save eigenvectors to file, false not to save them */
       
        // objects calculated from input parameters
        double C;            /**< central charge. Calculated by the code */
        int maxp;            /**< p size of Kac table. Calculated by the code */ 
        int maxq;            /**< q size of Kac table Calculated by the code*/
        int maxl;            /**< total size of Kac table Calculated by the code*/
        int numspace;        /**< number of subspace found by analysing the system Calculated by the code*/
        int* maxNumStates;   /**< array containng the maximum number of states for a given level Calculated by the code*/

        // tensors are first created as dummies. They are then properly allocated through assignment in paramters(ifstream)
        tensor2D<double> Kac;           /**< Contains the Kac table, meaning the weight of primary fields, indexed with standard r,s indeces */
        tensor6D< std::complex<long double> > cpt;  /**< Structure constant of the model. Calculated by the code */
        tensor3D< std::complex<long double> > cptsmall; /**< Structure constant of the model, indexed with field index instead of r,s. Calculated by the code */
        tensor3D< int > vermadim;       /**< contains the number of states for each level of each verma module. The first two indeces are the r,s indeces of the primary field, the third is the level */
        tensor3D< int > vermadimcomul;  /**< contains the cumulative number of states up to each level of each verma module. The first two indeces are the r,s indeces of the primary field, the third is the level */
        tensor2D< int > finalsubspace;  /**< contains the label of the verma module in each subspace. Calculated by the code */
        tensor2D< int > intMatrix;      /**< contains the sketch of the interaction matrix, where each element is a verma module. Calculated by the code */
        int* subdim;                    /**< contains the number of vermamodule in a give symmetry subspace. Calculated by the code  */
        int* leveldim;                  /**< Contains the */ 
        tensor2D<int> subdimStates;     /**< contains the predicted number of predicted non-chiral states in each subspace, and at given momentum */

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
        void analyse_subspaces();
        /** Destructor
         */
        ~parameters();
    private:
        /** Default constructor is private, to prevent instantiation of invalid object.
         */
        parameters();
};

#endif    
