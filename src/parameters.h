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
#include "sparse_tensor.h"
#include <complex>
#include <vector>
#include <string>
#include <map>

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
        std::pair<double,double>* deltas;      /**< array containing chiral conformal dimensions of each perturbation */
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
        int r_op;           /**< r index of anti holomorphic component of the operator to calculate expectation values */
        int s_op;           /**< s index of anti holomorphic component of the operator to calculate expectation values */
        bool onlyDiagExpV;   /**< True to calculate only diagonal expectation values, false for the full matrix */
        bool saveEigenvec;   /**< True to save eigenvectors to file, false not to save them */
        std::string fields_file;  /**< file name of the fields file */
        std::string constants_file; /**< file name of the structure constants file */

        // objects calculated from input parameters
        double C;            /**< central charge. Calculated by the code */
        int maxp;            /**< p size of Kac table. Calculated by the code */ 
        int maxq;            /**< q size of Kac table Calculated by the code*/
        int maxl;            /**< total size of Kac table Calculated by the code*/
        int n_fields;        /**< number of non-chiral fields read from file */
        int numspace;        /**< number of subspace found by analysing the system Calculated by the code*/
        int currentSubspace; /**< contains the subspace currently being analysed*/
        int* maxNumStates;   /**< array containng the maximum number of states for a given level Calculated by the code*/

        std::vector<std::pair< std::pair< int, int >, std::pair< int, int> > > fields; /** vector that contains pairs of pairs of indexes, such the a non chiral field is indexed by a pair (r,s) (r',s') for the chiral and anti chiral component respectively, as read from input file */
        std::map< int, std::pair< double, double> > fields_dimensions; /** map that contains the non chiral fields conformal dimesions. Given the index of a non chiral field, as from the input file, it return a pair of double, the chiral and anti chiral conformal dimension of the corresponding field */

        // tensors are first created as dummies. They are then properly allocated through assignment in paramters(ifstream)
        tensor2D<double> Kac;           /**< Contains the Kac table, meaning the weight of primary fields, indexed with standard r,s indeces */
        sparse_tensor6D<double> cpt;  /**< Structure constant of the model. Calculated by the code */
        //tensor3D< std::complex<long double> > cptsmall; /**< Structure constant of the model, indexed with field index instead of r,s. Calculated by the code */
        sparse_tensor3D<double> cptsmall; /**< Structure constants of the model, indexed with field index. Read from input file */
        std::vector< std::map < std::pair<int, int>, int > > PairToIndex; /**< Vector whose elements are map of pairs of indexes of fields (r,s) and the corresponding linear index of the chiral field, foreach subspace */
        std::vector< std::map < int, std::pair<int, int> > > IndexToPair; /**< Vector whose elements are map of the linear index of the chiral field and the correspoding pairs of indexes of fields (r,s), foreach subspace */
        tensor3D< int > vermadim;       /**< contains the number of states for each level of each verma module. The first two indeces are the r,s indeces of the primary field, the third is the level */
        tensor3D< int > vermadimcomul;  /**< contains the cumulative number of states up to each level of each verma module. The first two indeces are the r,s indeces of the primary field, the third is the level */
        tensor2D<int> finalsubspace;  /**< contains the index of the non chiral fields in each subspace. Calculated by the code */
        tensor2D< int > intMatrix;      /**< contains the sketch of the interaction matrix, where each element is a the tensor product of left and right verma module. Calculated by the code */
        int* subdim;                    /**< contains the number of vermamodule in a give symmetry subspace. Calculated by the code  */
        int* leveldim;                  /**< Contains the */ 
        tensor2D<int> subdimStates;     /**< contains the predicted number of predicted non-chiral states in each subspace, and at given momentum */

        // class methods
        /** Constructor
         * @param[in] The name of the model input file
         * @param[in] The name of the nrg input file
         */
        parameters(const char*, const char*);
        
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
