/** @file chirME.h 
 * Class that contains the chiral matrix elements of perturbing operators and operator to calculate expectation values 
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l., Robert Konik - BNL
 * last revision: December 2016
 */

#ifndef CHIRME_H
#define CHIRME_H
#include "parameters.h"
#include "states.h"
#include "tensors.h"
#include "filenames.h"

/** Class that contains the chiral matrix elements of perturbing operators and of the operator to calculate expectation values
 */
class chirME
{
    public:
        /** Constructor 
          * @param[in] p Reference to a paramters object
          * @param[in] st Reference to a states object
          * @param[in] subspace Index of the subspace to work on
          */
        chirME(parameters& p, states& st, int subspace);
        /** Destructor 
          */
        ~chirME();
        /** Routine that either calculates or reads from file chiral matrix elemment
          * @param[in] p Reference to a paramters object
          * @param[in] fnames Reference to a filenames object
          * @param[in] st Reference to a states object
          * @param[in] subspace Index of the subspace to work on
          */

        void calculate_or_read_chirME(parameters& p, filenames& fnames, states& st, int subspace);
 
        /** Routine that calculates chiral matrix element for the operator whose expectation values are asked 
          * @param[in] p Reference to a paramters object
          * @param[in] st Reference to a states object
          * @param[in] subspace Index of the subspace to work on
          */
        void calculate_chirME_for_ExpVal(parameters& p, states& st, int subspace);
 
        /** Routine that calculates the expectation value on the eigenvectors of the perturbed system of the given operator
          * @param[in] p Reference to a paramters object
          * @param[in] st Reference to a states object
          * @param[in] coeff Matrix containing the eigenvectors of the perturbed system in rows
          * @param[in] size Number of components of each eigenvectors (size of the computational basis)
          * @param[in] Momentum Momentum subspace
          * @param[in] fname Output file name
          * @param[in] R value of the radius R
          */
        void calculate_ExpVal(parameters& p, states& st, double* coeff, int size, int Momentum, std::string fname, double R);

        tensor5D< std::pair<long double, long double> > me_chi_ortho; /**< Contains the chiral matrix element for each perturbing operator. Indexed as (perturbation_number, left Verma mod, right Verma mod, left states number, right states number) */
        tensor4D< std::pair<long double, long double> > me_chi_ortho_op; /**< Contains the chiral matrix element for the operator whose expecation values are asked. Indexed as (left Verma mod, right Verma mod, left state number, right state number) */

    private:
        /** Default constructor is private, to prevent instantiation of invalid object.
         */
        chirME();
};
#endif  
