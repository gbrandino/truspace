/** @file matrix_elements.h 
 * Class that contains the matrix elements of perturbing operators and operator to calculate expectation values 
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l., Robert Konik - BNL
 * last revision: October 2019 
 */

#ifndef MATRIX_ELEMENTS_H
#define MATRIX_ELEMENTS_H
#include "custom_parameters.h"
#include "states.h"
#include "tensors.h"
#include "filenames.h"

/** Class that contains the chiral matrix elements of perturbing operators and of the operator to calculate expectation values
 */
class matrix_elements 
{
    public:
        /** Constructor 
          * @param[in] p Reference to a paramters object
          * @param[in] st Reference to a states object
          */
        matrix_elements(custom_parameters& p, states& st);

        /** Destructor 
          */
        ~matrix_elements();

        /** Routine that either calculates or reads from file chiral matrix elemment
          * @param[in] p Reference to a paramters object
          * @param[in] fnames Reference to a filenames object
          * @param[in] st Reference to a states object
          */
        void calculate(custom_parameters& p, filenames& fnames, states& st);
 
        /** Routine that calculates chiral matrix element for the operator whose expectation values are asked 
          * @param[in] p Reference to a paramters object
          * @param[in] st Reference to a states object
          */
        void calculate_for_ExpVal(custom_parameters& p, states& st);
 
        /** Routine that calculates the expectation value on the eigenvectors of the perturbed system of the given operator
          * @param[in] p Reference to a paramters object
          * @param[in] st Reference to a states object
          * @param[in] coeff Matrix containing the eigenvectors of the perturbed system in rows
          * @param[in] size Number of components of each eigenvectors (size of the computational basis)
          * @param[in] fname Output file name
          * @param[in] value of the system size 
          */
        void calculate_ExpVal(parameters* p, states& st, double* coeff, int size, std::string fname, double R);

        tensor3D< double > me_unperturbed; /**< Contains the matrix element for the each of the perturbation */
        tensor2D< double > me_unperturbed_for_expval; /**< Contains the matrix element for the operator whose expecation values are asked */

    private:
        /** Default constructor is private, to prevent instantiation of invalid object.
         */
        matrix_elements();
};
#endif  
