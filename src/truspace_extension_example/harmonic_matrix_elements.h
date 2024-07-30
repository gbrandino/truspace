/** @file harmonic_matrix_elements.h 
 * Class that contains the matrix elements of perturbing operators and operator to calculate expectation values 
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l., Robert Konik - BNL
 * last revision: October 2019 
 */

#ifndef HARMONIC_MATRIX_ELEMENTS_H
#define HARMONIC_MATRIX_ELEMENTS_H
#include "custom_parameters.h"
#include "matrix_elements.h"
#include "harmonic_states.h"
#include "tensors.h"
#include "harmonic_filenames.h"

/** Class that contains the chiral matrix elements of perturbing operators and of the operator to calculate expectation values
 * All the methods needs to be implemented, since they are virtual in the base class
 */
class harmonic_matrix_elements : public matrix_elements
{
    public:
        /** Constructor 
          * @param[in] p Reference to a paramters object
          * @param[in] st Reference to a states object
          */
        harmonic_matrix_elements(custom_parameters& p, harmonic_states& st);

        /** Destructor 
          */
        ~harmonic_matrix_elements();

        /** Routine that either calculates or reads from file chiral matrix elemment
          * @param[in] p Reference to a paramters object
          * @param[in] fnames Reference to a filenames object
          * @param[in] st Reference to a states object
          */
        void calculate(custom_parameters& p, harmonic_filenames& fnames, harmonic_states& st);
 
        /** Routine that calculates chiral matrix element for the operator whose expectation values are asked 
          * @param[in] p Reference to a paramters object
          * @param[in] st Reference to a states object
          */
        void calculate_for_ExpVal(custom_parameters& p, harmonic_states& st);
 
        /** Routine that calculates the expectation value on the eigenvectors of the perturbed system of the given operator
          * @param[in] p Reference to a paramters object
          * @param[in] st Reference to a states object
          * @param[in] coeff Matrix containing the eigenvectors of the perturbed system in rows
          * @param[in] size Number of components of each eigenvectors (size of the computational basis)
          * @param[in] fname Output file name
          * @param[in] R Value of the system size 
          */
        void calculate_ExpVal(parameters& p, states& st, double* coeff, int size, std::string fname, double R);

        tensor3D< double > me_unperturbed; /**< Contains the matrix element for the each of the perturbation */
        tensor2D< double > me_unperturbed_for_expval; /**< Contains the matrix element for the operator whose expecation values are asked */

    private:
        /** Default constructor is private, to prevent instantiation of invalid object.
         */
        harmonic_matrix_elements();
};
#endif  
