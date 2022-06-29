/** @file matrix_elements.h 
 * Class that contains the chiral matrix elements of perturbing operators and operator to calculate expectation values 
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l., Robert Konik - BNL
 * last revision: December 2016
 */

#ifndef MATRIX_ELEMENTS_H
#define MATRIX_ELEMENTS_H
#include "parameters.h"
#include "states.h"

/** Class that contains the chiral matrix elements of perturbing operators and of the operator to calculate expectation values
 */
class matrix_elements
{
    public:
        /** Constructor 
          * @param[in] p Reference to a paramters object
          * @param[in] st Reference to a states object
          * @param[in] subspace Index of the subspace to work on
          */
        matrix_elements(parameters& p, states& st, int subspace);
        /** Destructor 
          */
        virtual ~matrix_elements();

        /** Routine that calculates the expectation value on the eigenvectors of the perturbed system of the given operator
          * @param[in] p Reference to a paramters object
          * @param[in] st Reference to a states object
          * @param[in] coeff Matrix containing the eigenvectors of the perturbed system in rows
          * @param[in] size Number of components of each eigenvectors (size of the computational basis)
          * @param[in] fname Output file name
          * @param[in] R value of the radius R
          */
        virtual void calculate_ExpVal(parameters& p, states& st, double* coeff, int size, std::string fname, double R) = 0;


    protected:
        /** Default constructor is protected, to prevent instantiation of invalid object.
         */
        matrix_elements();
};
#endif  
