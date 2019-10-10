/** @file states.h
 * \brief Class that contains and calculates states 
 *
 * \author Giuseppe Piero Brandino - eXact-lab s.r.l., Robert Konik - BNL
 * \date October 2019
 */

#ifndef STATES_H
#define STATES_H
#include "parameters.h"
#include "filenames.h"

/** Class that contains and calculates chiral states by removing null vectors from each verma module
 */
class states
{
    public:
        /** Constructor 
         * @param[in] p Reference to a parameters object
         */
        states(parameters& p);
        /** Destructor 
         */
        ~states();
   
        /** Routine to calculate chiral states, stored internally in the object itself.
         * @param[in] p Reference to a parameters object
         * @param[in] subspace Subspace to work on
         * @param[in] fname Reference to a filenames object
         */ 
        void calculate(parameters& p);        
        /** Routine to calculate nonchiral states, stored internally in the object itself.
         * @param[in] p Reference to a parameters object
         * @param[in] subspace Subspace to work on
         * @param[in] mom Momentum subspace to work on
        */
        
            int nstate;    /**< Size of the truncated Hilbert Space **/
            double* state_en; /**< Array contaning the unperturbed energy of every state */
            int* state_quantum_number; /**< Array contaning the quantum number that identifies the state */

    private:
        states();

};
#endif 
