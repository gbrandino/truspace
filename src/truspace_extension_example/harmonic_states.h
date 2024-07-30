/** @file harmonic_states.h
 * \brief Class that contains and calculates harmonic_states 
 *
 * \author Giuseppe Piero Brandino - eXact-lab s.r.l., Robert Konik - BNL
 * \date October 2019
 */

#ifndef HARMONIC_STATES_H
#define HARMONIC_STATES_H
#include "states.h"
#include "custom_parameters.h"
#include "harmonic_filenames.h"


/** Class that contains and calculates chiral states by removing null vectors from each verma module
 * All the methods needs to be implemented, since they are virtual in the base class
 */
class harmonic_states : public states
{
    public:
        /** Constructor 
         * @param[in] p Reference to a parameters object
         */
        harmonic_states(custom_parameters& p);
        /** Destructor 
         */
        ~harmonic_states();
   
        /** Routine to calculate chiral states, stored internally in the object itself.
         * @param[in] p Reference to a parameters object
         */ 
        void calculate(custom_parameters& p);        
        /** Routine to calculate nonchiral states, stored internally in the object itself.
         * @param[in] p Reference to a parameters object
         * @param[in] subspace Subspace to work on
         * @param[in] mom Momentum subspace to work on
        */
        
        double* state_en; /**< Array contaning the unperturbed energy of every state */
        double* state_quantum_number; /**< Array contaning the quantum number that identifies the state */

    private:
        harmonic_states();

};
#endif 
