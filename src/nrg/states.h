/** @file states.h
 * \brief Class that contains and calculates chiral states by removing null vectors from each verma module
 *
 * \author Giuseppe Piero Brandino - eXact-lab s.r.l., Robert Konik - BNL
 * \date November 2017
 */

#ifndef STATES_H
#define STATES_H
#include "parameters.h"

/** Class that contains states 
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
        virtual ~states();

    	int nstate;       /**< The total number of state in the trucated Hilbert space*/


    protected:

        states();       
};
#endif 
