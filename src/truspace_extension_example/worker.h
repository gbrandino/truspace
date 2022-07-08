/** @file worker.h 
 * \brief Contains driver Class for states, matrix element and NRG calculation. 
 *  It takes care of 
 *    - calling NRG routine for each value of the changing paramters ( system size L in this case) 
 *    - calculating Matrix element on the eigenstates found 
 *
 * \author Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
 */

#ifndef WORKER_H
#define WORKER_H

#include "custom_parameters.h"

/** Driver Class for states, matrix element and NRG calculation. 
 *  It takes care of 
 *    - calling NRG routine for each value of the changing paramters ( system size L in this case) 
 *    - calculating Matrix element on the eigenstates found 
 *
 * \author Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date October 2019
 */
class worker
{
    public:
        /** Default constructor 
        */
        worker();
        /** Default destructor 
        */
        ~worker();    
        /** 
         * @param[in] p paramerts object
         */    
        void run(custom_parameters& p);

};
#endif
