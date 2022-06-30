/** @file harmonic_filenames.h
 * \brief Output filenames helper class. 
 *   - Contains the filenames helper class 
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date October 2019 
 */

#ifndef HARMONIC_FILENAMES_H
#define HARMONIC_FILENAMES_H
#include "custom_parameters.h"
#include "filenames.h"
#include <string>

/** Output filenames helper class. 
 *  Initializes and holds all the output file names, based on input paramters
*/
class harmonic_filenames : public filenames
{
    public:
     
        /** Constructor
         * @param[in] A reference to a parameters object
         */
        harmonic_filenames(custom_parameters&);
        ~harmonic_filenames();

    private:
        /** Default constructor is private, to prevent instantiation of invalid object.
         */
        harmonic_filenames();

};
#endif
