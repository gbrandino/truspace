/** @file minimal_model_filenames.h
 * \brief Output filenames helper class specilized for minimal models. 
 *   - Contains the filenames helper class 
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
 */

#ifndef MINIMAL_MODEL_FILENAMES_H
#define MINIMAL_MODEL_FILENAMES_H
#include "filenames.h"
#include "minimal_model_parameters.h"

/** Output filenames helper class. 
 *  Initializes and holds all the output file names, based on input paramters
*/
class minimal_model_filenames: public filenames
{
    public:
     
        /** Constructor
         * @param[in] p A reference to a parameters object
         */
        minimal_model_filenames(minimal_model_parameters& p);
        ~minimal_model_filenames();

        std::string* MeChirName;        /**< File name(s) of chiral matrix element file */
    private:
        /** Default constructor is private, to prevent instantiation of invalid object.
         */
        minimal_model_filenames();

};
#endif
