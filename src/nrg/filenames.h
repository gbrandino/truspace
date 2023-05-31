/** @file filenames.h
 * \brief Output filenames helper class. 
 *   - Contains the filenames helper class 
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
 */

#ifndef FILENAMES_H
#define FILENAMES_H
#include "parameters.h"

/** Output filenames helper class. 
 *  Initializes and holds all the output file names, based on input paramters
*/
class filenames
{
    public:
     
        /** Constructor
         * @param[in] p A reference to a parameters object
         */
        filenames(parameters& p);
        ~filenames();

        // outfiles and dump files
        std::string EigenName;          /**< File name of eigenvalues file */
        std::string EigenNameFinal;     /**< File name of the eigenvalues for last NRG iteration */
        std::string EigenVectName;      /**< File name of eigenvectors file */
        std::string MeOpName;           /**< File name of Expectation value on the final eigenstates*/
    protected:
        /** Default constructor is protected, to prevent instantiation of invalid object.
         */
        filenames();

};
#endif
