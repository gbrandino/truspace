/** @file filenames.h
 * \brief Output filenames helper class. 
 *   - Contains the filenames helper class 
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
 */

#ifndef FILENAMES_H
#define FILENAMES_H
#include "minimal_model_parameters.h"

/** Output filenames helper class. 
 *  Initializes and holds all the output file names, based on input paramters
*/
class filenames
{
    public:
     
        /** Constructor
         * @param[in] A reference to a parameters object
         */
        filenames(minimal_model_parameters&, int);
        ~filenames();

        // outfiles and dump files
        std::string EigenName;          /**< File name of eigenvalues file */
        std::string EigenNameFinal;     /**< File name of the eigenvalues for last NRG iteration */
        std::string EigenVectName;      /**< File name of eigenvectors file */
        std::string MeOpName;           /**< File name of Expectation value on the final eigenstates*/
        std::string StateName;          /**< File name of states file */
        std::string OverlapName;        /**< File name of file holding the overlap of states */
        std::string* MeChirName;        /**< File name(s) of chiral matrix element file */
    private:
        /** Default constructor is private, to prevent instantiation of invalid object.
         */
        filenames();

};
#endif
