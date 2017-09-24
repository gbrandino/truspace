/* 
 * Class that contains file name for the entire code
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l.
 * last revision: December 2016
 */

#ifndef FILENAMES_H
#define FILENAMES_H
#include "parameters.h"

class filenames
{
    public:
    
        filenames(parameters&, int);
        ~filenames();

        // outfiles and dump files
        std::string EigenName;          // File name of eigenvalues file
        std::string EigenNameFinal;     // File name of the eigenvalues for last NRG iteration
        std::string EigenVectName;      // File name of eigenvectors file
        std::string MeOpName;           // File name
        std::string StoreName;          // File name
        std::string OverlapName;        // File name
        std::string StateName;          // File name of states file 
        std::string* MeChirName;        // File name(s) of chiral matrix element file
    private:
        filenames();

};
#endif
