/* 
 * Class that contains the chiral matrix elements of perturbing operators and operator to calculate expectation values 
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l., Robert Konik - BNL
 * last revision: December 2016
 */

#ifndef CHIRME_H
#define CHIRME_H
#include "parameters.h"
#include "states.h"
#include "tensors.h"
#include "filenames.h"

class chirME
{
    public:
        chirME(parameters&, states&, int);
        ~chirME();
        void calculate_or_read_chirME(parameters&, filenames&, states&, int);
        void calculate_chirME_for_ExpVal(parameters&, states&, int);
        void calculate_ExpVal(parameters&, states&, double*, int, int, std::string, double);

        tensor5D<long double> me_chi_ortho;
        tensor4D<long double> me_chi_ortho_op_l;
        tensor4D<long double> me_chi_ortho_op_r;

    private:
        chirME();
};
#endif  
