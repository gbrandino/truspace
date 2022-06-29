/** Hamiltonian class file. 
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date June 2022
 */


#include "Hamiltonian.h"

Hamiltonian::Hamiltonian(parameters& p_in, states& st_in, matrix_elements& me_in)
:
p(p_in),
st(st_in),
me(me_in)
{

}

Hamiltonian::~Hamiltonian()
{

}
