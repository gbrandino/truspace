/* Set of helper functions manage ginac computation 
 * 
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l.
 * last revision: October 2016
 *          
 */
#ifndef CHARACTERS_H
#define CHARACTERS_H
#include <ginac/ginac.h>
#include "tensors.h"
using namespace std;
using namespace GiNaC;

ex invphi(int, symbol);
ex ksymb1(int, int, int, int, int, symbol);
ex ksymb2(int, int, int, int, int, symbol);
ex chisymb(int, int, int, int,int, symbol);
void characexpan(int, int, double, double, int, int*);
int partitions(int);
int determine_partition(int, tensor2D<int>&);
void  part(int, int*, int, tensor2D<int>&, int&);


#endif
