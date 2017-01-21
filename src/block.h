/*
 *  Routine to identify symmetry subspace by the analysis of the interaction matrix   
 *
 * authors: Giuseppe Piero Brandino - eXact-lab s.r.l.
 * last revision: November 2016
 */

#include "tensors.h"

void  blockrec(int maxl,tensor2D< complex <long double> >& a, tensor2D< int>& finalsubspace,int *subdim,int &numspace);
