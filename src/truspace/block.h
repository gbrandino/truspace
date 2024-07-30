/** @file block.h
 *  \brief Helper function identify symmetry subspace by the analysis of the interaction matrix   
 *
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date October 2017
 * \copyright GNU Public License.
 */

#include "tensors.h"

/** Routine to identify symmetry subspace by the analysis of the interaction matrix
 *  @param[in] maxl  size of the interaction matrix a
 *  @param[in] a     matrix of type tensor2D<double>>containing the constant (structure costant) of the various verma module and perturbation, of which we want to identify subspace
 *  @param[out] finalsubspace tensor2D< int> that will contain in the row the indeces of the Verma module belonging the each subspace
 *  @param[out] subdim Array that contains the number of Verma module in each subspace
 *  @param[out] numspace Number of subspace found 
 */
void  blockrec(int maxl,tensor2D< double >& a, tensor2D< int>& finalsubspace, int *subdim,int &numspace);
