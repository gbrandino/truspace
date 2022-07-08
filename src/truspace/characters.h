/** @file characters.h 
 *  \brief Set of helper functions to calculate integers partitions and manage GiNaC computation of character expansion. 
 *
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date October 2017
 * \copyright GNU Public License.
 */
#ifndef CHARACTERS_H
#define CHARACTERS_H
#include <ginac/ginac.h>
#include "tensors.h"
using namespace std;
using namespace GiNaC;

/** Calculates the inverse denominator of the characters formula up to order k.
 *  @param[in] k maximum exponent included in the product in the denominator
 *  @param[in] x GiNac symbol object, contains the variable used for the symbolic calculation
 *  return GiNaC expression containing the symbolic inverse denominator up to order k
 */
ex invphi(int k, symbol x);

/** Calculates the first contribution of the numerator of the characters formula up to order k.
 *  @param[in]   r The first index defining the Verma Module
 *  @param[in]   s The second index defining the Verma Module
 *  @param[in]   p The first index defining the model M_(p,q)
 *  @param[in]   q The second index defining the model M_(p,q)
 *  @param[in]   k  maximum exponent included in the product in the denominator
 *  @param[in]   x GiNac symbol object, contains the variable used for the symbolic calculation
 *  return        GiNaC expression containing the symbolic first contribution to the numerator up to order k
 */
ex ksymb1(int r, int s, int p, int q, int k, symbol x);

/** Calculates the second contribution of the numerator of the characters formula up to order k.
 *  @param[in]   r The first index defining the Verma Module
 *  @param[in]   s The second index defining the Verma Module
 *  @param[in]   p The first index defining the model M_(p,q)
 *  @param[in]   q The second index defining the model M_(p,q)
 *  @param[in]   k  maximum exponent included in the product in the denominator
 *  @param[in]   x GiNac symbol object, contains the variable used for the symbolic calculation
 *  return        GiNaC expression containing the symbolic second contribution to the numerator up to order k
 */
ex ksymb2(int r, int s, int p, int q, int k, symbol x);


/** Calculates the symbolic expression of the characters formula for given Verma module and model, up to order k.
 *  @param[in]   r The first index defining the Verma Module
 *  @param[in]   s The second index defining the Verma Module
 *  @param[in]   p The first index defining the model M_(p,q)
 *  @param[in]   q The second index defining the model M_(p,q)
 *  @param[in]   l  maximum exponent included in the product in the denominator
 *  @param[in]   x GiNac symbol object, contains the variable used for the symbolic calculation
 *  return        GiNaC expression containing the characters for the given Verma module and model up to order k
 */
ex chisymb(int r, int s, int p, int q, int l, symbol x);

/** Calculates the characters expantion, to get the number of independent states at each level for a given Verma module.
 *  @param[in]   r The first index defining the Verma Module
 *  @param[in]   s The second index defining the Verma Module
 *  @param[in]   p The first index defining the model M_(p,q)
 *  @param[in]   q The second index defining the model M_(p,q)
 *  @param[in]   l Maximum level to perform the expansion
 *  @param[out] leveldim Array of lenght l+1, contains the number of indenpendent states at each level
 *
 */
void characexpan(int r, int s, double p, double q, int l, int* leveldim);


/** Calculates the number of partions of the given integer.
 *  @param[in]  k The integer to calculate the partion
 *  return  Number of partions of the given integer
 */
int partitions(int k);

/** Calculates the partions of integers up to the given integer 
 *  @param[in]   level integer to calculate the partions
 *  @param[out]  results Reference to tensor2D<int> that will contain the different partions of the given integer, one per row 
 *  return       The number of partions for the given integer
 */
int determine_partition(int level, tensor2D<int>& results);

/** Helper function, recursive component of determine_partition
 */
void  part(int, int*, int, tensor2D<int>&, int&);


#endif
