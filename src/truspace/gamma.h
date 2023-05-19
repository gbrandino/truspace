/** @file gamma.h 
 *  \brief Computation of gamma function.
 *  Algorithms and coefficient values from "Computation of Special Functions", Zhang and Jin, John Wiley and Sons, 1996.
 
 *  \authors C. Bond. All rights reserved.
 *  \date 2003
 */

/** Computation of gamma function.
 * @params[in] x Value to calculate gamma(x)
 * return gamma function of argument 'x'.
 * \par <em> NOTE: Returns 1e308 if argument is a negative integer or 0, or if argument exceeds 171. </em> 
 */
double gamma_custom(double x);
