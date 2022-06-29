/** @file element_computation.h
 * \brief Routines to calculate matrix element and overlap of conformal states
 
 * \author  Robert Konik - BNL, Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date    November 2016
 */
#ifndef ELEMENT_COMPUTATION_H
#define ELEMENT_COMPUTATION_H
#include "minimal_model_config.h"

/** Return the factorial of the given integer
 * @param[in] n The integer to calculate the factorial
 * return The factorial of n
 */ 
long double fac(int n);

/** Calculates the matrix element of a primary field phi_(hr,hs) on the conformal basis.
 * @param[in] C The central charge of the model
 * @param[in] grade Array of size MaxNumGrade (defined in minimal_model_config.h). Contains the sequence of indeces of lowering and rising operator that define the bra and ket states, and 0 for the primary field, read from left to right
 * @param[in] n Length of the string or rising and lowering operators, including the primary field
 * @param[in] m Position of the primary field in the string 
 * @param[in] hr Weight of the primary field from which the current conformal ket states is descending
 * @param[in] hl Weight of the primary field from which the current conformal bra states is descending
 * @param[in] coeff Coefficient to multiply the matrix element. Used internally in recursive calls
 * @param[in] Del The conformal weight of the field whose matrix element are computed. 0 indicates identity, so the overlaps of the conformal basis states are computed. 
 * return The matrix element
 */
long double element_computation(long double C, long double grade[MaxNumGrade],int n,int m,long double hr,long double hl,long double coeff,long double Del);
#endif
