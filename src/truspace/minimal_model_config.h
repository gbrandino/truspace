/** @file minimal_model_config.h
 *  \brief Configuration file holding version number, some hard-coded limits that effect array dimensions, constants and macros
 *
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date October 2017
 * \copyright GNU Public License.
 */


#ifndef MINIMAL_MODEL_CONFIG_H
#define MINIMAL_MODEL_CONFIG_H
#define VERSION "2.0.0"

/** \def MaxPart
 *  \brief Maximum number of partitions of integers ( so the maxium possible number of different strings of lowering opertaros ) in a chiral expectation value.
 *  It could be in principle be calculated, but for any practical purposes given value is enough.
 *  Increasing this will typically require to increase stack size.
 */


/** \def MaxNumGrade
 *  \brief Maximum number of operators ( both field and lowering operator ) in a chiral expectation value.
 *  It could be in principle be calculated, but for any practical purposes given value is enough.
 *  Increasing this will typically require to increase stack size.
 */


#define MaxPart  1500 
#define MaxNumGrade  50 
#define SIGN(a,b) ((b) >= 0.0 ? fabsl(a) : -fabsl(a))
#define SIGNd(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define pi 3.1415926535897 

#endif
