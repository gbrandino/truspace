/** @file custom_config.h
 *  \brief Configuration file holding version number, some hard-coded limits that effect array dimensions, constants and macros
 *
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date October 2019
 * \copyright GNU Public License.
 */


#ifndef CONFIG_H
#define CONFIG_H
#define VERSION "1.0.0"

#define SIGN(a,b) ((b) >= 0.0 ? fabsl(a) : -fabsl(a))
#define SIGNd(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define PI 3.1415926535897 
#ifdef MPI_PAR
#define NB 16  /* blocking factor for BLACS/SCALAPACK */
#endif
#endif
