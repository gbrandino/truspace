/** @file config.h
 *  \brief Configuration file holding scalapack compilation-time parameters 
 *
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date October 2017
 * \copyright GNU Public License.
 */


#ifndef CONFIG_H_SCALAPACK
#define CONFIG_H_SCALAPACK
#ifdef MPI_PAR
#define NB 16  /* blocking factor for BLACS/SCALAPACK */
#endif
#endif
