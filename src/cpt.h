/** @file cpt.h
 *  \brief Set of helper function to calculate structure constant of Minimal CFTs.
 *
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date October 2017
 * \copyright GNU Public License.
 */
#ifndef CPT_H
#define CPT_H 
#include<complex>

/** Component of the structure constant formula, used by std::complex<long double> dot
 */
long double hrs(long double r, long double s, long double p, long double q); 
/** Component of the structure constant formula, used by std::complex<long double> dot
 */
long double dotmu(int ll, int l, long double ro, long double ro1);
/** Component of the structure constant formula, used by std::complex<long double> dot
 */
long double dotc(int ll, int l, long double ro, long double ro1, long double kn1, long double kn, long double ks1, long double ks, long double kp1, long double kp);
/** Component of the structure constant formula, used by std::complex<long double> dot
 */
long double dota(int kn1, int kn, long double ro, long double ro1);

/** Routine to calculate the structure constant of a Minimal CFT  M_(p,q) for a given triad of fields phi_(r1,s1) phi_(r2,s2) phi(r3,s3)
 * @param[in] p  first index of the minimal model 
 * @param[in] q  second index of the minimal model
 * @param[in] r1 first Kac index of the first field
 * @param[in] s1 second Kac index of the first field
 * @param[in] r2 first Kac index of the second field
 * @param[in] s2 second Kac index of the second field
 * @param[in] r3 first Kac index of the third field
 * @param[in] r3 second Kac index of the third field
 * return std::complex<long double>  containg the structure constant for the given fields.
 * \par <em> Beware that in the function call the s index (second index) comes before of the r index (first index) </em> for every field
 */
std::complex<long double> dot(int p, int q, int s1, int r1, int s2, int r2, int s3, int r3);
#endif



