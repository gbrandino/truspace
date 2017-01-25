#ifndef CONFIG_H
#define CONFIG_H
#define VERSION "1.0.0"
// these constants encode the maximum number of partitions of integers ( so the possible number of different strings of lowering opertaros 
// and the maximum number of operators ( both field and lowering operator ) in a chiral expectation value. They could be in principle be calculated,
// but for any practical purpose the constant below are enough.
#define MaxPart  2000
#define MaxNumGrade  200 
#define SIGN(a,b) ((b) >= 0.0 ? fabsl(a) : -fabsl(a))
#define SIGNd(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define pi 3.1415926535897 
#endif
