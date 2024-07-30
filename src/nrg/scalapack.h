#ifndef SCALAPACK_H
#define SCALAPACK_H

extern "C" {
void pdsyev_( char *jobz, char *uplo, int *n, double *a, int *ia, int *ja, int *desca, double *w, double *z, int *iz, int *jz, int *descz, double *work, int *lwork, int *info );
void pdsyevd_( char *jobz, char *uplo, int *n, double *a, int *ia, int *ja, int *desca, double *w, double *z, int *iz, int *jz, int *descz, double *work, int *lwork, int *iwork, int* liwork, int *info );
void pdelget_(char *scope, char *top, double* alpha, double *a, int* ia, int* ja, int* desca); 
void pdelset_( double *a, int* ia, int* ja, int* desca, double* alpha);
int numroc_(int *n, int *nb, int *iproc, int *isproc, int *nprocs);
void descinit_(int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld, int* info);
void pdgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *A, int *IA, int *JA, int *DESCA, double *B, int *IB, int *JB, int *DESCB, double *BETA, double *C, int *IC, int *JC, int *DESCC);
void   Cblacs_pinfo( int* mypnum, int* nprocs);
void   Cblacs_get( int context, int request, int* value);
int    Cblacs_gridinit( int* context, const char * order, int np_row, int np_col);
void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
void   Cblacs_gridexit( int context);
void   Cblacs_exit( int error_code);
}

#endif
