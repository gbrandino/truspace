/** harmonic_Hamiltonian class file. 
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
 */

#include "harmonic_Hamiltonian_MPI.h"
#include "custom_parameters.h"
#include "tensors.h"
#include "harmonic_matrix_elements.h"
#include "harmonic_states.h"
#include "custom_config.h"
#include <mpi.h>
#include "scalapack.h"
#include <cmath>
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>

harmonic_Hamiltonian::~harmonic_Hamiltonian()
{
    delete [] Offset;
}


harmonic_Hamiltonian::harmonic_Hamiltonian(custom_parameters& p_in, harmonic_states& st_in, harmonic_matrix_elements& me_in)
:
Hamiltonian(p_in, st_in, me_in)
{
    // init blacs grid
    int nb=NB;
    int iam;
    int nprocs;
    int izero = 0;
    Cblacs_pinfo( &iam, &nprocs ) ;
    Cblacs_get( -1, 0, &ictxt );
    nprow = npcol = sqrt(nprocs);
    Cblacs_gridinit( &ictxt, "Row", nprow, npcol );
    Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);


    if (!p.NRG)
        {
        Msize = st.nstate;
        rowH    = numroc_( &Msize, &nb, &myrow, &izero, &nprow );
        colH    = numroc_( &Msize, &nb, &mycol, &izero, &npcol );
        tensor2D< double> H_tmp(rowH, colH);
        H = H_tmp;

        localOffsetRow = 0;
        for (int proc=0; proc < myrow; proc++)
            {
            localOffsetRow +=  numroc_( &Msize, &nb, &proc, &izero, &nprow );
            }
        localOffsetCol = 0;
        for (int proc=0; proc < mycol; proc++)
            {
            localOffsetCol +=  numroc_( &Msize, &nb, &proc, &izero, &npcol );
            }
        }
    else
        {
        Msize = p.Nsize + p.Nstep;
        rowH    = numroc_( &Msize, &nb, &myrow, &izero, &nprow );
        colH    = numroc_( &Msize, &nb, &mycol, &izero, &npcol );
        tensor2D< double> H_tmp(rowH, colH);
        H = H_tmp;

        localOffsetRow = 0;
        for (int proc=0; proc < myrow; proc++)
            {
            localOffsetRow +=  numroc_( &Msize, &nb, &proc, &izero, &nprow );
            }
        localOffsetCol = 0;
        for (int proc=0; proc < mycol; proc++)
            {
            localOffsetCol +=  numroc_( &Msize, &nb, &proc, &izero, &npcol );
            }
        }
 
    // Offset contains the offsets for each process. The last element is the size of H
    Offset = new int[nprocs_mpi+1];
    int* localSizes=new int[nprocs_mpi];
    for (int i=0; i<nprocs_mpi; ++i)
        {
        if( i < Msize%nprocs_mpi )
            localSizes[i] = Msize/nprocs_mpi + 1;
        else
            localSizes[i] = Msize/nprocs_mpi;
        }
    
    int sumSizes = 0;
    Offset[0] = 0;   
    for (int i=1; i<nprocs_mpi; ++i)
        {
        sumSizes += localSizes[i-1];
        Offset[i] = sumSizes;
        }
    Offset[nprocs_mpi] = Msize;

    delete [] localSizes;

    // checking printing size 
    if ( p.maxEigen == 0 )
        {
        if ( p.NRG )
            {
            p.maxEigen = p.Nsize;
            }
        else
            {
            p.maxEigen = Msize;
            }
        }

    if ( p.NRG )
        {
        if (p.maxEigen > p.Nsize)
            {
            p.maxEigen = p.Nsize;
            }
        }
    else
        {
        if (p.maxEigen > Msize)
            {
            p.maxEigen = Msize;
            }
        }

}

void harmonic_Hamiltonian::init()
{
    custom_parameters& p_der = dynamic_cast<custom_parameters&>(p);
    harmonic_states& st_der = dynamic_cast<harmonic_states&>(st);
    harmonic_matrix_elements& me_der = dynamic_cast<harmonic_matrix_elements&>(me);

    for (int i=0;i<rowH*colH;i++)
         {
         H.data[i]=0.;
         }

    int izero=0;
    int nb=NB;
    int itemp=rowH>1?rowH:1;
    int info;
    descinit_( descH,  &Msize, &Msize, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info );

    int localSlice = Offset[myrank_mpi+1] - Offset[myrank_mpi];

    double* tempH = new double[localSlice * Msize];
    for (int i=0; i < localSlice*Msize; ++i)
        tempH[i]=0.;

    MPI_Barrier(MPI_COMM_WORLD);

    for (int pert=0; pert < p_der.numPert; pert++)
        {
        for(int i=Offset[myrank_mpi]; i < Offset[myrank_mpi+1]; ++i)
            {
            for(int j=0; j < Msize; ++j)
                    {
                    tempH[ (i- Offset[myrank_mpi])*Msize + j] += R*R * me_der.me_unperturbed(pert,i,j);
                    }
            // diagonal elements
            tempH[ (i- Offset[myrank_mpi])*Msize + i] += st_der.state_en[i]/(R*R)/p_der.numPert;
            }
        } 

    //Offset[1] always contains the size of the largest slice (if any of them is larger)
    double* shareH = new double[Offset[1]*Msize];
    for (int n=0; n<nprocs_mpi; ++n)
        {
        int slice = Offset[n+1] - Offset[n];
        if (myrank_mpi == n)
            {
            memcpy(shareH, tempH, slice*Msize*sizeof(double));
            }
        MPI_Bcast(shareH, slice*Msize, MPI_DOUBLE, n, MPI_COMM_WORLD); 
        for(int i=0; i < slice; ++i)
            {
            for(int j=0; j < Msize; ++j)
                {
                int row=i + Offset[n] + 1;
                int col=j + 1;
                double alpha = shareH[i*Msize +j];
                pdelset_(H.data, &row, &col, descH, &alpha);
                }
            }
        }

    delete [] tempH;
    delete [] shareH;

    // Hermiticity check
    for(int i=0;i<Msize;++i) 
        {
        for(int j=0;j<Msize;++j)
            {
            double alphai, alphaj;
            char getcode = 'A';
            char top = ' ';
            int row=i+1;
            int col=j+1;
            pdelget_(&getcode,&top, &alphai, H.data, &row, &col, descH);
            pdelget_(&getcode,&top, &alphaj, H.data, &col, &row, descH);
            //if (myrow == 0 and mycol == 0)
            //if (myrank_mpi == 0)
                //std::cout  <<  alphai << " ";
            if (fabs(alphai-alphaj) > .000001)
                std::cout << " Warning - hermiticity problem i,j H(i,j) H(j,i) " << i << " " << j << " " << alphai << " " << alphaj  << " " << std::endl;
            }
        //if (myrow == 0 and mycol == 0)
        //if (myrank_mpi == 0)
            //std::cout << std::endl;
        } 

}

void harmonic_Hamiltonian::NRG_block_A(double* eigenvalues)
{
    // called first, clears the entire Hamiltonian, then stores the eigenvalues on the diagonal 
    // for the Nsize*Nsize submatrix
    for (int i=0;i<Msize;++i)   
        { 
        for (int j=0; j < Msize; ++j)
            {
            int row = i + 1;
            int col = j + 1;
            double alpha = 0.0;
            if ( i == j )
                alpha = eigenvalues[i];
            pdelset_(H.data, &row, &col, descH, &alpha);
            }
        }

}

void harmonic_Hamiltonian::NRG_block_BC(int iteration, double* coeff)
{
    custom_parameters& p_der = dynamic_cast<custom_parameters&>(p);
    harmonic_states& st_der = dynamic_cast<harmonic_states&>(st);
    harmonic_matrix_elements& me_der = dynamic_cast<harmonic_matrix_elements&>(me);

    // all ranks calculate block B-C. This could probably be parallelized, even though the speed up could be irrelevant
    int shift = p.Nsize + p.Nstep*iteration;
    double* me_needed = new double[p.Nstep*shift];
    double* results = new double[p.Nstep*p.Nsize];
    for (int pert=0; pert < p_der.numPert; pert++)
        {
        for(int j=0;j < p.Nstep;++j) 
            {
            int jshift = shift + j;
            int jshift1 = j*shift;
            for(int k=0; k < shift;++k) 
                {
                me_needed[jshift1+k]= R*R*me_der.me_unperturbed(pert,jshift,k);
                }
            }

        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,p.Nstep,p.Nsize,shift,1.0,me_needed,shift,coeff,shift,0.0,results,p.Nsize);

        for(int i=0;i<p.Nstep;++i)
            {
            int ishift = i + p.Nsize;
            int ishift1 = i*p.Nsize;
            for(int j=0;j< p.Nsize;++j)
                {
                //H(ishift, j) += results[ishift1+j];
                //H(j,ishift) += results[ishift1+j];
                int row = ishift + 1;
                int col = j + 1;
                double alpha = results[ishift1+j];
                pdelset_(H.data, &row, &col, descH, &alpha);
                pdelset_(H.data, &col, &row, descH, &alpha);
                }
            }
        }

    delete [] results;    
    delete [] me_needed;
}

void harmonic_Hamiltonian::NRG_block_D(int iteration)
{
    custom_parameters& p_der = dynamic_cast<custom_parameters&>(p);
    harmonic_states& st_der = dynamic_cast<harmonic_states&>(st);
    harmonic_matrix_elements& me_der = dynamic_cast<harmonic_matrix_elements&>(me);

    double* tempD = new double[p.Nstep*p.Nstep];
    for(int i=0; i < p.Nstep*p.Nstep; ++i)
        tempD[i] = 0.0;

    // all ranks calculate block D
    for (int pert=0; pert < p_der.numPert; pert++)
        {
        int shift = p.Nsize + p.Nstep*iteration;
        for(int i=0; i < p.Nstep; ++i)
            {
            int ishift = shift + i; 

            for(int j=0; j < p.Nstep; ++j)
                {
                int jshift = shift + j;
                tempD[i*p.Nstep + j] += R*R * me_der.me_unperturbed(0,ishift, jshift);
                }
            // diagonal elements
            tempD[i*p.Nstep + i] += st_der.state_en[shift + i]/(R*R)/p_der.numPert; 
            }
        }

    // setting block D in the distributed H
    for(int i=0; i < p.Nstep; ++i)
        {
        for(int j=0; j < p.Nstep; ++j)
            {
            int row = i + p.Nsize + 1;
            int col = j + p.Nsize + 1;
            double alpha = tempD[i*p.Nstep + j];
            pdelset_(H.data, &row, &col, descH, &alpha);
            }
        }

    delete [] tempD;
}
 

void harmonic_Hamiltonian::sweep_block_BC(int iteration, double* coeff, double* coeffSweep, int NRGspacesize, int NRGstoresize)
{
    custom_parameters& p_der = dynamic_cast<custom_parameters&>(p);
    harmonic_states& st_der = dynamic_cast<harmonic_states&>(st);
    harmonic_matrix_elements& me_der = dynamic_cast<harmonic_matrix_elements&>(me);

    // all ranks calculate block B-C. This could probably be parallelized, the speed up in sweep could be useful
    double* H_part = new double[p.Nsize * NRGstoresize ];
    for (int i=0; i<p.Nsize*NRGstoresize; i++)
       H_part[i] = 0.0;

    for (int pert=0; pert < p_der.numPert; pert++)
        {

        for(int i=0;i < p.Nsize;++i)
            {
            //for(int j=0;j < p.Nstep;++j)
            //    {
                for(int k=0; k < NRGspacesize;++k )
                    {
                    for (int  k1=0; k1 < NRGstoresize; ++k1)
                        {
                        H_part[i*NRGstoresize + k1] += R*R * coeff[i*NRGspacesize + k] *me_der.me_unperturbed(pert,k,k1);
                        }
                    }
                }
            //}
        }

    double* H_part2 = new double[p.Nsize*p.Nstep ];

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,p.Nsize,p.Nstep,NRGstoresize,1.0,H_part,NRGstoresize,coeffSweep,NRGstoresize,0.0,H_part2,p.Nstep);

    delete [] H_part;

    //unperturbed contribution
    double* coefftemp= new double[p.Nsize*NRGstoresize];
    for (int i=0;i < p.Nsize;++i)
        {
        for (int k=0; k < NRGstoresize;k++)
            {
            coefftemp[i*NRGstoresize + k]=coeff[i*NRGspacesize+k]*st_der.state_en[k]/(R*R);
            }
        }

    double* H_part3=new double[p.Nsize*p.Nstep];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,p.Nsize,p.Nstep,NRGstoresize,1.0,coefftemp,NRGstoresize,coeffSweep,NRGstoresize,0.0,H_part3,p.Nstep);

    delete [] coefftemp;

    for (int i=0; i < p.Nsize;i++)
        {
        for (int j=0; j < p.Nstep;j++)
            {
            H_part2[i*p.Nstep + j] += H_part3[i*p.Nstep+j];
            }
        }
    delete [] H_part3;

    for (int i=0; i < p.Nsize;i++)
        {
        for (int j=0; j < p.Nstep;j++)
            {
            //H(i,j + p.Nsize) = H_part2[i*p.Nstep +j];
            //H(j +p.Nsize,i) = H_part2[i*p.Nstep +j];
            int row = i + 1;
            int col = j + p.Nsize + 1;
            double alpha = H_part2[i*p.Nstep +j];
            pdelset_(H.data, &row, &col, descH, &alpha);
            pdelset_(H.data, &col, &row, descH, &alpha);
            }
        }

    delete [] H_part2;
}



void harmonic_Hamiltonian::print_eigenValues(int iteration, int sweep, double* eigenValues, std::string fname, int maxState)
{

    std::ofstream eigen_output;
    eigen_output.open(fname.c_str(), std::ofstream::out | std::ofstream::app);

    eigen_output << iteration  << " " << sweep << " " << std::setiosflags(std::ios::fixed) << std::setprecision(14) << " " << R << " ";

    for (int i=0; i < p.maxEigen; ++i)
        {
        eigen_output  << eigenValues[i] << " ";
        }
    eigen_output << std::endl;
    

}

void harmonic_Hamiltonian::print_eigenVectors(double* eigenVectors, std::string fname, int size)
{

    if(myrank_mpi == 0)
        {
        std::ofstream eigenVec_output;
        eigenVec_output.open(fname.c_str(), std::ofstream::out | std::ofstream::app);

        for (int i=0; i < p.maxEigen; ++i)
            {
            for (int j=0; j < size; ++j)
                {
                double alpha = 0;
                char getcode = 'A';
                char top = ' ';
                int row=i+1;
                int col=j+1;
                pdelget_(&getcode,&top, &alpha, H.data, &row, &col, descH); 
                eigenVec_output  << alpha  << " ";
                }
            eigenVec_output << std::endl;
            }
        }
    MPI_Barrier(MPI_COMM_WORLD);
}
