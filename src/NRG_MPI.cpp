/** Class that handles exact diagonalization, NRG and Sweep in a model-agnostic way ( as much as possible ) - MPI Parallel Version using scalapack/pblas 
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
 */

#include "NRG_MPI.h"
#include "parameters.h"
#include "tensors.h"
#include "filenames.h"
#include "config.h"
#include "scalapack.h"
#include <cblas.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <mpi.h>

NRG::NRG()
{

}

NRG::~NRG()
{

}

void NRG::run(parameters& p, Hamiltonian& H, filenames& fnames)
{
    int myrank_mpi,nprocs_mpi;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);

    H.init();
    char jobz;
    char uplo='L';
    int n = H.Msize;
    int ia, ja;
    int iz,jz;
    int descZ[9];
    int izero=0;
    int ione=1;
    int lwork;
    double* work;
    int info;

    double* eigen = new double[H.Msize];

    if ( p.NRG or p.calcExpVal or p.saveEigenvec)
        {
        jobz = 'V';
        tensor2D<double> EigenVectors_tmp(H.rowH, H.colH);
        EigenVectors = EigenVectors_tmp;
        for(int i=0; i<9; i++)
            descZ[i]=H.descH[i];
        }
    else
        {
        jobz = 'N';
        }

    NSweeps= 0;
    if (!p.NRG)
         {
         NIters = 0;       
         NSweeps = 0; 
         }
    else
         {
         NIters = (H.st.nstate-p.Nsize-p.Nstep)/p.Nstep;
         if (p.nSweep>0)
            NSweeps = p.nSweep;
         }

    if ( p.NRG )
        {
        if (myrank_mpi == 0)
            {
            std::cout << " Performing " << NIters << " iterations " << std::endl;
            std::cout << " NRG iteration: 0 " << std::flush;
            }            
        }
    
    lwork = -1;
    work = new double[2];

    //
 
    int liwork = 7*n + 8*H.npcol +2;
    int* iwork = new int[liwork];
    jobz = 'V';
        tensor2D<double> EigenVectors_tmp(H.rowH, H.colH);
        EigenVectors = EigenVectors_tmp;
        for(int i=0; i<9; i++)
            descZ[i]=H.descH[i];
    //
 
    //pdsyev_( &jobz, &uplo, &n, H.H.data, &ione, &ione, H.descH, eigen, EigenVectors.data, &ione, &ione, descZ, work, &lwork, &info ); 
    pdsyevd_( &jobz, &uplo, &n, H.H.data, &ione, &ione, H.descH, eigen, EigenVectors.data, &ione, &ione, descZ, work, &lwork, iwork, &liwork, &info );    

    lwork = (int)work[0];
    delete [] work;  
    work = new double[lwork];

    //pdsyev_( &jobz, &uplo, &n, H.H.data, &ione, &ione, H.descH, eigen, EigenVectors.data, &ione, &ione, descZ, work, &lwork, &info );
    pdsyevd_( &jobz, &uplo, &n, H.H.data, &ione, &ione, H.descH, eigen, EigenVectors.data, &ione, &ione, descZ, work, &lwork, iwork, &liwork, &info );
    if (info > 0 and myrank_mpi == 0)
        {
        std::cerr << " SCALapack problem. Error code: " << info << std::endl;
        }

    if (myrank_mpi == 0)    
        H.print_eigenValues(0,0,eigen,fnames.EigenName, H.Msize);
    MPI_Barrier(MPI_COMM_WORLD);

    if ( p.saveEigenvec and !p.NRG )
        {
        H.print_eigenVectors(EigenVectors.data, fnames.EigenVectName, H.Msize);
        }
    if ( p.calcExpVal and !p.NRG)
        {
        H.me.calculate_ExpVal(p, H.st, EigenVectors.data, H.Msize, H.Momentum, fnames.MeOpName, H.R);
        }


    // starting iterations 
    double * coeff; 
    if ( p.NRG )
        {
        coeff=new double[p.Nsize*H.Msize];
        for(int i=0;i < p.Nsize;++i)
            {
            int ishift = i*H.Msize;
            for(int k=0; k < H.Msize;++k)
                 {
                 coeff[ishift+k] = EigenVectors(i,k);  //the i-ith eigenvector is stored as the i-th row of X
                 }
            }

        }

    std::fstream storeCoeff;
    std::fstream storeEigen;
    double* coeffSweep;

    for ( int sweep = 0; sweep <= NSweeps; ++sweep)
        {
        if (sweep ==  0 and NSweeps > 0)
            {
            //Storing to file eigenvectors and eigenvalues of 'step' part of hamiltonian for iteration 0
            std::stringstream storeCoeffname;
            storeCoeffname << "storeCoeff_0" << ".bin";
            storeCoeff.open(storeCoeffname.str().c_str(),std::ios::out | std::ios::binary | std::ios::trunc );
            for (int i= p.Nsize; i< H.Msize;i++)
                {
                for (int j=0;j< H.Msize;j++)
                    {
                    storeCoeff.write(reinterpret_cast<char*>( &EigenVectors.data[i*(H.Msize)+j] ),sizeof(double));
                    }
                }
            storeCoeff.close();

            std::stringstream storeEigenname;
            storeEigenname << "storeEigen_0"<< ".bin";
            storeEigen.open(storeEigenname.str().c_str(), std::ios::out| std::ios::binary | std::ios::trunc);
            for (int i= p.Nsize; i< H.Msize;i++)
                {
                storeEigen.write(reinterpret_cast<char*>( &eigen[i] ),sizeof(double));
                }
            storeEigen.close();
            }
        if ( sweep >0 )
            std::cout << " Sweep " << sweep << " : "<< std::flush;


        // starting NRG iterations
        for ( int iter = 0; iter <= NIters; ++iter)
            {
            // iteration 0 of sweep 0 has already been performed
            if ( sweep == 0 and iter == 0)
                continue;


            std::cout << iter << " " << std::flush;
            
            // imagine Hamiltonian in block form:
            //
            //       H = (A B)
            //           (C D)
            //
            if (sweep == 0)
                {
                for (int i = 0; i< H.Msize*H.Msize; ++i)
                    H.H.data[i]=0.;

                H.NRG_block_A(eigen);
                H.NRG_block_D(iter);
                H.NRG_block_BC(iter, coeff);
                }
            else
                {
                for (int i = 0; i< H.Msize*H.Msize; ++i)
                    H.H.data[i]=0.;

                H.NRG_block_A(eigen);

                // BLOCK D: while sweeping, block D simply contains the eigenvalues between the correspoding eigenstates of the previous iteration

                std::stringstream storeEigenname;
                storeEigenname << "storeEigen_"<< iter << ".bin";
                storeEigen.open(storeEigenname.str().c_str(), std::ios::in | std::ios::binary );
                storeEigen.seekg(0, std::ios::beg);
                for (int i = p.Nsize; i < H.Msize;i++)
                    {
                    storeEigen.read(reinterpret_cast<char*>( &H.H.data[i*H.Msize + i] ),sizeof(double));
                    }
                storeEigen.close();
                // BLOCK B and C: while sweeping, block B and C require the eigenvectors calculated at the previous iterations

                std::stringstream storeCoeffname;
                storeCoeffname << "storeCoeff_" << iter << ".bin";
                storeCoeff.open(storeCoeffname.str().c_str(),std::ios::in | std::ios::binary);
                storeCoeff.seekg(0, std::ios::beg);
                int NRGstoresize;
                int NRGspacesize;
                if ( sweep == 1)
                    {
                    NRGstoresize = p.Nsize + p.Nstep*(iter+1);
                    }
                else
                    {
                    NRGstoresize = p.Nsize + p.Nstep*(NIters+1);
                    }
                NRGspacesize = p.Nsize + p.Nstep*(NIters+1);
                coeffSweep = new double[p.Nstep*NRGstoresize];
                       
                for (int i=0; i < p.Nstep;i++)
                    {
                    for(int j=0;j < NRGstoresize;j++)
                        {
                        storeCoeff.read(reinterpret_cast<char*>( &coeffSweep[i*NRGstoresize+j] ),sizeof(double));
                        }
                    }
                storeCoeff.close();

                H.sweep_block_BC(iter, coeff, coeffSweep, NRGspacesize, NRGstoresize );

                }

            ///diagonalize Hamiltonian and sort for relevant states
            //dsyevr_(&jobz,&range,&uplo,&n,H.H.data,&lda,&vl,&vu,&il,&iu,&abstol,&ms,eigen,EigenVectors.data,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);

            int NRGsize;
            int newSize;
            if ( sweep == 0)
                {
                NRGsize = p.Nsize + p.Nstep*(iter);
                newSize = NRGsize + p.Nstep;
                }
            else
                {
                NRGsize = p.Nsize + p.Nstep*(NIters+1);
                newSize = NRGsize;
                }

            H.print_eigenValues(sweep,iter,eigen,fnames.EigenName, newSize);
            if ( sweep== NSweeps and iter==NIters)
                {
                H.print_eigenValues(sweep,iter,eigen,fnames.EigenNameFinal, newSize);
                }
            //compute coefficients of states in new basis

            // Determine how the coefficients of the first 'Nsize + step*iter'
            // basis states have changed

            double* coeff_part;
            coeff_part=new double[H.Msize*NRGsize];
 
            double* EigenVectorsPart = new double[H.Msize*p.Nsize];

            for(int i=0; i < H.Msize;++i)
                {
                int ishift = i*p.Nsize;
                int ishift1 = i*H.Msize;
                double* pt = &EigenVectorsPart[ishift];
                double* pt_part = &EigenVectors.data[ishift1];
                memcpy(pt,pt_part,p.Nsize*sizeof(double));
                }
    
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,H.Msize,NRGsize,p.Nsize,1.0,EigenVectorsPart,p.Nsize,coeff,NRGsize,0.0,coeff_part,NRGsize);

            delete [] coeff;
            delete [] EigenVectorsPart;

            coeff=new double[p.Nsize*newSize];

            if ( sweep == 0 )
                {
                for(int i=0; i < p.Nsize;++i)
                    {
                    double* pt = &coeff[i*newSize];
                    double* pt_part = &coeff_part[i*NRGsize];
                    memcpy(pt,pt_part,NRGsize*sizeof(double));
                    }

                // Determine how the coefficients of the last 'Nstep'
                // basis states have changed
                for(int i=0; i < p.Nsize;++i)
                    {
                    int ishift = i*newSize + NRGsize;
                    int ishift1 = i*H.Msize + p.Nsize;
                    double* pt = &coeff[ishift];
                    double* pt_part = &EigenVectors.data[ishift1];
                    memcpy(pt,pt_part,p.Nstep*sizeof(double));
                    }
                // if sweeping, save last 'Nstep' part of the eigenvector, on the conformal basis
                if ( NSweeps > 0 )
                    {
                    double* coeff_part2= new double[p.Nstep*newSize];

                    for(int i=0; i < p.Nstep;++i)
                        {
                        int ishift = (i+p.Nsize)*NRGsize; 
                        double* pt = &coeff_part2[i*newSize];
                        double* pt_part = &coeff_part[ishift];
                        memcpy(pt,pt_part,NRGsize*sizeof(double));
                        }

                    for(int i=0; i < p.Nstep;++i)
                        {
                        int ishift = i*newSize + NRGsize;
                        int ishift1 = (i + p.Nsize)*H.Msize + p.Nsize;
                        double* pt = &coeff_part2[ishift];
                        double* pt_part = &EigenVectors.data[ishift1];
                        memcpy(pt,pt_part,p.Nstep*sizeof(double));
                        }            
                    std::stringstream storeCoeffname;
                    storeCoeffname << "storeCoeff_" << iter << ".bin";
                    storeCoeff.open(storeCoeffname.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc );
                    for (int i = 0; i < p.Nstep;i++)
                        {
                        for (int j=0; j < newSize;j++)
                            {
                            storeCoeff.write(reinterpret_cast<char*>( &coeff_part2[i*newSize + j] ),sizeof(double));
                            }
                        }
                    storeCoeff.close();
                    
                    delete [] coeff_part2;
                    }
                delete [] coeff_part;
                }
            else
                {
                // while sweeping, the last Nstep coefficient of the eigenvectors are expressed on the basis of the eigenvectors found at the previous iteration
                double* EigenVectorsPart2 = new double[H.Msize*p.Nstep];
                for(int i=0; i < H.Msize;++i)
                    {
                    int ishift = i*p.Nstep;
                    int ishift1 = i*H.Msize + p.Nsize;
                    double* pt = &EigenVectorsPart2[ishift];
                    double* pt_part = &EigenVectors.data[ishift1];
                    memcpy(pt,pt_part,p.Nstep*sizeof(double));                  
                    }
                int NRGstoresize;
                if ( sweep == 1)
                    {
                    NRGstoresize = p.Nsize + p.Nstep*(iter+1);
                    }
                else
                    {
                    NRGstoresize = NRGsize;
                    }

                double* coeff_part2=new double[H.Msize*NRGstoresize];
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,H.Msize,NRGstoresize,p.Nstep,1.0,EigenVectorsPart2,p.Nstep,coeffSweep,NRGstoresize,0.0,coeff_part2,NRGstoresize);
                delete [] coeffSweep;
                for ( int i = 0; i < H.Msize; ++i)
                    {
                    for (int  j = 0; j < NRGstoresize; ++j)
                        {
                        coeff_part[i*NRGsize +j] += coeff_part2[i*NRGstoresize +j];
                        }
                    }
                for(int i=0; i < p.Nsize;++i)
                    {
                    int ishift = i*NRGsize;
                    double* pt = &coeff[ishift];
                    double* pt_part = &coeff_part[ishift];
                    memcpy(pt,pt_part,NRGsize*sizeof(double));
                    }

                std::stringstream storeCoeffname;
                storeCoeffname << "storeCoeff_" << iter << ".bin";
                storeCoeff.open(storeCoeffname.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc );
                for (int i = p.Nsize; i < H.Msize;i++)
                    {
                    for (int j=0; j < NRGsize;j++)
                        {
                        storeCoeff.write(reinterpret_cast<char*>( &coeff_part[i*NRGsize+j] ),sizeof(double));
                        }
                    }
                storeCoeff.close();

                delete [] coeff_part;
                delete [] coeff_part2;
                delete [] EigenVectorsPart2;
                }
            //Storing to file eigenvalues of 'step' part of hamiltonian
            if (NSweeps > 0)
                {

                std::stringstream storeEigenname;
                storeEigenname << "storeEigen_" << iter << ".bin";
                storeEigen.open(storeEigenname.str().c_str(), std::ios::out | std::ios::binary| std::ios::trunc);
                for (int j= p.Nsize; j < H.Msize;j++)
                    {
                    storeEigen.write(reinterpret_cast<char*>( &eigen[j] ),sizeof(double));
                    }

                storeEigen.close();
                }

            if ( sweep== NSweeps and iter==NIters)
                {
                if ( p.saveEigenvec )
                    H.print_eigenVectors(coeff, fnames.EigenVectName, newSize);

                if ( p.calcExpVal )
                    H.me.calculate_ExpVal(p, H.st, coeff, newSize, H.Momentum, fnames.MeOpName, H.R);
    
                }

            } // closing loop over NRG iterations
        if (myrank_mpi == 0)
            std::cout << std::endl;
        } // closing loop over sweeps


    if (p.NRG)
        delete [] coeff;

    delete [] work;
    delete [] iwork;
    delete [] eigen;
}
