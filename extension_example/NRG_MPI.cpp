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
    std::cout << "ssdfs " << lwork << std::endl;
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
        if (myrank_mpi == 0)
            H.me.calculate_ExpVal(p, H.st, EigenVectors.data, H.Msize, fnames.MeOpName, H.R);
        MPI_Barrier(MPI_COMM_WORLD);
        }


    // starting iterations 
    // Each MPI process holds a copy of the eigenvectors of the previous iteration
    double * coeff; 
    if ( p.NRG )
        {
        coeff=new double[p.Nsize*H.Msize];
        for(int i=0; i < p.Nsize;++i)
            {
            int ishift = i*H.Msize;
            for(int k=0; k < H.Msize;++k)
                 {
                 char getcode = 'A';
                 char top = ' ';
                 int row=i+1;
                 int col=k+1;
                 double alpha = 0.0;
                 pdelget_(&getcode,&top, &alpha, EigenVectors.data, &col, &row, descZ); //the i-th eigenvector is stored as the i-th col of EigenVectors.data. That's why we switch col and row
                 coeff[ishift + k] = alpha;  
                 }
            }

        }
    std::fstream storeCoeff;
    std::fstream storeEigen;
    double* coeffSweep;

    for ( int sweep = 0; sweep <= NSweeps; ++sweep)
        {
        if (sweep ==  0 and NSweeps > 0 )
            {
            //Storing to file eigenvectors and eigenvalues of 'step' part of hamiltonian for iteration 0
            // Only rank 0 menages I/O
            if (myrank_mpi == 0)
                {
                std::stringstream storeCoeffname;
                storeCoeffname << "storeCoeff_0" << ".bin";
                storeCoeff.open(storeCoeffname.str().c_str(),std::ios::out | std::ios::binary | std::ios::trunc );
                }
            for (int i= p.Nsize; i< H.Msize;i++)
                {
                for (int j=0;j< H.Msize;j++)
                    {
                    char getcode = 'A';
                    char top = ' ';
                    int row = i+1;
                    int col = j+1;
                    double alpha;
                    pdelget_(&getcode,&top, &alpha, EigenVectors.data, &col, &row, H.descH);  //the i-th eigenvector is stored as the i-th col of EigenVectors.data. That's why we switch col and row
                    if (myrank_mpi == 0)
                        storeCoeff.write(reinterpret_cast<char*>( &alpha ),sizeof(double));
                    }
                }
            if (myrank_mpi == 0)
                storeCoeff.close();

            if (myrank_mpi == 0)
                {
                std::stringstream storeEigenname;
                storeEigenname << "storeEigen_0"<< ".bin";
                storeEigen.open(storeEigenname.str().c_str(), std::ios::out| std::ios::binary | std::ios::trunc);
                for (int i= p.Nsize; i< H.Msize;i++)
                    {
                    storeEigen.write(reinterpret_cast<char*>( &eigen[i] ),sizeof(double));
                    }
                storeEigen.close();
                }
            }
        MPI_Barrier(MPI_COMM_WORLD);

        if ( sweep >0 and myrank_mpi == 0)
            std::cout << " Sweep " << sweep << " : "<< std::flush;


        // starting NRG iterations
        for ( int iter = 0; iter <= NIters; ++iter)
            {
            // iteration 0 of sweep 0 has already been performed
            if ( sweep == 0 and iter == 0)
                continue;

            if ( myrank_mpi == 0 )
                std::cout << iter << " " << std::flush;
            
            // imagine Hamiltonian in block form:
            //
            //       H = (A B)
            //           (C D)
            //
            if (sweep == 0)
                {
                //for (int i = 0; i< H.rowH*H.colH; ++i)
                //    H.H.data[i]=0.;

                H.NRG_block_A(eigen);
                H.NRG_block_D(iter);
                H.NRG_block_BC(iter, coeff);
                }
            else
                {
//                for (int i = 0; i< H.Msize*H.Msize; ++i)
//                    H.H.data[i]=0.;

                H.NRG_block_A(eigen);

                // BLOCK D: while sweeping, block D simply contains the eigenvalues between the correspoding eigenstates of the previous iteration
                double* tempD = new double[p.Nstep];
                if ( myrank_mpi == 0 )  
                    {
                    std::stringstream storeEigenname;
                    storeEigenname << "storeEigen_"<< iter << ".bin";
                    storeEigen.open(storeEigenname.str().c_str(), std::ios::in | std::ios::binary );
                    storeEigen.seekg(0, std::ios::beg);
                    for (int i = 0; i < p.Nstep; i++)
                        {
                        storeEigen.read(reinterpret_cast<char*>( &tempD[i] ),sizeof(double));
                        }
                    storeEigen.close();
                    }
                MPI::COMM_WORLD.Bcast(tempD, p.Nstep, MPI_DOUBLE, 0);

                for(int i=0; i < p.Nstep; ++i)
                    {
                    for(int j=0; j < p.Nstep; ++j)
                        {
                        int row = i + p.Nsize + 1;
                        int col = j + p.Nsize + 1;
                        double alpha = 0.0;
                        if ( i == j )
                            alpha = tempD[i];
                        pdelset_(H.H.data, &row, &col, H.descH, &alpha);
                        }
                    }
                delete [] tempD;
                // BLOCK B and C: while sweeping, block B and C require the eigenvectors calculated at the previous iterations
        
                if (myrank_mpi == 0)
                    {
                    std::stringstream storeCoeffname;
                    storeCoeffname << "storeCoeff_" << iter << ".bin";
                    storeCoeff.open(storeCoeffname.str().c_str(),std::ios::in | std::ios::binary);
                    storeCoeff.seekg(0, std::ios::beg);
                    }

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

                if (myrank_mpi == 0)
                    {                       
                    for (int i=0; i < p.Nstep;i++)
                        {
                        for(int j=0;j < NRGstoresize;j++)
                            {
                            storeCoeff.read(reinterpret_cast<char*>( &coeffSweep[i*NRGstoresize+j] ),sizeof(double));
                            }
                        }
                    storeCoeff.close();
                    }
                int Bcast_size = p.Nstep*NRGstoresize;
                MPI::COMM_WORLD.Bcast(coeffSweep, Bcast_size, MPI_DOUBLE, 0);            

                H.sweep_block_BC(iter, coeff, coeffSweep, NRGspacesize, NRGstoresize );

                }// closing else of if (sweep == 0)

            ///diagonalize Hamiltonian and sort for relevant states
            //dsyevr_(&jobz,&range,&uplo,&n,H.H.data,&lda,&vl,&vu,&il,&iu,&abstol,&ms,eigen,EigenVectors.data,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
            pdsyevd_( &jobz, &uplo, &n, H.H.data, &ione, &ione, H.descH, eigen, EigenVectors.data, &ione, &ione, descZ, work, &lwork, iwork, &liwork, &info );

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


            if (myrank_mpi == 0)
                H.print_eigenValues(sweep,iter,eigen,fnames.EigenName, newSize);

            if ( sweep== NSweeps and iter==NIters and myrank_mpi == 0)
                {
                H.print_eigenValues(sweep,iter,eigen,fnames.EigenNameFinal, newSize);
                }
            MPI_Barrier(MPI_COMM_WORLD);
            //compute coefficients of states in new basis

            // Determine how the coefficients of the first 'Nsize + step*iter'
            // basis states have changed

            double* coeff_part;
            coeff_part=new double[H.Msize*NRGsize];
 
            double* EigenVectorsPart = new double[H.Msize*p.Nsize];

            for(int i=0; i < H.Msize;++i)
                {
                int ishift = i*p.Nsize;
                //int ishift1 = i*H.Msize;
                //double* pt = &EigenVectorsPart[ishift];
                //double* pt_part = &EigenVectors.data[ishift1];
                //memcpy(pt,pt_part,p.Nsize*sizeof(double));
                for (int j=0; j < p.Nsize; ++j)
                    {
                    char getcode = 'A';
                    char top = ' ';
                    double alpha;
                    int row = i + 1;
                    int col = j + 1;
                    pdelget_(&getcode, &top, &alpha, EigenVectors.data, &col, &row, descZ);  //the i-th eigenvector is stored as the i-th col of EigenVectors.data. That's why we switch col and row
                    EigenVectorsPart[ishift + j] = alpha;
                    }       
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
                    //int ishift1 = i*H.Msize + p.Nsize;
                    //double* pt = &coeff[ishift];
                    //double* pt_part = &EigenVectors.data[ishift1];
                    //memcpy(pt,pt_part,p.Nstep*sizeof(double));
                    for (int j=0; j < p.Nstep; ++j)
                        {
                        char getcode = 'A';
                        char top = ' ';
                        double alpha;
                        int row = i + 1;
                        int col = j + p.Nsize + 1;
                        pdelget_(&getcode,&top, &alpha, EigenVectors.data, &col, &row, descZ);  //the i-th eigenvector is stored as the i-th col of EigenVectors.data. That's why we switch col and row
                        coeff[ishift + j] = alpha;
                        }
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
                        //int ishift1 = (i + p.Nsize)*H.Msize + p.Nsize; 
                        //double* pt = &coeff_part2[ishift];
                        //double* pt_part = &EigenVectors.data[ishift1];
                        //memcpy(pt,pt_part,p.Nstep*sizeof(double));

                        for (int j=0; j < p.Nstep; ++j)
                            {
                            char getcode = 'A';
                            char top = ' ';
                            double alpha;
                            int row = i + p.Nsize + 1;
                            int col = j + p.Nsize + 1;
                            pdelget_(&getcode,&top, &alpha, EigenVectors.data, &col, &row, descZ);  //the i-th eigenvector is stored as the i-th col of EigenVectors.data. That's why we switch col and row
                            coeff_part2[ishift + j] = alpha;
                            }
                        }

                    if (myrank_mpi == 0)
                        {            
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
                        }
                    MPI_Barrier(MPI_COMM_WORLD);
                    
                    delete [] coeff_part2;
                    }//closing if ( NSweeps > 0 )
                delete [] coeff_part;
                }
            else
                {
                // while sweeping, the last Nstep coefficient of the eigenvectors are expressed on the basis of the eigenvectors found at the previous iteration
                double* EigenVectorsPart2 = new double[H.Msize*p.Nstep];
                for(int i=0; i < H.Msize;++i)
                    {
                    char getcode = 'A';
                    char top = ' ';
                    double alpha;
                    int ishift = i*p.Nstep;
                    //int ishift1 = i*H.Msize + p.Nsize;
                    //double* pt = &EigenVectorsPart2[ishift];
                    //double* pt_part = &EigenVectors.data[ishift1];
                    //memcpy(pt,pt_part,p.Nstep*sizeof(double));                 
                    for (int j=0; j < p.Nstep; ++j) 
                        {
                        int row = i + 1;
                        int col = j + p.Nsize + 1;
                        pdelget_(&getcode,&top, &alpha, EigenVectors.data, &col, &row, descZ);  //the i-th eigenvector is stored as the i-th col of EigenVectors.data. That's why we switch col and row
                        EigenVectorsPart2[ishift + j] = alpha;
                        }
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

                if (myrank_mpi == 0)
                    {
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
                    }
                MPI_Barrier(MPI_COMM_WORLD);

                delete [] coeff_part;
                delete [] coeff_part2;
                delete [] EigenVectorsPart2;
                }
            //Storing to file eigenvalues of 'step' part of hamiltonian
            if (NSweeps > 0)
                {
                if (myrank_mpi == 0)
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
                MPI_Barrier(MPI_COMM_WORLD);
                }

            if ( sweep== NSweeps and iter==NIters)
                {
                if ( p.saveEigenvec )
                    H.print_eigenVectors(coeff, fnames.EigenVectName, newSize);

                if ( p.calcExpVal)
                    {
                    if (myrank_mpi == 0)                
                        H.me.calculate_ExpVal(p, H.st, coeff, newSize, fnames.MeOpName, H.R);
                    MPI_Barrier(MPI_COMM_WORLD);
                    }
    
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
