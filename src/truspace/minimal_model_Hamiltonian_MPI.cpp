/** Hamiltonian class file. 
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
 */

#include "minimal_model_Hamiltonian_MPI.h"
#include "minimal_model_parameters.h"
#include "tensors.h"
#include "chirME.h"
#include "minimal_model_states.h"
#include "minimal_model_config.h"
#include <mpi.h>
#include "scalapack.h"
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#endif
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>

minimal_model_Hamiltonian::~minimal_model_Hamiltonian()
{
    delete [] Offset;
}


minimal_model_Hamiltonian::minimal_model_Hamiltonian(minimal_model_parameters& p_in, minimal_model_states& st_in, chirME& me_in)
:
Hamiltonian(p_in, st_in, me_in)
{
    minimal_model_states& st_der = dynamic_cast<minimal_model_states&>(st);
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
        Msize = st_der.nstate;
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

void minimal_model_Hamiltonian::init()
{
    minimal_model_parameters& p_der = dynamic_cast<minimal_model_parameters&>(p);
    minimal_model_states& st_der = dynamic_cast<minimal_model_states&>(st);
    chirME& me_der = dynamic_cast<chirME&>(me);

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
        //effective strength of each perturbation
        double c = p_der.pConst[pert]*2.*PI*pow(R/2./PI,1.-(p_der.deltas[pert].first + p_der.deltas[pert].second)); 

        //for(int i=0; i < Msize; ++i) 
        for(int i=Offset[myrank_mpi]; i < Offset[myrank_mpi+1]; ++i)
            {
            int st_vmod_left = st_der.state_vmod_left[i];
            int st_vmod_right = st_der.state_vmod_right[i];
            int st_left = st_der.state_left[i];
            int st_right  = st_der.state_right[i];
            std::pair <int,int> vmod_left_i  = p_der.IndexToPair[p_der.currentSubspace][st_vmod_left];
            std::pair <int,int> vmod_right_i = p_der.IndexToPair[p_der.currentSubspace][st_vmod_right];

            //off-diagonal elements
            for(int j=0; j < Msize; ++j)
            //for(int j=localOffsetCol; j < localOffsetCol + colH; ++j) 
                {
                if( ( st_der.cstate_lev(st_der.state_vmod_left[i],st_der.state_left[i]) - \
                      st_der.cstate_lev(st_der.state_vmod_right[i],st_der.state_right[i]) )== \
                    ( st_der.cstate_lev(st_der.state_vmod_left[j],st_der.state_left[j]) - \
                      st_der.cstate_lev(st_der.state_vmod_right[j],st_der.state_right[j]) ) )
                    {
                    int st_vmodj_left = st_der.state_vmod_left[j];
                    int st_vmodj_right = st_der.state_vmod_right[j];
                    std::pair <int,int> vmod_left_j  = p_der.IndexToPair[p_der.currentSubspace][st_vmodj_left];
                    std::pair <int,int> vmod_right_j = p_der.IndexToPair[p_der.currentSubspace][st_vmodj_right];
                    int field_i = -1;
                    int field_j = -1;
                    for (int f = 0; f < p_der.n_fields; f++)
                    {
                        if (p_der.fields[f] == std::make_pair(vmod_left_i,vmod_right_i) ){
                            field_i = f+1;
                        }
                        if (p_der.fields[f] == std::make_pair(vmod_left_j,vmod_right_j) ){
                            field_j = f+1;
                        }
                    }
                    if (field_i == -1){
                        std::cout << " Error, field not found! " << field_i << " (" << vmod_left_i.first << "," << vmod_left_i.second << ") (" << vmod_right_i.first << "," << vmod_right_i.second << ")"<< std::endl;
                    }
                    if (field_j == -1){
                        std::cout << " Error, field not found! " << field_j << " (" << vmod_left_j.first << "," << vmod_left_j.second << ") (" << vmod_right_j.first << "," << vmod_right_j.second <<  ")"<< std::endl;
                    }
                    long double SCtmp = p_der.cptsmall(field_i, p_der.r[pert], field_j);

                    double temp = me_der.me_chi_ortho(pert, st_vmod_left, st_vmodj_left, st_left, st_der.state_left[j]).first;
                    temp *= me_der.me_chi_ortho(pert, st_vmod_right, st_vmodj_right, st_right, st_der.state_right[j]).second;
                    tempH[ (i- Offset[myrank_mpi])*Msize + j] += SCtmp * c * temp;
                    }
                if ( i==j )
                    tempH[ (i-Offset[myrank_mpi])*Msize + i] += 2.*PI/R*(st_der.state_en[i]-p_der.C/12.) / p_der.numPert;
                }
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
            //    std::cout  <<  alphai << " ";
            if (fabs(alphai-alphaj) > .000001)
                std::cout << " Warning - hermiticity problem i,j H(i,j) H(j,i) " << i << " " << j << " " << alphai << " " << alphaj  << " " << std::endl;
            }
        //if (myrow == 0 and mycol == 0)
        //        std::cout << std::endl;
        } 

}

void minimal_model_Hamiltonian::NRG_block_A(double* eigenvalues)
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

void minimal_model_Hamiltonian::NRG_block_BC(int iteration, double* coeff)
{
    minimal_model_parameters& p_der = dynamic_cast<minimal_model_parameters&>(p);
    minimal_model_states& st_der = dynamic_cast<minimal_model_states&>(st);
    chirME& me_der = dynamic_cast<chirME&>(me);
    // all ranks calculate block B-C. This could probably be parallelized, even though the speed up could be irrelevant
    int shift = p.Nsize + p.Nstep*iteration;
    double* me_needed = new double[p.Nstep*shift];
    double* results = new double[p.Nstep*p.Nsize];
    for (int pert=0; pert < p_der.numPert; pert++)
        {
        double c = p_der.pConst[pert]*2.*PI*pow(R/2./PI,1.-(p_der.deltas[pert].first + p_der.deltas[pert].second));
        
        for(int j=0;j < p.Nstep;++j) 
            {
            int jshift = shift + j;
            int jshift1 = j*shift;
            for(int k=0; k < shift;++k) 
                {
                double temp;
                int st_vmod_left = st_der.state_vmod_left[jshift];
                int st_vmod_right = st_der.state_vmod_right[jshift];
                int st_left = st_der.state_left[jshift];
                int st_right = st_der.state_right[jshift];
                int st_vmodk_left = st_der.state_vmod_left[k];
                int st_vmodk_right = st_der.state_vmod_right[k];

                if( (st_der.cstate_lev(st_vmod_left,st_left) - st_der.cstate_lev(st_vmod_right, st_right)) == \
                    (st_der.cstate_lev(st_vmodk_left,st_der.state_left[k]) - st_der.cstate_lev(st_vmodk_right, st_der.state_right[k])))
                    {
                    std::pair <int,int> vmod_left_j  = p_der.IndexToPair[p_der.currentSubspace][st_vmod_left];
                    std::pair <int,int> vmod_right_j = p_der.IndexToPair[p_der.currentSubspace][st_vmod_right];
                    std::pair <int,int> vmod_left_k  = p_der.IndexToPair[p_der.currentSubspace][st_vmodk_left];
                    std::pair <int,int> vmod_right_k = p_der.IndexToPair[p_der.currentSubspace][st_vmodk_right];

                    int field_j = -1;
                    int field_k = -1;
                    for (int f = 0; f < p_der.n_fields; f++)
                    {
                        if( p_der.fields[f] == std::make_pair(vmod_left_j, vmod_right_j) ){
                            field_j = f+1;
                        }
                        if( p_der.fields[f] == std::make_pair(vmod_left_k, vmod_right_k) ){
                            field_k = f+1;
                        }
                    }
                    if (field_j == -1 || field_k == -1){
                        std::cout << " Error, fields not found! " << field_j << " " << field_k << std::endl;
                    }

                    long double SCtmp = p_der.cptsmall(field_j, p_der.r[pert], field_k);

                    temp = me_der.me_chi_ortho(pert,st_vmod_left,st_vmodk_left,st_left,st_der.state_left[k]).first;
                    temp *= SCtmp * me_der.me_chi_ortho(pert,st_vmod_right,st_vmodk_right,st_right,st_der.state_right[k]).second;
                    }
                else
                    {
                    temp=0.;
                    }
                me_needed[jshift1+k]= c*temp;
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

void minimal_model_Hamiltonian::NRG_block_D(int iteration)
{
    minimal_model_parameters& p_der = dynamic_cast<minimal_model_parameters&>(p);
    minimal_model_states& st_der = dynamic_cast<minimal_model_states&>(st);
    chirME& me_der = dynamic_cast<chirME&>(me);

    double* tempD = new double[p.Nstep*p.Nstep];
    for(int i=0; i < p.Nstep*p.Nstep; ++i)
        tempD[i] = 0.0;

    // all ranks calculate block D
    for (int pert=0; pert < p_der.numPert; pert++)
        {
        //effective strength of each perturbation
        double c = p_der.pConst[pert]*2.*PI*pow(R/2./PI,1.- (p_der.deltas[pert].first + p_der.deltas[pert].second) );
        int shift = p.Nsize + p.Nstep*iteration;
        for(int i=0; i < p.Nstep; ++i)
            {
            int ishift = shift + i; 
            int st_vmod_left = st_der.state_vmod_left[ishift];
            int st_vmod_right = st_der.state_vmod_right[ishift];
            int st_left = st_der.state_left[ishift];
            int st_right= st_der.state_right[ishift];
            std::pair <int,int> vmod_left_i  = p_der.IndexToPair[p_der.currentSubspace][st_vmod_left];
            std::pair <int,int> vmod_right_i = p_der.IndexToPair[p_der.currentSubspace][st_vmod_right];

            //off-diagonal elements
            for(int j=0; j < p.Nstep; ++j)
                {
                int jshift = shift + j;
                if( ( st_der.cstate_lev(st_vmod_left, st_left) - st_der.cstate_lev(st_vmod_right, st_right) )== \
                    ( st_der.cstate_lev(st_der.state_vmod_left[jshift],st_der.state_left[jshift]) - \
                      st_der.cstate_lev(st_der.state_vmod_right[jshift],st_der.state_right[jshift]) ) )
                    {
                    int st_vmodj_left = st_der.state_vmod_left[jshift];
                    int st_vmodj_right = st_der.state_vmod_right[jshift];
                    std::pair <int,int> vmod_left_j  = p_der.IndexToPair[p_der.currentSubspace][st_vmodj_left];
                    std::pair <int,int> vmod_right_j = p_der.IndexToPair[p_der.currentSubspace][st_vmodj_right];
                    int field_i = -1;
                    int field_j = -1;
                    for (int f = 0; f < p_der.n_fields; f++)
                    {
                        if (p_der.fields[f] == std::make_pair(vmod_left_i,vmod_right_i) ){
                            field_i = f+1;
                        }
                        if (p_der.fields[f] == std::make_pair(vmod_left_j,vmod_right_j) ){
                            field_j = f+1;
                        }
                    }
                    if (field_i == -1){
                        std::cout << " Error, field not found! " << field_i << " (" << vmod_left_i.first << "," << vmod_left_i.second << ") (" << vmod_right_i.first << "," << vmod_right_i.second << ")"<< std::endl;
                    }
                    if (field_j == -1){
                        std::cout << " Error, field not found! " << field_j << " (" << vmod_left_j.first << "," << vmod_left_j.second << ") (" << vmod_right_j.first << "," << vmod_right_j.second <<  ")"<< std::endl;
                    }
                    long double SCtmp = p_der.cptsmall(field_i, p_der.r[pert], field_j);

                    double temp = me_der.me_chi_ortho(pert, st_vmod_left, st_vmodj_left, st_left, st_der.state_left[jshift]).first;
                    temp *= me_der.me_chi_ortho(pert, st_vmod_right, st_vmodj_right, st_right, st_der.state_right[jshift]).second;
                    tempD[i*p.Nstep + j] +=  (SCtmp*c*temp); 
                    }
                }

            // diagonal elements (check 2pi factor)
            tempD[i*p.Nstep + i] += 2.*PI/R*(st_der.state_en[shift + i]-p_der.C/12.)/p_der.numPert;
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
 

void minimal_model_Hamiltonian::sweep_block_BC(int iteration, double* coeff, double* coeffSweep, int NRGspacesize, int NRGstoresize)
{
    minimal_model_parameters& p_der = dynamic_cast<minimal_model_parameters&>(p);
    minimal_model_states& st_der = dynamic_cast<minimal_model_states&>(st);
    chirME& me_der = dynamic_cast<chirME&>(me);

    double* H_part = new double[p.Nsize * NRGstoresize ];   
    for (int i=0; i<p.Nsize*NRGstoresize; i++)
       H_part[i] = 0.0;

    for (int pert=0; pert < p_der.numPert; pert++)
        {
        double c = p_der.pConst[pert]*2.*PI*pow(R/2./PI,1.- (p_der.deltas[pert].first + p_der.deltas[pert].second));

        for(int i=0;i < p.Nsize;++i)
            {
            //for(int j=0;j < p.Nstep;++j)
            //    {
                for(int k=0; k < NRGspacesize;++k )
                    {
                    for (int  k1=0; k1 < NRGstoresize; ++k1)
                        {
                        double temp;
                        int st_vmod_left = st_der.state_vmod_left[k];
                        int st_vmod_right = st_der.state_vmod_right[k];
                        int st_left = st_der.state_left[k];
                        int st_right = st_der.state_right[k];
                        int st_vmodk_left = st_der.state_vmod_left[k1];
                        int st_vmodk_right = st_der.state_vmod_right[k1];
                        if( (st_der.cstate_lev(st_vmod_left,st_left) - st_der.cstate_lev(st_vmod_right, st_right)) == \
                            (st_der.cstate_lev(st_vmodk_left,st_der.state_left[k1]) - st_der.cstate_lev(st_vmodk_right, st_der.state_right[k1])))
                            {
                    // Kac_vmod_{right,left}_{i,j} contains the Kac indexes of the chiral and antichiral components, as oppposed to st_vmod_{left,right} which contain the LINEAR index in the SUBSPACE
                    // of the chiral and antichiral component
                            std::pair <int,int> vmod_left_j  = p_der.IndexToPair[p_der.currentSubspace][st_vmod_left];
                            std::pair <int,int> vmod_right_j = p_der.IndexToPair[p_der.currentSubspace][st_vmod_right];
                            std::pair <int,int> vmod_left_k  = p_der.IndexToPair[p_der.currentSubspace][st_vmodk_left];
                            std::pair <int,int> vmod_right_k = p_der.IndexToPair[p_der.currentSubspace][st_vmodk_right];

                            int field_j = -1;
                            int field_k = -1;
                            for (int f = 0; f < p_der.n_fields; f++)
                                {
                                if( p_der.fields[f] == std::make_pair(vmod_left_j, vmod_right_j) ){
                                    field_j = f+1;
                                }
                                if( p_der.fields[f] == std::make_pair(vmod_left_k, vmod_right_k) ){
                                    field_k = f+1;
                                }
                            }
                            if (field_j == -1 || field_k == -1){
                                std::cout << " Error, fields not found! " << field_j << " " << field_k << std::endl;
                            }

                            long double SCtmp = p_der.cptsmall(field_j, p_der.r[pert], field_k);
                            temp = me_der.me_chi_ortho(pert, st_vmod_left, st_vmodk_left, st_left, st_der.state_left[k1]).first;
                            temp *= SCtmp * me_der.me_chi_ortho(pert, st_vmod_right, st_vmodk_right, st_right, st_der.state_right[k1]).second;
                            }
                        else
                            {
                            temp=0.;
                            }
                        H_part[i*NRGstoresize + k1] += c * coeff[i*NRGspacesize + k] *temp;
                        }
                    }
                }
            //}
        }

    double* H_part2 = new double[p.Nsize*p.Nstep ];

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,p.Nsize,p.Nstep,NRGstoresize,1.0,H_part,NRGstoresize,coeffSweep,NRGstoresize,0.0,H_part2,p.Nstep);

    delete [] H_part;
    // basic version till here
 
    //conformal energy contribution
    double* coefftemp= new double[p.Nsize*NRGstoresize];
    for (int i=0;i < p.Nsize;++i)
        {
        for (int k=0; k < NRGstoresize;k++)
            {
            coefftemp[i*NRGstoresize + k]=coeff[i*NRGspacesize+k]*(2.*PI/R*(st_der.state_en[k]-p_der.C/12.));
            }
        }

    double* H_part3=new double[p.Nsize*p.Nstep];
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,p.Nsize,p.Nstep,NRGstoresize,1.0,coefftemp,NRGstoresize,coeffSweep,NRGstoresize,0.0,H_part3,p.Nstep);

    delete [] coefftemp;

    //sum up perturbation and conformal contribution
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



void minimal_model_Hamiltonian::print_eigenValues(int iteration, int sweep, double* eigenValues, std::string fname, int maxState)
{
    minimal_model_states& st_der = dynamic_cast<minimal_model_states&>(st);

    std::ofstream eigen_output;
    if ( Momentum == -1)
        {
        eigen_output.open(fname.c_str(), std::ofstream::out | std::ofstream::app);
        }
    else
        {
        std::stringstream temp;
        temp << fname << "_Mom_" << Momentum;
        eigen_output.open(temp.str().c_str(), std::ofstream::out | std::ofstream::app);
        }

    eigen_output << iteration  << " " << sweep << " " << std::setiosflags(std::ios::fixed) << std::setprecision(14) << st_der.state_eff_en[maxState] << " " << R << " ";

    for (int i=0; i < p.maxEigen; ++i)
        {
        eigen_output  << eigenValues[i] << " ";
        }
    eigen_output << std::endl;
    

}

void minimal_model_Hamiltonian::print_eigenVectors(double* eigenVectors, std::string fname, int size)
{

    if(myrank_mpi == 0)
        {
        std::ofstream eigenVec_output;
        if ( Momentum == -1)
            {
            eigenVec_output.open(fname.c_str(), std::ofstream::out | std::ofstream::app);
            }
        else
            {
            std::stringstream temp;
            temp << fname << "_Mom_" << Momentum; 
            eigenVec_output.open(temp.str().c_str(), std::ofstream::out | std::ofstream::app);
            }

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
