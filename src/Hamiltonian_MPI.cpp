/** Hamiltonian class file. 
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
 */

#include "Hamiltonian_MPI.h"
#include "parameters.h"
#include "tensors.h"
#include "chirME.h"
#include "states.h"
#include "config.h"
#include <mpi.h>
#include "scalapack.h"
#include <cblas.h>
#include <iostream>
#include <fstream>
#include <iomanip>

Hamiltonian::~Hamiltonian()
{
    delete [] Offset;
}


Hamiltonian::Hamiltonian(parameters& p_in, states& st_in, chirME& me_in, int sub)
:
p(p_in),
st(st_in),
me(me_in),
subspace(sub)
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
    
    int sumSizes;
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

void Hamiltonian::init()
{

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

    for (int pert=0; pert < p.numPert; pert++)
        {
        //effective strength of each perturbation
        double c = p.pConst[pert]*2.*pi*pow(R/2./pi,1.-2.*p.deltas[pert]); 

        //for(int i=0; i < Msize; ++i) 
        for(int i=Offset[myrank_mpi]; i < Offset[myrank_mpi+1]; ++i)
            {
            int st_vmod = st.state_vmod[i];
            int st_left = st.state_left[i];
            int st_right  = st.state_right[i];
            //off-diagonal elements
            for(int j=0; j < Msize; ++j)
            //for(int j=localOffsetCol; j < localOffsetCol + colH; ++j) 
                {
                if( ( st.cstate_lev(st.state_vmod[i],st.state_left[i]) - st.cstate_lev(st.state_vmod[i],st.state_right[i]) )==( st.cstate_lev(st.state_vmod[j],st.state_left[j]) - st.cstate_lev(st.state_vmod[j],st.state_right[j]) ) )
                    {
                    int st_vmodj = st.state_vmod[j];
                    double temp = me.me_chi_ortho(pert, st_vmod, st_vmodj, st_left, st.state_left[j]);
                    temp *= me.me_chi_ortho(pert, st_vmod, st_vmodj, st_right, st.state_right[j]);
                    tempH[ (i- Offset[myrank_mpi])*Msize + j] += c*temp;
                    }
                if ( i==j )
                    tempH[ (i-Offset[myrank_mpi])*Msize + i] += 2.*pi/R*(st.state_en[i]-p.C/12.)/p.numPert;
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
        MPI::COMM_WORLD.Bcast(shareH, slice*Msize, MPI_DOUBLE, n); 
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

void Hamiltonian::NRG_block_A(double* eigenvalues)
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

void Hamiltonian::NRG_block_BC(int iteration, double* coeff)
{
    // all ranks calculate block B-C. This could probably be parallelized, even though the speed up could be irrelevant
    int shift = p.Nsize + p.Nstep*iteration;
    double* me_needed = new double[p.Nstep*shift];
    double* results = new double[p.Nstep*p.Nsize];
    for (int pert=0; pert < p.numPert; pert++)
        {
        double c = p.pConst[pert]*2.*pi*pow(R/2./pi,1.-2.*p.deltas[pert]);
        
        for(int j=0;j < p.Nstep;++j) 
            {
            int jshift = shift + j;
            int jshift1 = j*shift;
            for(int k=0; k < shift;++k) 
                {
                double temp;
                int st_vmod = st.state_vmod[jshift];
                int st_left = st.state_left[jshift];
                int st_right = st.state_right[jshift];
                int st_vmodk = st.state_vmod[k];
                if( (st.cstate_lev(st_vmod,st_left) - st.cstate_lev(st_vmod, st_right)) == (st.cstate_lev(st_vmodk,st.state_left[k]) - st.cstate_lev(st_vmodk, st.state_right[k])))
                    {
                    temp = me.me_chi_ortho(pert,st_vmod,st_vmodk,st_left,st.state_left[k]);
                    temp *= me.me_chi_ortho(pert,st_vmod,st_vmodk,st_right,st.state_right[k]);
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

void Hamiltonian::NRG_block_D(int iteration)
{
    double* tempD = new double[p.Nstep*p.Nstep];
    for(int i=0; i < p.Nstep*p.Nstep; ++i)
        tempD[i] = 0.0;

    // all ranks calculate block D
    for (int pert=0; pert < p.numPert; pert++)
        {
        //effective strength of each perturbation
        double c = p.pConst[pert]*2.*pi*pow(R/2./pi,1.-2.*p.deltas[pert]);
        int shift = p.Nsize + p.Nstep*iteration;
        for(int i=0; i < p.Nstep; ++i)
            {
            int ishift = shift + i; 
            int st_vmod = st.state_vmod[ishift];
            int st_left = st.state_left[ishift];
            int st_right  = st.state_right[ishift];
            //off-diagonal elements
            for(int j=0; j < p.Nstep; ++j)
                {
                int jshift = shift + j;
                if( ( st.cstate_lev(st_vmod, st_left) - st.cstate_lev(st_vmod, st_right) )==( st.cstate_lev(st.state_vmod[jshift],st.state_left[jshift]) - st.cstate_lev(st.state_vmod[jshift],st.state_right[jshift]) ) )
                    {
                    int st_vmodj = st.state_vmod[jshift];
                    double temp = me.me_chi_ortho(pert, st_vmod, st_vmodj, st_left, st.state_left[jshift]);
                    temp *= me.me_chi_ortho(pert, st_vmod, st_vmodj, st_right, st.state_right[jshift]);
                    tempD[i*p.Nstep + j] +=  (c*temp); 
                    }
                }

            // diagonal elements (check 2pi factor)
            tempD[i*p.Nstep + i] += 2.*pi/R*(st.state_en[shift + i]-p.C/12.)/p.numPert;
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
 

void Hamiltonian::sweep_block_BC(int iteration, double* coeff, double* coeffSweep, int NRGspacesize, int NRGstoresize)
{

    int subsizeMax = 0;
    for ( int i = 0; i<p.subdim[subspace]; ++i)
        {
        for (int j = 0; j<=p.maxLev; ++j)
            {
            if ( st.nstates_lev_each( i, j) > subsizeMax)
                subsizeMax = st.nstates_lev_each(i, j);
            }
        }

    tensor6D<double> m1(p.numPert, p.subdim[subspace], p.subdim[subspace], p.maxLev+1, p.maxLev+1, subsizeMax*subsizeMax);
    tensor6D<double> m2(p.numPert, p.subdim[subspace], p.subdim[subspace], p.maxLev+1, p.maxLev+1, subsizeMax*subsizeMax);

    //storing relevant chiral matrix elements
    for ( int pert=0; pert<p.numPert; ++pert)
        {
        for ( int vmod2=0; vmod2 < p.subdim[subspace]; ++vmod2)
            {
            for ( int lev2=p.minLev; lev2 <=p.maxLev; ++lev2)
                {
                 for ( int vmod=0; vmod < p.subdim[subspace]; ++vmod)
                    {
                    for ( int lev=p.minLev; lev <=p.maxLev; ++lev)
                        {
                        int subsize = st.nstates_lev_each(vmod,lev);
                        int subsize2 = st.nstates_lev_each(vmod2,lev2);
                        if (subsize!=0 and subsize2!=0)
                            {
                            for (int ii=0;ii<subsize2;ii++)
                                {
                                for (int jj=0;jj<subsize;jj++)
                                    {
                                    int ind1 = st.state_left_lev(vmod2,lev2,ii);
                                    int ind2 = st.state_left_lev(vmod,lev,jj);
                                    m1(pert, vmod2, vmod, lev2, lev, ii*subsize+jj)=me.me_chi_ortho(pert, vmod2, vmod, ind1, ind2);
                                        
                                    ind1 = st.state_right_lev(vmod2, lev2, ii);
                                    ind2 = st.state_right_lev(vmod, lev, jj);
                                    m2(pert, vmod2, vmod, lev2, lev, ii*subsize+jj) = me.me_chi_ortho(pert, vmod2, vmod, ind1, ind2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
   

    double* coeffSweep_j;
    double* results;
    double* results2;
    double* H_part= new double[p.Nstep*NRGspacesize];

    for (int i=0; i<p.Nstep*NRGspacesize; ++i)
        H_part[i]=0.;        

    for (int state=0;state < p.Nstep;state++)
        {
        //double temp=0.; 
        for ( int vmod2=0; vmod2 < p.subdim[subspace]; ++vmod2)
            {
            for (int lev2 = p.minLev; lev2 <= p.maxLev; ++lev2)
                {
                int subsize2 = st.nstates_lev_each(vmod2,lev2);
                for ( int vmod=0; vmod < p.subdim[subspace]; ++vmod)
                    {
                    for(int lev = p.minLev ; lev<= p.maxLev; ++lev)
                        {
                        int subsize = st.nstates_lev_each(vmod,lev);
                        if (subsize!=0 and subsize2!=0)
                            {
                            coeffSweep_j=new double[subsize*subsize];
                            for (int ii=0;ii<subsize;ii++)
                                {
                                for (int jj=0;jj<subsize;jj++)
                                    {
                                    if (st.state_nonchiral_lev(vmod, lev, ii, jj) < NRGstoresize)
                                        {
                                        coeffSweep_j[ii*subsize+jj] = coeffSweep[ state*NRGstoresize + st.state_nonchiral_lev(vmod,lev,ii,jj) ];
                                        }
                                    else
                                        {
                                        coeffSweep_j[ii*subsize+jj]=0.;
                                        }
                                    }
                                }


                            for(int pert=0; pert < p.numPert; ++pert)
                                {
                                results=new double[subsize*subsize2];
                                results2=new double[subsize2*subsize2];

                                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,subsize2,subsize,subsize,1.0,&m1(pert,vmod2,vmod,lev2,lev,0),subsize,coeffSweep_j,subsize,0.0,results,subsize);
                                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,subsize2,subsize2,subsize,1.0,results,subsize,&m2(pert,vmod2,vmod,lev2,lev,0),subsize,0.0,results2,subsize2);

                                double c = p.pConst[pert]*2.*pi*pow(R/2./pi,1.-2.*p.deltas[pert]);
                                for (int ii=0;ii<subsize2;ii++)
                                    {
                                    for (int jj=0;jj<subsize2;jj++)
                                        {
                                        int ind1 = st.state_nonchiral_lev(vmod2, lev2, ii, jj);
                                        if (ind1 < NRGspacesize)
                                            {
                                            H_part[state*NRGspacesize+ind1] += c*results2[ii*subsize2+jj];
                                            }
                                        }
                                    }
                                delete [] results;
                                delete [] results2;
                                }//closing loop over perturbation
                            delete [] coeffSweep_j;
                            }//closing if
                        } //closing loop over level -->lev
                    }// closing loop over verma modules -->vmod
                } // closing loop over lev -->lev2
            } // closing loop over verma modules --> vmod2
        } //closing loop over state 
     

    double* H_part2 = new double[p.Nsize*p.Nstep ];

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,p.Nsize,p.Nstep,NRGspacesize,1.0,coeff,NRGspacesize,H_part,NRGspacesize,0.0,H_part2,p.Nstep);    

    delete [] H_part;

    double* coefftemp= new double[p.Nsize*NRGstoresize];
    for (int i=0;i < p.Nsize;++i)
        {
        for (int k=0; k < NRGstoresize;k++)
            {
            coefftemp[i*NRGstoresize + k]=coeff[i*NRGspacesize+k]*(2.*pi/R*(st.state_en[k]-p.C/12.));
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



void Hamiltonian::print_eigenValues(int iteration, int sweep, double* eigenValues, std::string fname, int maxState)
{

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

    eigen_output << iteration  << " " << sweep << " " << std::setiosflags(std::ios::fixed) << std::setprecision(14) << st.state_eff_en[maxState] << " " << R << " ";

    for (int i=0; i < p.maxEigen; ++i)
        {
        eigen_output  << eigenValues[i] << " ";
        }
    eigen_output << std::endl;
    

}

void Hamiltonian::print_eigenVectors(double* eigenVectors, std::string fname, int size)
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
