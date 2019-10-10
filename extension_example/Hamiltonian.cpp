/** Hamiltonian class file. 
 * \authors Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date November 2017
 */

#include "Hamiltonian.h"
#include "parameters.h"
#include "tensors.h"
#include "matrix_elements.h"
#include "states.h"
#include "config.h"
#include <cblas.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

Hamiltonian::~Hamiltonian()
{


}


Hamiltonian::Hamiltonian(parameters& p_in, states& st_in, matrix_elements& me_in)
:
p(p_in),
st(st_in),
me(me_in)
{
    if (!p.NRG)
        {
        Msize = st.nstate;
        tensor2D< double> H_tmp(Msize, Msize);
        H = H_tmp;
        }
    else
        {
        Msize = p.Nsize + p.Nstep;
        tensor2D< double> H_tmp(Msize, Msize);
        H = H_tmp;
        }

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

    for (int i=0;i<Msize*Msize;i++)
         {
         H.data[i]=0.;
         }

    for (int pert=0; pert < p.numPert; pert++)
        {
        #pragma omp parallel for 
        for(int i=0; i < Msize; ++i) 
            {
            for(int j=0; j < Msize; ++j)
                    {
                    H(i,j) += R*R * me.me_unperturbed(pert,i,j);
                    }
            // diagonal elements (check 2pi factor)
            H(i,i) += st.state_en[i]/(R*R*R);
    
            }
        }
    //check for hermiticity
    for(int i=0;i<Msize;++i) 
        {
        for(int j=0;j<Msize;++j) 
            {
            if (fabs(H(i,j)-H(j,i)) > .000001)
                std::cout << " Warning - hermiticity problem " << std::endl;
            }
        } 



}

void Hamiltonian::NRG_block_A(double* eigenvalues)
{
    // called first, clears the entire Hamiltonian, then stores the eigenvalues on the diagonal 
    // for the Nsize*Nsize submatrix
    for (int i=0;i<Msize*Msize;i++)
         {
         H.data[i]=0.;
         }
    for( int i=0; i < p.Nsize;++i)
        H(i,i) = eigenvalues[i];

}

void Hamiltonian::NRG_block_BC(int iteration, double* coeff)
{
    int shift = p.Nsize + p.Nstep*iteration;
    double* me_needed = new double[p.Nstep*shift];
    double* results = new double[p.Nstep*p.Nsize];
    for (int pert=0; pert < p.numPert; pert++)
        {
        for(int j=0;j < p.Nstep;++j) 
            {
            int jshift = shift + j;
            int  jshift1 = j*shift;
            for(int k=0; k < shift;++k) 
                {
                me_needed[jshift1+k]= R*R*me.me_unperturbed(pert,jshift,k) ;
                }
            }

        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,p.Nstep,p.Nsize,shift,1.0,me_needed,shift,coeff,shift,0.0,results,p.Nsize);

        for(int i=0;i<p.Nstep;++i)
            {
            int ishift = i + p.Nsize;
            int ishift1 = i*p.Nsize;
            for(int j=0;j< p.Nsize;++j)
                {
                H(ishift, j) = results[ishift1+j];
                H(j,ishift) = results[ishift1+j];
                }
            }
        }

    delete [] results;    
    delete [] me_needed;
}

void Hamiltonian::NRG_block_D(int iteration)
{
   for (int pert=0; pert < p.numPert; pert++)
        {
        //effective strength of each perturbation
        int shift = p.Nsize + p.Nstep*iteration;
        for(int i=0; i < p.Nstep; ++i)
            {
            int ishift = shift + i; 

            //off-diagonal elements
            for(int j=0; j < p.Nstep; ++j)
                {
                int jshift = shift + j;
                H(i + p.Nsize, j + p.Nsize) += R*R * me.me_unperturbed(0,ishift, jshift);
                }

            // diagonal elements (check 2pi factor)
            H(i + p.Nsize, i + p.Nsize ) += st.state_en[shift + i]/(R*R)/p.numPert;

            }
        }


}
 

void Hamiltonian::sweep_block_BC(int iteration, double* coeff, double* coeffSweep, int NRGspacesize, int NRGstoresize)
{


    double* H_part = new double[p.Nsize * NRGstoresize ];   
    for (int i=0; i<p.Nsize*NRGstoresize; i++)
       H_part[i] = 0.0;

    for (int pert=0; pert < p.numPert; pert++)
        {

        for(int i=0;i < p.Nsize;++i)
            {
            //for(int j=0;j < p.Nstep;++j)
            //    {
                for(int k=0; k < NRGspacesize;++k )
                    {
                    for (int  k1=0; k1 < NRGstoresize; ++k1)
                        {
                        H_part[i*NRGstoresize + k1] += R*R * coeff[i*NRGspacesize + k] *me.me_unperturbed(pert,k,k1);
                        }
                    }
                }
            //}
        }

    double* H_part2 = new double[p.Nsize*p.Nstep ];

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,p.Nsize,p.Nstep,NRGstoresize,1.0,H_part,NRGstoresize,coeffSweep,NRGstoresize,0.0,H_part2,p.Nstep);

    delete [] H_part;
    // basic version till here
 
    //unperturbed contribution
    double* coefftemp= new double[p.Nsize*NRGstoresize];
    for (int i=0;i < p.Nsize;++i)
        {
        for (int k=0; k < NRGstoresize;k++)
            {
            coefftemp[i*NRGstoresize + k]=coeff[i*NRGspacesize+k]*st.state_en[k]/(R*R);
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
            H(i,j + p.Nsize) = H_part2[i*p.Nstep +j];
            H(j +p.Nsize,i) = H_part2[i*p.Nstep +j];
            }
        }

    delete [] H_part2;
}



void Hamiltonian::print_eigenValues(int iteration, int sweep, double* eigenValues, std::string fname, int maxState)
{

    std::ofstream eigen_output;
    eigen_output.open(fname.c_str(), std::ofstream::out | std::ofstream::app);

    eigen_output << iteration  << " " << sweep << " " << std::setiosflags(std::ios::fixed) << std::setprecision(14) << R << " ";

    for (int i=0; i < p.maxEigen; ++i)
        {
        eigen_output  << eigenValues[i] << " ";
        }
    eigen_output << std::endl;
    

}

void Hamiltonian::print_eigenVectors(double* eigenVectors, std::string fname, int size)
{
    std::ofstream eigenVec_output;
    eigenVec_output.open(fname.c_str(), std::ofstream::out | std::ofstream::app);

    for (int i=0; i < p.maxEigen; ++i)
        {
        for (int j=0; j < size; ++j)
            {
            eigenVec_output  << eigenVectors[i*size + j]  << " ";
            }
        eigenVec_output << std::endl;
        }

}
